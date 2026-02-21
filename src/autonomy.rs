//! Autonomous spacecraft control — trajectory models and correction decisions.

use crate::orbit::SpacecraftState;
use crate::comm::ModelDifferential;

/// Autonomy level classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AutonomyLevel {
    GroundControlled = 0,
    Supervised = 1,
    Conditional = 2,
    HighAutonomy = 3,
    FullAutonomy = 4,
}

/// Polynomial trajectory model: position as f(t).
#[derive(Debug, Clone)]
pub struct TrajectoryModel {
    /// Coefficients in groups of 3 (x, y, z): \[x0, y0, z0, vx, vy, vz, ax, ay, az, ...\]
    pub coefficients: Vec<f64>,
    pub epoch_ns: u64,
    pub valid_duration_ns: u64,
}

impl TrajectoryModel {
    pub fn new(coefficients: Vec<f64>, epoch_ns: u64, valid_ns: u64) -> Self {
        Self { coefficients, epoch_ns, valid_duration_ns: valid_ns }
    }

    /// Evaluate position at time t_ns.
    /// x(t) = c\[0\] + c\[3\]\*dt + c\[6\]\*dt², etc.
    pub fn evaluate(&self, t_ns: u64) -> [f64; 3] {
        let dt = (t_ns.saturating_sub(self.epoch_ns)) as f64 / 1e9;
        let mut pos = [0.0; 3];

        // Constant terms (c[0], c[1], c[2])
        for (i, p) in pos.iter_mut().enumerate() {
            if i < self.coefficients.len() {
                *p = self.coefficients[i];
            }
        }
        // Linear terms (c[3], c[4], c[5])
        for (i, p) in pos.iter_mut().enumerate() {
            if 3 + i < self.coefficients.len() {
                *p += self.coefficients[3 + i] * dt;
            }
        }
        // Quadratic terms (c[6], c[7], c[8])
        let dt2 = dt * dt;
        for (i, p) in pos.iter_mut().enumerate() {
            if 6 + i < self.coefficients.len() {
                *p += self.coefficients[6 + i] * dt2;
            }
        }
        pos
    }

    /// Check if model is valid at given time.
    #[inline]
    pub fn is_valid_at(&self, t_ns: u64) -> bool {
        t_ns >= self.epoch_ns && t_ns <= self.epoch_ns + self.valid_duration_ns
    }
}

/// Apply a model differential to update trajectory coefficients.
pub fn apply_differential(model: &mut TrajectoryModel, diff: &ModelDifferential) {
    for &(param_hash, value) in &diff.param_updates {
        // Use param_hash as coefficient index (low bits)
        let idx = (param_hash & 0xFF) as usize;
        if idx < model.coefficients.len() {
            model.coefficients[idx] = value;
        }
    }
}

/// A control decision (thrust correction).
#[derive(Debug, Clone)]
pub struct ControlDecision {
    pub action_hash: u64,
    pub thrust_vector: [f64; 3],
    pub burn_duration_s: f64,
    pub confidence: f64,
    pub timestamp_ns: u64,
}

/// Compute correction to align current state with model prediction.
pub fn compute_correction(current: &SpacecraftState, model: &TrajectoryModel, t_ns: u64) -> ControlDecision {
    let predicted = model.evaluate(t_ns);
    let error = [
        predicted[0] - current.position_km[0],
        predicted[1] - current.position_km[1],
        predicted[2] - current.position_km[2],
    ];
    let error_magnitude = (error[0] * error[0] + error[1] * error[1] + error[2] * error[2]).sqrt();

    // Thrust proportional to error, normalized — precompute reciprocal to avoid 3 divisions
    let thrust = if error_magnitude > 1e-10 {
        let rcp = 1.0 / error_magnitude;
        [error[0] * rcp, error[1] * rcp, error[2] * rcp]
    } else {
        [0.0, 0.0, 0.0]
    };

    // Confidence inversely proportional to error
    let confidence = (1.0 / (1.0 + error_magnitude * 0.01)).clamp(0.0, 1.0);

    // Burn duration proportional to error
    let burn_duration = (error_magnitude * 0.001).min(300.0);

    ControlDecision {
        action_hash: 0,
        thrust_vector: thrust,
        burn_duration_s: burn_duration,
        confidence,
        timestamp_ns: t_ns,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orbit::SpacecraftState;

    #[test]
    fn trajectory_evaluate_constant() {
        let model = TrajectoryModel::new(vec![100.0, 200.0, 300.0], 0, 1_000_000_000);
        let pos = model.evaluate(500_000_000);
        assert!((pos[0] - 100.0).abs() < 1e-10);
        assert!((pos[1] - 200.0).abs() < 1e-10);
        assert!((pos[2] - 300.0).abs() < 1e-10);
    }

    #[test]
    fn trajectory_evaluate_linear() {
        // x(t) = 0 + 10*t (vx=10 km/s)
        let model = TrajectoryModel::new(vec![0.0, 0.0, 0.0, 10.0, 0.0, 0.0], 0, 10_000_000_000);
        let pos = model.evaluate(1_000_000_000); // t = 1 second
        assert!((pos[0] - 10.0).abs() < 1e-6);
    }

    #[test]
    fn trajectory_validity() {
        let model = TrajectoryModel::new(vec![0.0; 3], 1000, 5000);
        assert!(model.is_valid_at(1000));
        assert!(model.is_valid_at(6000));
        assert!(!model.is_valid_at(999));
        assert!(!model.is_valid_at(6001));
    }

    #[test]
    fn correction_zero_error() {
        let state = SpacecraftState {
            position_km: [100.0, 200.0, 300.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            timestamp_ns: 0,
            fuel_kg: 100.0,
        };
        let model = TrajectoryModel::new(vec![100.0, 200.0, 300.0], 0, 1_000_000_000);
        let decision = compute_correction(&state, &model, 0);
        assert!(decision.burn_duration_s < 1e-6);
        assert!(decision.confidence > 0.99);
    }

    #[test]
    fn correction_with_error() {
        let state = SpacecraftState {
            position_km: [0.0, 0.0, 0.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            timestamp_ns: 0,
            fuel_kg: 100.0,
        };
        let model = TrajectoryModel::new(vec![1000.0, 0.0, 0.0], 0, 1_000_000_000);
        let decision = compute_correction(&state, &model, 0);
        assert!(decision.thrust_vector[0] > 0.0); // should thrust toward target
        assert!(decision.burn_duration_s > 0.0);
    }

    #[test]
    fn confidence_bounds() {
        let state = SpacecraftState {
            position_km: [0.0, 0.0, 0.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            timestamp_ns: 0,
            fuel_kg: 50.0,
        };
        let model = TrajectoryModel::new(vec![99999.0, 0.0, 0.0], 0, 1_000_000_000);
        let d = compute_correction(&state, &model, 0);
        assert!(d.confidence >= 0.0 && d.confidence <= 1.0);
    }
}
