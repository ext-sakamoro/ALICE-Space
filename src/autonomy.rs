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

/// Fault type classification for autonomous fault detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaultType {
    /// No fault detected.
    Nominal,
    /// Sensor reading out of expected range.
    SensorAnomaly,
    /// Thruster underperformance or failure.
    ThrusterDegradation,
    /// Communication link degradation.
    CommLinkLoss,
    /// Power subsystem fault (solar panel, battery).
    PowerAnomaly,
    /// Attitude control divergence.
    AttitudeError,
}

/// A node in a decision tree for autonomous fault response.
#[derive(Debug, Clone)]
pub struct DecisionNode {
    /// Threshold for the evaluated metric.
    pub threshold: f64,
    /// Fault type to assign if metric exceeds threshold.
    pub fault: FaultType,
    /// Severity weight (0.0 = informational, 1.0 = critical).
    pub severity: f64,
}

/// Evaluate a decision tree against a vector of sensor readings.
///
/// Each `DecisionNode` checks whether the corresponding reading exceeds its
/// threshold. Returns the highest-severity fault detected, or `FaultType::Nominal`.
pub fn evaluate_decision_tree(readings: &[f64], tree: &[DecisionNode]) -> (FaultType, f64) {
    let mut worst_fault = FaultType::Nominal;
    let mut worst_severity = 0.0_f64;

    for (reading, node) in readings.iter().zip(tree.iter()) {
        let abs_reading = reading.abs();
        if abs_reading > node.threshold {
            if node.severity > worst_severity {
                worst_severity = node.severity;
                worst_fault = node.fault;
            }
        }
    }
    (worst_fault, worst_severity)
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
    fn fault_type_nominal_default() {
        let readings = [0.1, 0.2, 0.05];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::SensorAnomaly, severity: 0.5 },
            DecisionNode { threshold: 1.0, fault: FaultType::ThrusterDegradation, severity: 0.7 },
            DecisionNode { threshold: 1.0, fault: FaultType::PowerAnomaly, severity: 0.9 },
        ];
        let (fault, sev) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::Nominal);
        assert!(sev < 1e-10);
    }

    #[test]
    fn fault_type_single_trigger() {
        let readings = [0.1, 5.0, 0.05];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::SensorAnomaly, severity: 0.5 },
            DecisionNode { threshold: 1.0, fault: FaultType::ThrusterDegradation, severity: 0.7 },
            DecisionNode { threshold: 1.0, fault: FaultType::PowerAnomaly, severity: 0.9 },
        ];
        let (fault, sev) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::ThrusterDegradation);
        assert!((sev - 0.7).abs() < 1e-10);
    }

    #[test]
    fn fault_type_highest_severity_wins() {
        let readings = [5.0, 5.0, 5.0];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::SensorAnomaly, severity: 0.5 },
            DecisionNode { threshold: 1.0, fault: FaultType::ThrusterDegradation, severity: 0.7 },
            DecisionNode { threshold: 1.0, fault: FaultType::PowerAnomaly, severity: 0.9 },
        ];
        let (fault, sev) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::PowerAnomaly);
        assert!((sev - 0.9).abs() < 1e-10);
    }

    #[test]
    fn fault_type_negative_reading() {
        let readings = [-5.0];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::AttitudeError, severity: 0.8 },
        ];
        let (fault, _) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::AttitudeError);
    }

    #[test]
    fn fault_type_empty_inputs() {
        let (fault, sev) = evaluate_decision_tree(&[], &[]);
        assert_eq!(fault, FaultType::Nominal);
        assert!(sev < 1e-10);
    }

    #[test]
    fn fault_type_mismatched_lengths() {
        let readings = [5.0, 5.0, 5.0, 5.0, 5.0];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::CommLinkLoss, severity: 0.6 },
        ];
        let (fault, sev) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::CommLinkLoss);
        assert!((sev - 0.6).abs() < 1e-10);
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

    #[test]
    fn trajectory_evaluate_quadratic() {
        // x(t) = 0 + 0*t + 5*t^2, y=z=0
        let model = TrajectoryModel::new(
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0],
            0,
            10_000_000_000,
        );
        // At t = 2s: x = 5 * 4 = 20
        let pos = model.evaluate(2_000_000_000);
        assert!((pos[0] - 20.0).abs() < 1e-6, "quadratic pos[0]={}", pos[0]);
        assert!(pos[1].abs() < 1e-10);
        assert!(pos[2].abs() < 1e-10);
    }

    #[test]
    fn trajectory_before_epoch_clamps_to_zero_dt() {
        let model = TrajectoryModel::new(
            vec![100.0, 200.0, 300.0, 10.0, 20.0, 30.0],
            1_000_000_000,
            5_000_000_000,
        );
        // Before epoch: saturating_sub yields 0 → dt = 0 → constant terms only
        let pos = model.evaluate(0);
        assert!((pos[0] - 100.0).abs() < 1e-10);
        assert!((pos[1] - 200.0).abs() < 1e-10);
        assert!((pos[2] - 300.0).abs() < 1e-10);
    }

    #[test]
    fn trajectory_empty_coefficients() {
        let model = TrajectoryModel::new(vec![], 0, 1_000_000_000);
        let pos = model.evaluate(500_000_000);
        assert_eq!(pos, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn trajectory_partial_coefficients() {
        // Only x constant, no y or z
        let model = TrajectoryModel::new(vec![42.0], 0, 1_000_000_000);
        let pos = model.evaluate(500_000_000);
        assert!((pos[0] - 42.0).abs() < 1e-10);
        assert!(pos[1].abs() < 1e-10);
        assert!(pos[2].abs() < 1e-10);
    }

    #[test]
    fn apply_differential_updates_coefficients() {
        let mut model = TrajectoryModel::new(vec![1.0, 2.0, 3.0, 4.0, 5.0], 0, 1_000_000_000);
        let mut diff = crate::comm::ModelDifferential::new(1, 100);
        // FNV-1a hash of "x" modulo 256 gives index into coefficients
        let _hash = crate::fnv1a("x".as_bytes());
        // Manually push with a known index that fits
        diff.param_updates.push((0, 99.0)); // idx=0
        apply_differential(&mut model, &diff);
        assert!((model.coefficients[0] - 99.0).abs() < 1e-10);
    }

    #[test]
    fn apply_differential_out_of_range_ignored() {
        let mut model = TrajectoryModel::new(vec![1.0, 2.0, 3.0], 0, 1_000_000_000);
        let mut diff = crate::comm::ModelDifferential::new(1, 100);
        // Index 255 is out of range for a 3-element vec → should be ignored
        diff.param_updates.push((255, 99.0));
        apply_differential(&mut model, &diff);
        assert!((model.coefficients[0] - 1.0).abs() < 1e-10);
        assert!((model.coefficients[1] - 2.0).abs() < 1e-10);
        assert!((model.coefficients[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn correction_thrust_unit_vector() {
        let state = SpacecraftState {
            position_km: [0.0, 0.0, 0.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            timestamp_ns: 0,
            fuel_kg: 100.0,
        };
        let model = TrajectoryModel::new(vec![300.0, 400.0, 0.0], 0, 1_000_000_000);
        let d = compute_correction(&state, &model, 0);
        let mag = (d.thrust_vector[0].powi(2) + d.thrust_vector[1].powi(2) + d.thrust_vector[2].powi(2)).sqrt();
        assert!((mag - 1.0).abs() < 1e-10, "Thrust should be unit vector, mag={}", mag);
    }

    #[test]
    fn correction_burn_duration_capped_at_300() {
        let state = SpacecraftState {
            position_km: [0.0, 0.0, 0.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            timestamp_ns: 0,
            fuel_kg: 100.0,
        };
        // Very large error → burn should cap at 300s
        let model = TrajectoryModel::new(vec![1e9, 0.0, 0.0], 0, 1_000_000_000);
        let d = compute_correction(&state, &model, 0);
        assert!((d.burn_duration_s - 300.0).abs() < 1e-10, "Burn capped at 300, got {}", d.burn_duration_s);
    }

    #[test]
    fn fault_type_exact_threshold_no_trigger() {
        // Reading exactly at threshold should NOT trigger (> not >=)
        let readings = [1.0];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::SensorAnomaly, severity: 0.5 },
        ];
        let (fault, _) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::Nominal);
    }

    #[test]
    fn fault_type_just_above_threshold() {
        let readings = [1.0001];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::SensorAnomaly, severity: 0.5 },
        ];
        let (fault, sev) = evaluate_decision_tree(&readings, &tree);
        assert_eq!(fault, FaultType::SensorAnomaly);
        assert!((sev - 0.5).abs() < 1e-10);
    }

    #[test]
    fn autonomy_level_ordering() {
        assert_eq!(AutonomyLevel::GroundControlled as u8, 0);
        assert_eq!(AutonomyLevel::Supervised as u8, 1);
        assert_eq!(AutonomyLevel::Conditional as u8, 2);
        assert_eq!(AutonomyLevel::HighAutonomy as u8, 3);
        assert_eq!(AutonomyLevel::FullAutonomy as u8, 4);
    }

    #[test]
    fn fault_type_equal_severity_first_wins() {
        // Two faults with identical severity: the first one evaluated wins
        let readings = [5.0, 5.0];
        let tree = vec![
            DecisionNode { threshold: 1.0, fault: FaultType::SensorAnomaly, severity: 0.8 },
            DecisionNode { threshold: 1.0, fault: FaultType::CommLinkLoss, severity: 0.8 },
        ];
        let (fault, _) = evaluate_decision_tree(&readings, &tree);
        // Second cannot surpass the first with > comparison, so first wins
        assert_eq!(fault, FaultType::SensorAnomaly);
    }
}
