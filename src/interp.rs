//! Interpolation of tabulated trajectories: Lagrange (positions) and Hermite
//! (positions with velocities), the two methods a CCSDS OEM names in its
//! `INTERPOLATION` keyword, also used for SP3 precise orbits.
//!
//! Both work on a window of samples around the requested time: `degree + 1`
//! samples for Lagrange, `(degree + 1) / 2` samples (rounded up) for Hermite, each
//! carrying a position and a velocity. The window is the one whose nodes are
//! centred on the time as closely as the table allows (it slides inward at the
//! ends of the table), which is how CCSDS 502.0 describes the use of the degree.
//! Times are taken relative to the window's first node before any product is
//! formed, so a table far from the origin of its time scale keeps its precision.
//!
//! The formulas use `+ − × ÷` only, so a result is the same on every target.
//!
//! Author: Moroya Sakamoto

/// Why a value could not be interpolated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InterpError {
    /// The table has fewer samples than the method and degree need
    /// (needed, available).
    NotEnoughSamples(usize, usize),
    /// The time is outside the span of the table (the time).
    OutOfRange(f64),
    /// The sample times are not strictly increasing (the index of the first
    /// sample that is not after its predecessor), or one is not finite.
    UnorderedTimes(usize),
    /// The times and the values have different lengths.
    LengthMismatch,
}

impl core::fmt::Display for InterpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NotEnoughSamples(n, have) => {
                write!(f, "interpolation needs {n} samples, the table has {have}")
            }
            Self::OutOfRange(t) => write!(f, "time {t} is outside the table"),
            Self::UnorderedTimes(i) => write!(f, "sample time {i} is not after the previous one"),
            Self::LengthMismatch => write!(f, "times and values differ in length"),
        }
    }
}

impl std::error::Error for InterpError {}

/// A position and a velocity (any consistent units).
pub type State = ([f64; 3], [f64; 3]);

fn check_times(times: &[f64]) -> Result<(), InterpError> {
    for (i, t) in times.iter().enumerate() {
        if !t.is_finite() || (i > 0 && *t <= times[i - 1]) {
            return Err(InterpError::UnorderedTimes(i));
        }
    }
    Ok(())
}

/// The first index of the `count` nodes centred on `t` within `times`.
fn window_start(times: &[f64], t: f64, count: usize) -> usize {
    // index of the last node at or before t
    let k = times.partition_point(|&x| x <= t).saturating_sub(1);
    let half = (count - 1) / 2;
    k.saturating_sub(half).min(times.len() - count)
}

/// Lagrange interpolation of degree `degree` of the 3-vectors `values` given at
/// `times`, at `t`: the value and its time derivative (the derivative of the
/// interpolating polynomial).
///
/// # Errors
///
/// [`InterpError`] when the table is too short for the degree, unordered, of
/// mismatched length, or `t` is outside it.
pub fn lagrange(
    times: &[f64],
    values: &[[f64; 3]],
    degree: usize,
    t: f64,
) -> Result<State, InterpError> {
    if times.len() != values.len() {
        return Err(InterpError::LengthMismatch);
    }
    let count = degree + 1;
    if times.len() < count {
        return Err(InterpError::NotEnoughSamples(count, times.len()));
    }
    check_times(times)?;
    if !(times[0]..=times[times.len() - 1]).contains(&t) {
        return Err(InterpError::OutOfRange(t));
    }
    let s = window_start(times, t, count);
    let x: Vec<f64> = times[s..s + count]
        .iter()
        .map(|&ti| ti - times[s])
        .collect();
    let y = &values[s..s + count];
    let u = t - times[s];

    let mut p = [0.0; 3];
    let mut dp = [0.0; 3];
    for j in 0..count {
        // basis L_j(u) = prod_{m != j} (u - x_m) / (x_j - x_m) and its derivative
        let mut denom = 1.0;
        for m in 0..count {
            if m != j {
                denom *= x[j] - x[m];
            }
        }
        let mut num = 1.0;
        let mut dnum = 0.0;
        for m in 0..count {
            if m != j {
                dnum = dnum * (u - x[m]) + num;
                num *= u - x[m];
            }
        }
        let (l, dl) = (num / denom, dnum / denom);
        for k in 0..3 {
            p[k] += l * y[j][k];
            dp[k] += dl * y[j][k];
        }
    }
    Ok((p, dp))
}

/// Hermite interpolation of the states `(positions, velocities)` given at
/// `times`, at `t`, matching position and velocity at each node of the window:
/// `degree` is the degree of the interpolating polynomial, so the window has
/// `(degree + 1) / 2` nodes (rounded up). Returns the position and the velocity.
///
/// # Errors
///
/// [`InterpError`] as for [`lagrange`].
pub fn hermite(
    times: &[f64],
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    degree: usize,
    t: f64,
) -> Result<State, InterpError> {
    if times.len() != positions.len() || times.len() != velocities.len() {
        return Err(InterpError::LengthMismatch);
    }
    let count = (degree + 2) / 2;
    if times.len() < count {
        return Err(InterpError::NotEnoughSamples(count, times.len()));
    }
    check_times(times)?;
    if !(times[0]..=times[times.len() - 1]).contains(&t) {
        return Err(InterpError::OutOfRange(t));
    }
    let s = window_start(times, t, count);
    let u = t - times[s];
    // doubled nodes z_{2i} = z_{2i+1} = x_i, Newton divided differences
    let n2 = 2 * count;
    let z: Vec<f64> = (0..n2).map(|i| times[s + i / 2] - times[s]).collect();
    let mut p = [0.0; 3];
    let mut dp = [0.0; 3];
    for k in 0..3 {
        let mut q: Vec<f64> = (0..n2).map(|i| positions[s + i / 2][k]).collect();
        // first differences: velocity at a doubled node, slope otherwise
        let mut col: Vec<f64> = (1..n2)
            .map(|i| {
                if i % 2 == 1 {
                    velocities[s + i / 2][k]
                } else {
                    (q[i] - q[i - 1]) / (z[i] - z[i - 1])
                }
            })
            .collect();
        let mut coef = vec![q[0], col[0]];
        for level in 2..n2 {
            let next: Vec<f64> = (0..col.len() - 1)
                .map(|i| (col[i + 1] - col[i]) / (z[i + level] - z[i]))
                .collect();
            coef.extend([next[0]]);
            col = next;
        }
        q.clear();
        // Horner on the Newton form, with the derivative
        let mut v = coef[n2 - 1];
        let mut dv = 0.0;
        for i in (0..n2 - 1).rev() {
            dv = dv * (u - z[i]) + v;
            v = v * (u - z[i]) + coef[i];
        }
        p[k] = v;
        dp[k] = dv;
    }
    Ok((p, dp))
}
