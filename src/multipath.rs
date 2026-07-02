//! Multipath detection via Code-Minus-Carrier (`CMC`).
//!
//! Multipath is the arrival of the same satellite signal via reflected
//! paths in addition to the direct line-of-sight ray. It biases the code
//! pseudorange by tens of centimetres to several metres, while the
//! carrier phase is affected by at most a few centimetres. Differencing
//! code and carrier therefore isolates the multipath contribution:
//!
//! ```text
//! CMC = P - Φ - 2·I - offsets
//! ```
//!
//! where `P` is the pseudorange, `Φ` is the wavelength-scaled carrier
//! phase, and `I` is the ionospheric delay (approximated by half of the
//! `GF` LC on a dual-frequency receiver). The residual variability of
//! `CMC` over a short window (typically 20–60 s) is a direct multipath
//! indicator.
//!
//! # References
//!
//! - Estey, L. H. & Meertens, C. M. (1999), "TEQC: The multi-purpose
//!   toolkit for GPS/GLONASS data", GPS Solutions, 3(1), 42-49.
//! - Kaplan, E. D. & Hegarty, C. J. (2017), §7.3 Multipath and Errors.

use crate::cycle_slip::{L1_HZ, L2_HZ};

const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

fn wavelength_m(freq_hz: f64) -> f64 {
    SPEED_OF_LIGHT_M_S / freq_hz
}

// ---------------------------------------------------------------------------
// CmcSample
// ---------------------------------------------------------------------------

/// One dual-frequency observation used to compute a CMC sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CmcInput {
    /// Reception epoch in seconds.
    pub epoch_s: f64,
    /// L1 carrier phase in cycles.
    pub l1_cycles: f64,
    /// L2 carrier phase in cycles.
    pub l2_cycles: f64,
    /// L1 pseudorange in metres.
    pub p1_m: f64,
}

impl CmcInput {
    /// Compute the L1 code-minus-carrier value in metres.
    ///
    /// `CMC = P1 − λ1·L1 − 2·(λ1·L1 − λ2·L2) · f2² / (f1² − f2²)`
    ///
    /// The last term removes the ionospheric contribution using the
    /// geometry-free LC.
    #[must_use]
    pub fn cmc_l1_m(&self) -> f64 {
        let lambda1 = wavelength_m(L1_HZ);
        let lambda2 = wavelength_m(L2_HZ);
        let f1sq = L1_HZ * L1_HZ;
        let f2sq = L2_HZ * L2_HZ;
        let gf = lambda1 * self.l1_cycles - lambda2 * self.l2_cycles;
        let iono_l1 = 2.0 * gf * f2sq / (f1sq - f2sq);
        self.p1_m - lambda1 * self.l1_cycles - iono_l1
    }
}

// ---------------------------------------------------------------------------
// Detector
// ---------------------------------------------------------------------------

/// Sliding-window multipath detector.
///
/// The detector maintains up to `window_size` CMC samples and reports
/// the running standard deviation. When the std exceeds
/// [`Detector::threshold_m`], the epoch is flagged as multipath-degraded.
#[derive(Debug, Clone)]
pub struct Detector {
    /// Maximum number of samples retained in the window.
    pub window_size: usize,
    /// Std threshold in metres. Typical benign environments give < 0.1 m.
    pub threshold_m: f64,
    samples: Vec<f64>,
}

impl Detector {
    /// Construct a new detector with the given window and threshold.
    #[must_use]
    pub fn new(window_size: usize, threshold_m: f64) -> Self {
        Self {
            window_size,
            threshold_m,
            samples: Vec::with_capacity(window_size),
        }
    }

    /// Push a new CMC sample and return a report for the current window.
    pub fn push(&mut self, sample_m: f64) -> Report {
        if self.samples.len() == self.window_size && !self.samples.is_empty() {
            self.samples.remove(0);
        }
        self.samples.push(sample_m);
        let mean = self.samples.iter().sum::<f64>() / (self.samples.len() as f64);
        let var = if self.samples.len() > 1 {
            self.samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                / ((self.samples.len() - 1) as f64)
        } else {
            0.0
        };
        let std = var.sqrt();
        Report {
            window_len: self.samples.len(),
            mean_m: mean,
            std_m: std,
            multipath_flag: std > self.threshold_m,
        }
    }

    /// Reset the detector (typically called after a satellite reacquisition).
    pub fn reset(&mut self) {
        self.samples.clear();
    }

    /// Current window length.
    #[must_use]
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Whether the window is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

/// Report emitted by [`Detector::push`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Report {
    /// Number of samples in the running window.
    pub window_len: usize,
    /// Sample mean of the window in metres.
    pub mean_m: f64,
    /// Sample standard deviation of the window in metres.
    pub std_m: f64,
    /// Whether the std exceeds the configured threshold.
    pub multipath_flag: bool,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn quiet_sample(epoch_s: f64) -> CmcInput {
        CmcInput {
            epoch_s,
            l1_cycles: 100_000_000.0 + epoch_s * 0.1,
            l2_cycles: 77_000_000.0 + epoch_s * 0.1 * L2_HZ / L1_HZ,
            p1_m: 20_200_000.0 + epoch_s * wavelength_m(L1_HZ) * 0.1,
        }
    }

    #[test]
    fn cmc_is_finite_for_quiet_input() {
        assert!(quiet_sample(0.0).cmc_l1_m().is_finite());
    }

    #[test]
    fn quiet_stream_stays_below_threshold() {
        let mut det = Detector::new(20, 0.5);
        for i in 0..30 {
            let s = quiet_sample(f64::from(i)).cmc_l1_m();
            let r = det.push(s);
            if r.window_len >= 5 {
                assert!(!r.multipath_flag, "false alarm at i={i}, std={}", r.std_m);
            }
        }
    }

    #[test]
    fn noisy_stream_triggers_flag() {
        let mut det = Detector::new(10, 0.1);
        for i in 0..10 {
            let noise = if i % 2 == 0 { 1.0 } else { -1.0 };
            let s = quiet_sample(f64::from(i)).cmc_l1_m() + noise;
            det.push(s);
        }
        // After 10 alternating ±1 m samples, std should exceed 0.1 m.
        let last = det.push(quiet_sample(11.0).cmc_l1_m() + 1.0);
        assert!(last.std_m > 0.1);
        assert!(last.multipath_flag);
    }

    #[test]
    fn window_size_is_respected() {
        let mut det = Detector::new(3, 100.0);
        for i in 0..10 {
            det.push(f64::from(i));
        }
        assert_eq!(det.len(), 3);
    }

    #[test]
    fn empty_detector_reports_zero_std() {
        let mut det = Detector::new(10, 0.5);
        let r = det.push(1.23);
        assert_eq!(r.window_len, 1);
        assert!((r.std_m - 0.0).abs() < 1e-12);
    }

    #[test]
    fn reset_clears_samples() {
        let mut det = Detector::new(5, 0.5);
        det.push(1.0);
        det.push(2.0);
        det.reset();
        assert!(det.is_empty());
    }

    #[test]
    fn mean_matches_analytic_value() {
        let mut det = Detector::new(3, 100.0);
        let r = det.push(1.0);
        assert!((r.mean_m - 1.0).abs() < 1e-12);
        det.push(2.0);
        let r = det.push(3.0);
        assert!((r.mean_m - 2.0).abs() < 1e-12);
    }

    #[test]
    fn threshold_boundary() {
        let mut det = Detector::new(2, 0.5);
        det.push(0.0);
        let r = det.push(0.0);
        assert!(!r.multipath_flag);
    }
}
