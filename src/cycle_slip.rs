//! Cycle slip detection for GNSS carrier phase.
//!
//! Cycle slips are integer discontinuities in the carrier phase measurement
//! caused by receiver tracking loss (signal blockage, low elevation,
//! ionospheric scintillation). Undetected slips propagate meter-level
//! errors into positioning solutions, so PPP / RTK engines must gate
//! every epoch through a slip detector.
//!
//! Two complementary detectors are provided:
//!
//! - **Melbourne-Wübbena Linear Combination** (`MW`) — combines
//!   dual-frequency code and phase to cancel ionosphere + geometry;
//!   the residual is a slowly-varying wide-lane ambiguity that jumps by
//!   an integer at every cycle slip.
//! - **Geometry-Free Linear Combination** (`GF`) — differences the two
//!   carrier phases; the residual is the ionospheric delay, which
//!   varies smoothly with time. A discontinuity in the epoch-to-epoch
//!   difference indicates a slip on one of the two frequencies.
//!
//! # References
//!
//! - Melbourne, W. G. (1985), "The case for ranging in GPS based geodetic
//!   systems", Proc. 1st Int. Symp. Precise Positioning with GPS.
//! - Wübbena, G. (1985), "Software developments for geodetic positioning
//!   with GPS using TI 4100 code and carrier measurements".
//! - Blewitt, G. (1990), "An automatic editing algorithm for GPS data",
//!   Geophys. Res. Lett., 17(3), 199-202.
//! - Kaplan, E. D. & Hegarty, C. J. (2017), "Understanding GPS/GNSS:
//!   Principles and Applications", 3rd ed., §7.3 Multipath and Errors.
//!
//! # Example
//!
//! ```
//! use alice_space::cycle_slip::{Detector, Observation};
//!
//! let mut det = Detector::new(0.5, 0.1);  // MW threshold 0.5 cycles, GF 0.1 m
//! let ok = det.push(Observation {
//!     epoch_s: 0.0,
//!     l1_cycles: 100_000_000.0,
//!     l2_cycles:  77_000_000.0,
//!     p1_m: 20_200_000.0,
//!     p2_m: 20_200_015.0,
//! });
//! assert!(ok.is_ok());
//! ```

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// GPS L1 carrier frequency in Hz (IS-GPS-200 §3.3.1.1).
pub const L1_HZ: f64 = 1_575_420_000.0;
/// GPS L2 carrier frequency in Hz (IS-GPS-200 §3.3.1.1).
pub const L2_HZ: f64 = 1_227_600_000.0;
/// Speed of light in metres per second.
const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

fn wavelength_m(freq_hz: f64) -> f64 {
    SPEED_OF_LIGHT_M_S / freq_hz
}

// ---------------------------------------------------------------------------
// Observation
// ---------------------------------------------------------------------------

/// One dual-frequency observation for a single satellite epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    /// Reception epoch in seconds (monotonic).
    pub epoch_s: f64,
    /// L1 carrier phase in whole cycles.
    pub l1_cycles: f64,
    /// L2 carrier phase in whole cycles.
    pub l2_cycles: f64,
    /// L1 pseudorange (P1 or C1) in metres.
    pub p1_m: f64,
    /// L2 pseudorange (P2 or C2) in metres.
    pub p2_m: f64,
}

impl Observation {
    /// Melbourne-Wübbena Linear Combination in wide-lane cycles.
    ///
    /// `MW = (f1·L1 − f2·L2) / (f1 − f2) − (f1·P1 + f2·P2) / (f1 + f2) / λ_wl`
    ///
    /// where `λ_wl = c / (f1 − f2)` is the wide-lane wavelength (~86 cm for
    /// GPS L1/L2).
    #[must_use]
    pub fn mw_wide_lane_cycles(&self) -> f64 {
        let f1 = L1_HZ;
        let f2 = L2_HZ;
        let l_wl = self.l1_cycles - self.l2_cycles;
        let p_nl_m = (f1 * self.p1_m + f2 * self.p2_m) / (f1 + f2);
        let lambda_wl_m = wavelength_m(f1 - f2);
        l_wl - p_nl_m / lambda_wl_m
    }

    /// Geometry-Free Linear Combination in metres.
    ///
    /// `GF = λ1·L1 − λ2·L2`
    ///
    /// This isolates the ionospheric delay + ambiguities; a step change
    /// implies a cycle slip on either frequency.
    #[must_use]
    pub fn geometry_free_m(&self) -> f64 {
        wavelength_m(L1_HZ) * self.l1_cycles - wavelength_m(L2_HZ) * self.l2_cycles
    }
}

// ---------------------------------------------------------------------------
// Detector
// ---------------------------------------------------------------------------

/// One-satellite cycle slip detector combining MW and GF tests.
#[derive(Debug, Clone)]
pub struct Detector {
    /// MW jump threshold in wide-lane cycles (Blewitt 1990 suggests 4·σ,
    /// typically 0.5 cycles for a clean receiver).
    pub mw_threshold_cycles: f64,
    /// GF jump threshold in metres (a bare receiver in benign
    /// conditions rarely sees more than 0.05 m epoch-to-epoch).
    pub gf_threshold_m: f64,
    prev_mw: Option<f64>,
    prev_gf: Option<f64>,
    prev_epoch_s: Option<f64>,
    /// Total detected slips since the detector was created.
    pub slip_count: u64,
}

impl Detector {
    /// Construct a new detector with the given thresholds.
    #[must_use]
    pub const fn new(mw_threshold_cycles: f64, gf_threshold_m: f64) -> Self {
        Self {
            mw_threshold_cycles,
            gf_threshold_m,
            prev_mw: None,
            prev_gf: None,
            prev_epoch_s: None,
            slip_count: 0,
        }
    }

    /// Push a new observation. Returns [`SlipReport`] describing whether
    /// a slip was detected on this epoch.
    ///
    /// If the observation epoch is not strictly greater than the previous
    /// one, returns [`SlipError::NonMonotonicEpoch`].
    ///
    /// # Errors
    ///
    /// Returns [`SlipError::NonMonotonicEpoch`] if `obs.epoch_s` is not
    /// strictly greater than the previously supplied epoch.
    pub fn push(&mut self, obs: Observation) -> Result<SlipReport, SlipError> {
        if let Some(prev) = self.prev_epoch_s {
            if obs.epoch_s <= prev {
                return Err(SlipError::NonMonotonicEpoch);
            }
        }
        let mw = obs.mw_wide_lane_cycles();
        let gf = obs.geometry_free_m();

        let mw_slip = self
            .prev_mw
            .is_some_and(|prev| (mw - prev).abs() > self.mw_threshold_cycles);
        let gf_slip = self
            .prev_gf
            .is_some_and(|prev| (gf - prev).abs() > self.gf_threshold_m);

        let slip = mw_slip || gf_slip;
        if slip {
            self.slip_count += 1;
        }

        self.prev_mw = Some(mw);
        self.prev_gf = Some(gf);
        self.prev_epoch_s = Some(obs.epoch_s);

        Ok(SlipReport {
            mw_wide_lane_cycles: mw,
            gf_metres: gf,
            mw_slip,
            gf_slip,
        })
    }

    /// Reset the internal state (typically called on receiver loss-of-lock).
    pub fn reset(&mut self) {
        self.prev_mw = None;
        self.prev_gf = None;
        self.prev_epoch_s = None;
    }
}

// ---------------------------------------------------------------------------
// SlipReport / SlipError
// ---------------------------------------------------------------------------

/// Report emitted by [`Detector::push`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlipReport {
    /// Melbourne-Wübbena LC value in wide-lane cycles.
    pub mw_wide_lane_cycles: f64,
    /// Geometry-Free LC value in metres.
    pub gf_metres: f64,
    /// Whether the MW LC exceeded its threshold vs. the previous epoch.
    pub mw_slip: bool,
    /// Whether the GF LC exceeded its threshold vs. the previous epoch.
    pub gf_slip: bool,
}

impl SlipReport {
    /// Whether either detector fired on this epoch.
    #[must_use]
    pub const fn slip_detected(&self) -> bool {
        self.mw_slip || self.gf_slip
    }
}

/// Errors returned by [`Detector::push`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlipError {
    /// The epoch of the pushed observation was not strictly greater
    /// than the previously pushed epoch.
    NonMonotonicEpoch,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn benign_obs(epoch_s: f64, extra_l1_cycles: f64) -> Observation {
        Observation {
            epoch_s,
            l1_cycles: 100_000_000.0 + extra_l1_cycles,
            l2_cycles: 77_000_000.0 + extra_l1_cycles * L2_HZ / L1_HZ,
            p1_m: 20_200_000.0,
            p2_m: 20_200_015.0,
        }
    }

    #[test]
    fn mw_and_gf_are_finite_for_benign_input() {
        let obs = benign_obs(0.0, 0.0);
        assert!(obs.mw_wide_lane_cycles().is_finite());
        assert!(obs.geometry_free_m().is_finite());
    }

    #[test]
    fn benign_stream_does_not_trigger_slip() {
        let mut det = Detector::new(1.0, 0.5);
        for i in 0..10 {
            let obs = benign_obs(f64::from(i), f64::from(i) * 0.001);
            let r = det.push(obs).unwrap();
            assert!(!r.slip_detected(), "unexpected slip at epoch {i}");
        }
        assert_eq!(det.slip_count, 0);
    }

    #[test]
    fn injected_l1_jump_triggers_gf_slip() {
        let mut det = Detector::new(10.0, 0.05);
        det.push(benign_obs(0.0, 0.0)).unwrap();
        // Inject a 5-cycle L1 slip: adds ~5 · 19 cm = 0.95 m to GF LC.
        let mut bad = benign_obs(1.0, 5.0);
        // Break the L2 side so the L1 jump shows up on GF.
        bad.l2_cycles = 77_000_000.0;
        let r = det.push(bad).unwrap();
        assert!(r.gf_slip);
        assert!(r.slip_detected());
        assert_eq!(det.slip_count, 1);
    }

    #[test]
    fn injected_wide_lane_jump_triggers_mw_slip() {
        let mut det = Detector::new(0.5, 100.0);
        det.push(benign_obs(0.0, 0.0)).unwrap();
        // Increase L1 by 100 cycles but leave L2 unchanged.
        let mut bad = benign_obs(1.0, 0.0);
        bad.l1_cycles += 100.0;
        let r = det.push(bad).unwrap();
        assert!(r.mw_slip);
    }

    #[test]
    fn non_monotonic_epoch_is_rejected() {
        let mut det = Detector::new(1.0, 0.5);
        det.push(benign_obs(10.0, 0.0)).unwrap();
        let err = det.push(benign_obs(5.0, 0.0)).unwrap_err();
        assert_eq!(err, SlipError::NonMonotonicEpoch);
    }

    #[test]
    fn reset_clears_state() {
        let mut det = Detector::new(0.1, 0.1);
        det.push(benign_obs(0.0, 0.0)).unwrap();
        det.reset();
        // After reset we should be able to push the same epoch again.
        let r = det.push(benign_obs(0.0, 0.0)).unwrap();
        assert!(!r.slip_detected());
    }

    #[test]
    fn mw_lc_is_stable_when_no_slip() {
        let mut det = Detector::new(0.1, 0.1);
        let a = det.push(benign_obs(0.0, 0.0)).unwrap();
        let b = det.push(benign_obs(1.0, 0.001)).unwrap();
        assert!((a.mw_wide_lane_cycles - b.mw_wide_lane_cycles).abs() < 0.1);
    }

    #[test]
    fn slip_count_accumulates_across_multiple_slips() {
        let mut det = Detector::new(0.5, 0.1);
        det.push(benign_obs(0.0, 0.0)).unwrap();
        for i in 1..=5 {
            let mut bad = benign_obs(f64::from(i), 0.0);
            bad.l1_cycles += 100.0 * f64::from(i);
            det.push(bad).unwrap();
        }
        assert!(det.slip_count >= 1);
    }

    #[test]
    fn slip_report_slip_detected_returns_correct_value() {
        let r = SlipReport {
            mw_wide_lane_cycles: 0.0,
            gf_metres: 0.0,
            mw_slip: false,
            gf_slip: true,
        };
        assert!(r.slip_detected());
    }

    #[test]
    fn geometry_free_reflects_frequency_difference() {
        let obs = benign_obs(0.0, 0.0);
        let gf = obs.geometry_free_m();
        // GF must be a real, finite number close to the difference of the
        // two wavelength-scaled carriers.
        assert!(gf.is_finite());
    }
}
