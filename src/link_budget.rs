//! Friis free-space link budget calculator
//!
//! Computes received power, SNR, and link margin for deep-space
//! communication links using the Friis transmission equation.
//!
//! Author: Moroya Sakamoto

use crate::fnv1a;
use std::f64::consts::PI;

/// Speed of light in km/s.
const C_KM_S: f64 = 299792.458;

/// Boltzmann constant in dBW/K/Hz.
const K_BOLTZMANN_DBW: f64 = -228.6;

/// Friis free-space path loss in dB.
///
/// FSPL = 20*log10(4π d f / c)
///
/// where d is distance in km, f is frequency in GHz.
#[inline]
pub fn friis_path_loss_db(distance_km: f64, frequency_ghz: f64) -> f64 {
    // Convert to meters and Hz for the formula
    let d_m = distance_km * 1000.0;
    let f_hz = frequency_ghz * 1e9;
    let c_m_s = C_KM_S * 1000.0;
    20.0 * (4.0 * PI * d_m * f_hz / c_m_s).log10()
}

/// Link budget parameters.
#[derive(Debug, Clone, Copy)]
pub struct LinkBudget {
    /// Transmit power in dBW.
    pub tx_power_dbw: f64,
    /// Transmit antenna gain in dBi.
    pub tx_gain_dbi: f64,
    /// Receive antenna gain in dBi.
    pub rx_gain_dbi: f64,
    /// Carrier frequency in GHz.
    pub frequency_ghz: f64,
    /// Link distance in km.
    pub distance_km: f64,
    /// System noise temperature in Kelvin.
    pub system_noise_temp_k: f64,
    /// Data rate in bits per second.
    pub data_rate_bps: f64,
    /// Required Eb/N0 in dB for acceptable BER.
    pub required_eb_n0_db: f64,
    /// Implementation loss in dB.
    pub implementation_loss_db: f64,
}

impl Default for LinkBudget {
    fn default() -> Self {
        Self {
            tx_power_dbw: 20.0,        // 100 W
            tx_gain_dbi: 40.0,         // Deep-space high-gain antenna
            rx_gain_dbi: 70.0,         // DSN 70m dish
            frequency_ghz: 8.4,        // X-band
            distance_km: 384400.0,     // Earth-Moon
            system_noise_temp_k: 25.0, // Cryogenic LNA
            data_rate_bps: 1000.0,     // 1 kbps
            required_eb_n0_db: 3.0,    // BPSK threshold
            implementation_loss_db: 2.0,
        }
    }
}

/// Link budget computation result.
#[derive(Debug, Clone, Copy)]
pub struct LinkBudgetResult {
    /// Free-space path loss in dB.
    pub path_loss_db: f64,
    /// Effective isotropic radiated power in dBW.
    pub eirp_dbw: f64,
    /// Received power in dBW.
    pub received_power_dbw: f64,
    /// System noise power in dBW.
    pub noise_power_dbw: f64,
    /// Carrier-to-noise ratio in dB.
    pub cn0_db_hz: f64,
    /// Achieved Eb/N0 in dB.
    pub eb_n0_db: f64,
    /// Link margin in dB (positive = link closes).
    pub margin_db: f64,
    /// Does the link close? (margin >= 0)
    pub link_closes: bool,
    /// Deterministic content hash.
    pub content_hash: u64,
}

impl LinkBudget {
    /// Compute the full link budget.
    pub fn compute(&self) -> LinkBudgetResult {
        let path_loss = friis_path_loss_db(self.distance_km, self.frequency_ghz);
        let eirp = self.tx_power_dbw + self.tx_gain_dbi;
        let received_power = eirp + self.rx_gain_dbi - path_loss;

        // Noise power: N = k * T * B (in dB: k_dB + 10*log10(T) + 10*log10(B))
        let noise_power = K_BOLTZMANN_DBW + 10.0 * self.system_noise_temp_k.log10();

        // C/N0 = received_power - noise_power (dB-Hz)
        let cn0 = received_power - noise_power;

        // Eb/N0 = C/N0 - 10*log10(data_rate)
        let eb_n0 = cn0 - 10.0 * self.data_rate_bps.log10();

        // Margin = Eb/N0 - required_Eb/N0 - implementation_loss
        let margin = eb_n0 - self.required_eb_n0_db - self.implementation_loss_db;

        // Content hash
        let mut buf = [0u8; 32];
        buf[..8].copy_from_slice(&path_loss.to_bits().to_le_bytes());
        buf[8..16].copy_from_slice(&received_power.to_bits().to_le_bytes());
        buf[16..24].copy_from_slice(&margin.to_bits().to_le_bytes());
        buf[24..32].copy_from_slice(&self.distance_km.to_bits().to_le_bytes());
        let content_hash = fnv1a(&buf);

        LinkBudgetResult {
            path_loss_db: path_loss,
            eirp_dbw: eirp,
            received_power_dbw: received_power,
            noise_power_dbw: noise_power,
            cn0_db_hz: cn0,
            eb_n0_db: eb_n0,
            margin_db: margin,
            link_closes: margin >= 0.0,
            content_hash,
        }
    }

    /// Quick check: can this link achieve the required BER?
    pub fn can_close(&self) -> bool {
        self.compute().link_closes
    }

    /// Maximum data rate achievable with positive margin (binary search).
    pub fn max_data_rate_bps(&self) -> f64 {
        // Eb/N0 = C/N0 - 10*log10(R) >= required + impl_loss
        // => 10*log10(R) <= C/N0 - required - impl_loss
        // => R <= 10^((C/N0 - required - impl_loss)/10)
        let path_loss = friis_path_loss_db(self.distance_km, self.frequency_ghz);
        let eirp = self.tx_power_dbw + self.tx_gain_dbi;
        let received = eirp + self.rx_gain_dbi - path_loss;
        let noise = K_BOLTZMANN_DBW + 10.0 * self.system_noise_temp_k.log10();
        let cn0 = received - noise;
        let max_log_r = cn0 - self.required_eb_n0_db - self.implementation_loss_db;
        10.0_f64.powf(max_log_r / 10.0)
    }
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fspl_earth_moon_x_band() {
        // Earth-Moon at X-band (8.4 GHz, 384400 km)
        let loss = friis_path_loss_db(384400.0, 8.4);
        // Expected ~222.6 dB
        assert!((loss - 222.6).abs() < 1.0, "FSPL: {} dB", loss);
    }

    #[test]
    fn fspl_increases_with_distance() {
        let loss_near = friis_path_loss_db(1000.0, 8.4);
        let loss_far = friis_path_loss_db(100000.0, 8.4);
        assert!(loss_far > loss_near);
    }

    #[test]
    fn fspl_increases_with_frequency() {
        let loss_low = friis_path_loss_db(384400.0, 2.0);
        let loss_high = friis_path_loss_db(384400.0, 32.0);
        assert!(loss_high > loss_low);
    }

    #[test]
    fn link_budget_earth_moon_closes() {
        let lb = LinkBudget::default(); // Earth-Moon defaults
        let result = lb.compute();
        assert!(result.link_closes, "Margin: {} dB", result.margin_db);
        assert!(result.margin_db > 0.0);
    }

    #[test]
    fn link_budget_deep_space_may_not_close() {
        let lb = LinkBudget {
            distance_km: 778_500_000.0, // Jupiter distance
            data_rate_bps: 1_000_000.0, // 1 Mbps
            tx_power_dbw: 10.0,         // Low power
            tx_gain_dbi: 20.0,          // Small antenna
            rx_gain_dbi: 40.0,          // Moderate ground station
            ..Default::default()
        };
        let result = lb.compute();
        assert!(!result.link_closes, "Margin: {} dB", result.margin_db);
    }

    #[test]
    fn eirp_calculation() {
        let lb = LinkBudget {
            tx_power_dbw: 20.0,
            tx_gain_dbi: 40.0,
            ..Default::default()
        };
        let result = lb.compute();
        assert!((result.eirp_dbw - 60.0).abs() < 1e-10);
    }

    #[test]
    fn content_hash_deterministic() {
        let lb = LinkBudget::default();
        let r1 = lb.compute();
        let r2 = lb.compute();
        assert_eq!(r1.content_hash, r2.content_hash);
        assert_ne!(r1.content_hash, 0);
    }

    #[test]
    fn can_close_shortcut() {
        let lb = LinkBudget::default();
        assert!(lb.can_close());
    }

    #[test]
    fn max_data_rate_positive() {
        let lb = LinkBudget::default();
        let max_rate = lb.max_data_rate_bps();
        assert!(
            max_rate > lb.data_rate_bps,
            "Max rate {} should exceed configured {}",
            max_rate,
            lb.data_rate_bps
        );
    }

    #[test]
    fn higher_power_gives_more_margin() {
        let lb_low = LinkBudget {
            tx_power_dbw: 10.0,
            ..Default::default()
        };
        let lb_high = LinkBudget {
            tx_power_dbw: 30.0,
            ..Default::default()
        };
        let r_low = lb_low.compute();
        let r_high = lb_high.compute();
        assert!(r_high.margin_db > r_low.margin_db);
    }

    #[test]
    fn noise_power_reasonable() {
        let lb = LinkBudget::default();
        let result = lb.compute();
        // Noise power at 25K: k_dB + 10*log10(25) ≈ -228.6 + 14.0 = -214.6 dBW/Hz
        assert!((result.noise_power_dbw - (-214.6)).abs() < 1.0);
    }

    #[test]
    fn max_data_rate_exceeds_current_when_link_closes() {
        let lb = LinkBudget::default();
        assert!(lb.can_close());
        let max_rate = lb.max_data_rate_bps();
        assert!(max_rate > lb.data_rate_bps);
    }

    #[test]
    fn content_hash_changes_with_distance() {
        let lb1 = LinkBudget {
            distance_km: 384400.0,
            ..Default::default()
        };
        let lb2 = LinkBudget {
            distance_km: 778_500_000.0,
            ..Default::default()
        };
        let r1 = lb1.compute();
        let r2 = lb2.compute();
        assert_ne!(r1.content_hash, r2.content_hash);
    }

    #[test]
    fn fspl_doubles_distance_adds_6db() {
        // Doubling distance adds ~6 dB of path loss (20*log10(2) ≈ 6.02)
        let loss1 = friis_path_loss_db(10000.0, 8.4);
        let loss2 = friis_path_loss_db(20000.0, 8.4);
        assert!(
            (loss2 - loss1 - 6.02).abs() < 0.1,
            "6dB rule: diff={}",
            loss2 - loss1
        );
    }
}
