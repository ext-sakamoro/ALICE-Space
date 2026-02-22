# ALICE-Space

Deep-space communication, autonomous trajectory control, and model-differential protocol for spacecraft operating under extreme latency and minimal bandwidth constraints.

## Features

- Keplerian orbital mechanics: period, vis-viva velocity, Hohmann transfer delta-v, light delay
- RK4 numerical orbit propagator with energy conservation verification
- Two-body gravitational acceleration with reciprocal-optimized computation
- Walker Delta Pattern constellation generator (GPS, Iridium, and arbitrary T/P/F configurations)
- Friis free-space path loss and full link budget analysis (EIRP, C/N0, Eb/N0, margin)
- Model-differential protocol: compact parameter updates transmitted over bandwidth-constrained links
- Autonomous trajectory correction with unit-vector thrust and capped burn duration
- Decision-tree fault detection (sensor anomaly, thruster degradation, comm loss, power fault, attitude error)
- Append-only mission event log with delta-v accounting and phase filtering
- FNV-1a deterministic content hashing on all output types
- Zero external dependencies

## Module Overview

| Module | Key Types | Description |
|---|---|---|
| `orbit` | `OrbitalElements`, `SpacecraftState`, `CelestialBody`, `BodyId` | Keplerian elements, state vectors, and orbital mechanics primitives |
| `propagator` | `TwoBodyAccel`, `propagate_rk4`, `propagate_rk4_single` | 4th-order Runge-Kutta two-body orbit integration |
| `comm` | `CommLink`, `ModelDifferential` | Deep-space link latency, bandwidth windows, and compact parameter differentials |
| `autonomy` | `TrajectoryModel`, `ControlDecision`, `DecisionNode`, `FaultType`, `AutonomyLevel` | Polynomial trajectory models, correction computation, and fault detection |
| `constellation` | `WalkerConstellation`, `WalkerSatellite` | Walker Delta Pattern satellite constellation layout generator |
| `link_budget` | `LinkBudget`, `LinkBudgetResult` | Friis path loss, received power, Eb/N0, and link margin computation |
| `mission` | `MissionLog`, `MissionEvent`, `MissionPhase` | Append-only mission event log with phase and delta-v tracking |

## Quick Start

```rust
use alice_space::{
    orbital_period, orbital_velocity, delta_v_hohmann, light_delay_s,
    propagate_rk4, SpacecraftState,
    WalkerConstellation,
    LinkBudget,
    CommLink, ModelDifferential, can_transmit,
    TrajectoryModel, compute_correction,
    MissionLog, MissionPhase,
};
use std::f64::consts::PI;

const MU_EARTH: f64 = 398600.4418; // km³/s²

// Orbital mechanics
let period = orbital_period(6778.0, MU_EARTH);          // ISS: ~5543 s
let v_leo  = orbital_velocity(6778.0, 6778.0, MU_EARTH); // ~7.67 km/s
let (dv1, dv2) = delta_v_hohmann(6578.0, 42164.0, MU_EARTH); // LEO→GEO

// Light delay
let moon_delay  = light_delay_s(384_400.0);       // ~1.28 s
let mars_delay  = light_delay_s(225_000_000.0);   // ~750 s

// RK4 propagation
let initial = SpacecraftState {
    position_km:    [6778.0, 0.0, 0.0],
    velocity_km_s:  [0.0, (MU_EARTH / 6778.0).sqrt(), 0.0],
    timestamp_ns:   0,
    fuel_kg:        100.0,
};
let trajectory = propagate_rk4(&initial, MU_EARTH, 60.0, 90); // 90-min orbit

// Walker Delta constellation (GPS: 24/6/1)
let gps = WalkerConstellation::new(24, 6, 1, 26559.7, 55.0 * PI / 180.0);
let satellites = gps.generate(); // 24 WalkerSatellite entries

// Link budget (Earth-Moon X-band with DSN 70m dish)
let lb = LinkBudget::default();
let result = lb.compute();
assert!(result.link_closes); // margin_db > 0

// Model differential over deep-space link
let link = CommLink::new(1, 2, 225_000_000.0, 9600.0); // Earth-Mars, 9600 bps
let mut diff = ModelDifferential::new(1, 0);
diff.add_param("thrust_x", 0.05);
diff.add_param("thrust_y", 0.0);
diff.finalize();
let fits = can_transmit(&diff, &link, 60.0); // will it fit in a 60-second window?

// Autonomous trajectory correction
let model = TrajectoryModel::new(
    vec![6800.0, 0.0, 0.0, 0.0, 7.67, 0.0],
    0,
    5_600_000_000_000,
);
let decision = compute_correction(&initial, &model, 0);
// decision.thrust_vector is a unit vector; decision.burn_duration_s <= 300.0

// Mission log
let mut log = MissionLog::new();
log.log_event(MissionPhase::Launch,       0,    2.46, 500.0);
log.log_event(MissionPhase::TransferOrbit, 1000, 1.48, 450.0);
println!("Total dv: {} km/s", log.total_delta_v());
```

## Performance

All hot paths apply the following optimizations:

**Reciprocal pre-computation to avoid repeated division**

`TwoBodyAccel::acceleration` computes `-mu / (r_mag * r2)` as a single reciprocal
factor applied to all three components, replacing three divisions with one.
`compute_correction` similarly pre-computes `1.0 / error_magnitude` before
normalizing the thrust vector.

**`#[inline]` and `#[inline(always)]` on primitives**

`orbital_period`, `orbital_velocity`, `light_delay_s`, `CommLink::latency_s`,
`CommLink::bits_per_window`, `TrajectoryModel::is_valid_at`,
`friis_path_loss_db`, `TwoBodyAccel::acceleration`, and all internal RK4
helpers (`pack`, `deriv`, `scale`, `add`) are marked `#[inline]` or
`#[inline(always)]` to eliminate call overhead at integration sites.

**FNV-1a hash marked `#[inline(always)]`**

The shared `fnv1a` primitive used across `orbit`, `comm`, `constellation`,
`link_budget`, and `mission` is inlined unconditionally, keeping hashing
overhead at near-zero cost.

**Burn duration capped without branching overhead**

`compute_correction` uses `.min(300.0)` to cap burn duration, and `.clamp(0.0, 1.0)`
for confidence, both of which compile to conditional-move instructions on x86-64
and AArch64.

**RK4 state as fixed-size stack arrays**

The propagator represents the 6-element state as `[f64; 6]` (type alias `State6`),
keeping all intermediate RK4 stages on the stack with no heap allocation per step.
`propagate_rk4` pre-allocates the result `Vec` with `with_capacity(steps + 1)`.

**110 tests across 7 modules** (`orbit`: 16, `propagator`: 15, `autonomy`: 24,
`constellation`: 16, `link_budget`: 14, `comm`: 13, `mission`: 12).

## Release Profile

```toml
[profile.release]
opt-level = 3
lto = "fat"
codegen-units = 1
panic = "abort"
strip = true
```

## License

AGPL-3.0-only
