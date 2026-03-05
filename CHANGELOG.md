# Changelog

All notable changes to ALICE-Space will be documented in this file.

## [0.1.0] - 2026-02-23

### Added
- `orbit` — Keplerian orbital elements, `orbital_period`, `orbital_velocity`, `delta_v_hohmann`, `light_delay_s`, celestial body database
- `propagator` — RK4 two-body orbit propagation (`propagate_rk4`, `propagate_rk4_single`)
- `autonomy` — `TrajectoryModel`, `compute_correction`, `evaluate_decision_tree`, `AutonomyLevel`, `FaultType`
- `comm` — `CommLink`, `ModelDifferential`, `can_transmit` bandwidth check
- `constellation` — `WalkerConstellation` geometry (inclination, planes, phasing)
- `link_budget` — `LinkBudget` with Friis path-loss calculation and margin analysis
- `mission` — `MissionPhase` FSM, `MissionLog` event recording
- FNV-1a shared hash utility
- Zero runtime dependencies (proptest dev-dependency only)
- 122 tests (121 unit + 1 doc-test)
- Clippy pedantic + nursery 0 warnings
- `const fn` for all trivial constructors and accessors
- `mul_add` / `to_degrees` for numerically stable floating-point ops
- proptest property-based tests: orbit, propagator, link_budget, constellation

### Fixed
- Collapsible `if` in `evaluate_decision_tree` (clippy)
