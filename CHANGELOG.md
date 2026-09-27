# Changelog

All notable changes to ALICE-Space will be documented in this file.

## [Unreleased]

### Changed
- **License: `AGPL-3.0-only` → `AGPL-3.0-only OR LicenseRef-Commercial` (dual-licensed、2026-09-27)** AGPL 側の条件は変更なし (既存 AGPL 利用者への影響ゼロ)、商用という選択肢が追加されただけ SPDX が AGPL 単独だと cargo-deny / FOSSA / SBOM に「商用オプションなし」と見えるため宣言を dual に 変更点: SPDX / `LICENSE` → `LICENSE-AGPL` / `LICENSE-COMMERCIAL.md` (商用トリガー 6 条件 = クローズド製品・商用 SaaS・エッジ / ファームウェア配布・plugin 再配布・プラットフォーム NDA・保証、社内利用は AGPL 側で無償と明記) / README の選択肢表 商用窓口は法人 `contact@extoria.co.jp`

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
