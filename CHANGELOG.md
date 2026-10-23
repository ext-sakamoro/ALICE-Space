# Changelog

All notable changes to ALICE-Space will be documented in this file.

## [Unreleased]

### Added
- CI: `.github/workflows/sgp4-golden.yml` が SGP4 の bit golden (`the_walk_gives_the_same_bits_on_every_target`) と Vallado の oracle と TLE の小数の oracle を macOS / Windows / wasm32-wasip1 (wasmtime) で走らせる (src/sgp4.rs・src/tle.rs・試験・Cargo.toml を変える時だけ)、3 本のどれかが走らない・落ちる・通過 6 件未満なら失敗 `ci.yml` の通過件数の下限 `MIN_PASSED` を 346 から 359 に上げ `scripts/preflight.sh` を生成し直した
- `sgp4`: SGP4 / SDP4 (Vallado, Crawford, Hujsak, Kelso 2006 revision, AIAA 2006-6753) for TLE mean elements `Sgp4::new(&MeanElements, GravityModel, OpsMode)` と `Sgp4::propagate(minutes)` が TEME の位置 (km)・速度 (km/s) を返す `MeanElements::from_tle` で `tle::Tle` から作る 三角関数・逆三角関数・冪は全て alice-det-math (全 target で同じ bit) で、SDP4 の共鳴の積分は呼び出しごとに epoch から始めるので結果は呼び出しの順序に依らない 初期化は e が [0, 1) の外、平均運動が正でない、NaN・無限大の field を `Sgp4InitError` で拒む (参照実装は 0 除算か負数の平方根になる入力) 伝播の誤差は参照の code 1-6 (`Sgp4Error`、NaN も範囲外として扱う) 試験 `tests/sgp4_vallado.rs`: Vallado の検証 TLE 33 組 (`SGP4-VER.TLE`) を参照の C++ の出力 (`tcppver.out`) と 666 状態 × 7 field で比べ差は 2e-7 未満 (最大 1.15e-7)、誤差の列 1, 1, 6, 6, 4, 3, 6 も一致、2 file は python-sgp4 (MIT) から持ってきて SHA-256 で固定、全状態の bit の SHA-256 を golden に (x86_64 Linux の CI と arm64 macOS で確認、他の target は CI に lane が無い) `tests/sgp4_degenerate.rs`: 円軌道・範囲外の e・0 と負と NaN の平均運動・NaN の時刻・大きな抗力での落下 依存に alice-det-math 0.5.1 (git rev 9252dd2) を足した
- CI: `Test (public siblings)` job `cargo test --all-features` を走らせ、cargo の失敗、失敗した試験、通った試験が 346 件未満のどれでも失敗 `ci/**` の push でも走る `scripts/preflight.sh` は生成し直して同じ step を持つ

### Changed
- **License: `AGPL-3.0-only OR LicenseRef-Commercial` → `LicenseRef-Proprietary`** (all rights reserved、`publish = false`) 0.7.0 までに公開した版は公開時の license (`AGPL-3.0-only OR LicenseRef-Commercial`) のまま `LICENSE-AGPL` と `LICENSE-COMMERCIAL.md` を外し `LICENSE` を置いた README の license の節と TRADEMARK_NOTICE の参照を合わせた
- **License: `AGPL-3.0-only` → `AGPL-3.0-only OR LicenseRef-Commercial` (dual-licensed、2026-09-27)** AGPL 側の条件は変更なし (既存 AGPL 利用者への影響ゼロ)、商用という選択肢が追加されただけ SPDX が AGPL 単独だと cargo-deny / FOSSA / SBOM に「商用オプションなし」と見えるため宣言を dual に 変更点: SPDX / `LICENSE` → `LICENSE-AGPL` / `LICENSE-COMMERCIAL.md` (商用トリガー 6 条件 = クローズド製品・商用 SaaS・エッジ / ファームウェア配布・plugin 再配布・プラットフォーム NDA・保証、社内利用は AGPL 側で無償と明記) / README の選択肢表 商用窓口は法人 `contact@extoria.co.jp`
- CI: `.github/preflight-exact` を置いた `scripts/preflight.sh` が生成器の出力 (見出しの `generator rev` の版) と byte で一致しなければ送る前の確認が通らない
- 設定 file の comment を中身の説明だけにし、.gitignore を project の build 出力だけにした
- ライセンス全文を公式の AGPL v3 に置き換え
- TRADEMARK_NOTICE が指す licence の file を LICENSE-* の file にした

### Fixed
- `tle::parse`: BSTAR と平均運動の 2 階微分 (仮定小数の field、例 ` -22483-4`) を `仮数 × 10f64.powi(指数)` でなく小数の文字列 (`-0.22483e-4`) を 1 回 parse して作る `powi` の精度は規定されておらず、Vallado の検証 TLE と CelesTrak の 2 組の 70 field のうち 6 field が正しく丸めた値から 1 ulp ずれていた 試験 `tests/tle_assumed_decimal.rs` (70 field 全てが正しく丸めた値と bit 一致)

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
