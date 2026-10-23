# Changelog

All notable changes to ALICE-Space will be documented in this file.

## [Unreleased]

### Added
- `interp`: Lagrange (`lagrange`、位置と補間多項式の微分としての速度) と Hermite (`hermite`、位置と速度を節点で一致) 窓は時刻に中心を合わせた節点 (表の端では内側に寄せる、CCSDS 502.0 の次数の使い方) で、積は窓の最初の節点からの相対時刻で作る 試験 `tests/interp_oracle.rs`: 次数 n の Lagrange が n 次以下の多項式を値と微分で再現 (1 から 10 次)、k 節点の Hermite が 2k−1 次以下を再現、円軌道で中心を合わせた窓の誤差限界の内側、退化入力
- `ccsds_odm`: CCSDS 502.0-B-3 の OEM (`parse_oem_kvn`、segment ごとの metadata・使える時間窓・`INTERPOLATION` / `INTERPOLATION_DEGREE`・7 列か 10 列の行、covariance の block は読み飛ばす) と OPM (`parse_opm_kvn`、状態 vector と Kepler 要素) を KVN で読む `OemSegment::state_at` が message の補間で状態を返す epoch は整数の日と日内の秒 (`OdmEpoch`) で持ち segment 内の相対時刻は精度を落とさない UTC の segment の閏秒 (23:59:60) は拒む (時刻系の表が要る) 試験 `tests/ccsds_odm_oracle.rs`: 円軌道を 60 s ごとに書いた OEM (Lagrange 7 次 / Hermite 5 次、2 segment、covariance と COMMENT つき) が行の間で閉形式を誤差限界内で再現、1 行おきに抜いた表が抜いた行を再現、日付をまたぐ segment、OPM、壊れた入力 17 種
- `sp3`: `Sp3File::interpolate(prn, unix_seconds, degree)` が `interp::lagrange` で衛星の位置を補間する (位置が全て 0 の欠測の行は節点にしない) 試験 `tests/sp3_interp.rs`: GPS 型の円軌道を 15 分ごとに 1e-6 km で書いた SP3 を 9 次で補間して閉形式を 5e-6 km 以内で再現、欠測の epoch をまたいでも同じ
- `omm`: CCSDS 502.0-B-3 の OMM (Orbit Mean-Elements Message) を SGP4 の要素として読む `parse_kvn` (KVN、単位注記 `[deg]` 等と `COMMENT` を読み飛ばし、`EPOCH` の後の `OBJECT_NAME` / `EPOCH` で次の要素へ) / `parse_xml` (NDM XML の `<segment>` ごと) / `parse_json` (CelesTrak の GP JSON、配列か 1 object) / `parse_csv` (CelesTrak の GP CSV) が `Omm` を返し、`Omm::to_mean_elements` が `sgp4::MeanElements` を作る 単位の換算は TLE と同じ式、epoch は暦の時刻 (`YYYY-MM-DDThh:mm:ss.f` と通日の `YYYY-DDDThh:mm:ss.f`) を年と通日にして TLE と同じ epoch の算術に通すので、同じ要素の OMM と TLE は SGP4 の入力が bit 一致する `CENTER_NAME` / `REF_FRAME` / `TIME_SYSTEM` / `MEAN_ELEMENT_THEORY` は EARTH / TEME / UTC / SGP4 以外を拒む、必須 key の欠落・範囲外の値・閏秒 (60 秒) は key を名指す誤り 試験 `tests/omm_vs_tle.rs`: CelesTrak の MARIO (XML / CSV と、同じ値の KVN / JSON) と Vanguard 1 の JSON が TLE と全 field bit 一致、Vallado の検証 TLE 33 組を列ごとに KVN に書き写したものも 330 field bit 一致 `tests/omm_degenerate.rs`: 欠落 key 11 種、不正な epoch 8 種、非対応の metadata 4 種、範囲外・非数 7 種、UNKNOWN の id、2 組の KVN、各形式の壊れた文 手で入れた変異 7 本は全て red
- CI: `.github/workflows/sgp4-golden.yml` が SGP4 の bit golden (`the_walk_gives_the_same_bits_on_every_target`) と Vallado の oracle と TLE の小数の oracle を macOS / Windows / wasm32-wasip1 (wasmtime) で走らせる (src/sgp4.rs・src/tle.rs・試験・Cargo.toml を変える時だけ、push では main との共通祖先からの差分で判定し、手動でも起動できる、`tests/data/**` は改行を変換せずに checkout する)、3 本のどれかが走らない・落ちる・通過 6 件未満なら失敗 `ci.yml` の通過件数の下限 `MIN_PASSED` を 346 から 359 に上げ `scripts/preflight.sh` を生成し直した
- `sgp4`: SGP4 / SDP4 (Vallado, Crawford, Hujsak, Kelso 2006 revision, AIAA 2006-6753) for TLE mean elements `Sgp4::new(&MeanElements, GravityModel, OpsMode)` と `Sgp4::propagate(minutes)` が TEME の位置 (km)・速度 (km/s) を返す `MeanElements::from_tle` で `tle::Tle` から作る 三角関数・逆三角関数・冪は全て alice-det-math (全 target で同じ bit) で、SDP4 の共鳴の積分は呼び出しごとに epoch から始めるので結果は呼び出しの順序に依らない 初期化は e が [0, 1) の外、平均運動が正でない、NaN・無限大の field を `Sgp4InitError` で拒む (参照実装は 0 除算か負数の平方根になる入力) 伝播の誤差は参照の code 1-6 (`Sgp4Error`、NaN も範囲外として扱う) 試験 `tests/sgp4_vallado.rs`: Vallado の検証 TLE 33 組 (`SGP4-VER.TLE`) を参照の C++ の出力 (`tcppver.out`) と 666 状態 × 7 field で比べ差は 2e-7 未満 (最大 1.15e-7)、誤差の列 1, 1, 6, 6, 4, 3, 6 も一致、2 file は python-sgp4 (MIT) から持ってきて SHA-256 で固定、全状態の bit の SHA-256 を golden に (x86_64 Linux の CI と arm64 macOS で確認、他の target は CI に lane が無い) `tests/sgp4_degenerate.rs`: 円軌道・範囲外の e・0 と負と NaN の平均運動・NaN の時刻・大きな抗力での落下 依存に alice-det-math 0.5.1 (git rev 9252dd2) を足した
- CI: `Test (public siblings)` job `cargo test --all-features` を走らせ、cargo の失敗、失敗した試験、通った試験が 346 件未満のどれでも失敗 `ci/**` の push でも走る `scripts/preflight.sh` は生成し直して同じ step を持つ

### Changed
- `sgp4::epoch_days_since_1950`: 1 月 1 日までの日数を Julian date の式でなく整数の日数から作る (結果の bit は同じ、参照と同じ Julian date 経由の丸めは保つ) OMM の epoch も同じ関数を通す
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
