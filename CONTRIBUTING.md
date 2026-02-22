# Contributing to ALICE-Space

## Build

```bash
cargo build
```

## Test

```bash
cargo test
```

## Lint

```bash
cargo clippy -- -W clippy::all
cargo fmt -- --check
cargo doc --no-deps 2>&1 | grep warning
```

## Design Constraints

- **Zero external dependencies**: all orbital mechanics and link budget math are self-contained.
- **Model-differential protocol**: transmit mathematical model updates (coefficients) instead of raw telemetry.
- **Deterministic propagation**: RK4 integrator produces bit-exact results for the same inputs.
- **Fault detection**: decision tree evaluates sensor readings against thresholds — highest severity wins.
- **Deep-space latency tolerance**: all control decisions are autonomous-capable without ground link.
