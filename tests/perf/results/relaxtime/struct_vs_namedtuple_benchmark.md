# Struct vs NamedTuple Performance Benchmark Results

Generated: 2026-01-26T00:24:36.898

## Test Configuration

- Quark masses: u=1.52, d=1.52, s=3.04
- Chemical potentials: u=0.3, d=0.3, s=0.3
- Temperature: T=0.15
- Polyakov loops: Φ=0.5, Φbar=0.5
- Anisotropy: ξ=0.0

## Benchmark Results

| Function | Struct (ms) | NamedTuple (ms) | Ratio (S/NT) | Difference (%) | Status |
|----------|-------------|-----------------|--------------|----------------|--------|
| average_scattering_rate             |      0.3974 |          0.4020 |       0.9885 |          -1.15 | ✓ Pass     |
| relaxation_times                    |    105.4384 |        108.8276 |       0.9689 |          -3.11 | ✓ Pass     |

## Normalization Overhead

The inline normalization helpers have minimal overhead:
- QuarkParams struct normalization: 321.0 ns
- QuarkParams NamedTuple passthrough: 5.0 ns
- ThermoParams struct normalization: 5.0 ns
- ThermoParams NamedTuple passthrough: 10.0 ns

The normalization overhead is negligible compared to the actual computation time.

## Analysis

✓ **All benchmarks passed**: Performance is within 5% tolerance.

The struct migration successfully maintains performance requirements.

## Interpretation

- **Ratio < 1.0**: Struct version is faster than NamedTuple
- **Ratio ≈ 1.0**: Both versions have similar performance
- **Ratio > 1.0**: Struct version is slower than NamedTuple

Performance differences within ±5% are considered acceptable and may be due to:
- Measurement noise
- JIT compilation variations
- Memory layout differences

## Requirements Validation

This benchmark validates requirements:
- **Requirement 10.1**: Relaxation time calculations within 5% performance
- **Requirement 10.2**: Average scattering rate calculations within 5% performance

Note: Requirement 10.3 (cross-section calculations) requires additional work to handle
the A field properly in the struct normalization process. This will be addressed in a
future update.
