# mixed_check_fine_force_globalfallback_b2_strict_bbo

Temporary Mott-phase scan dataset for evaluating strict mixed-branch behavior under forced global fallback.

## Scope
- profile: `mott_phase_muB0_xi_scan`
- xi list: `[-0.3, -0.15, 0.0, 0.15, 0.3]`
- muB = `0 MeV`
- T window: `120..260 MeV`, step `5 MeV`
- mixed branch policy: `strict_sign_binding`
- force global fallback: `true`

## Optimizer Library (full name)
- fallback global optimizer: `BlackBoxOptim.jl`

## Key Files
- `mott_phase_scan_mixed_check_fine_force_globalfallback_b2_strict_bbo.csv`: raw scan rows
- `mott_phase_scan_mixed_check_fine_force_globalfallback_b2_strict_bbo_derived.csv`: derived indicators
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: script/config hash and run metadata
- `figures/`: mode A/B/AB plots exported from this run

## Purpose
This dataset is used as a focused baseline to compare behavior when forced global fallback is enabled with `BlackBoxOptim.jl`, especially for mixed eta/eta_prime branch continuity checks and anomaly triage.
