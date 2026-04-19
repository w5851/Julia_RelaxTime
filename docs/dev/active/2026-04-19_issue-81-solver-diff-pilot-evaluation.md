# Issue #81 Solver/Diff Pilot Evaluation

## Scope

- Script chain: `scripts/pnjl/calculate_derivatives.jl`
- Analysis chain: `scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl`
- Unified pilot service:
  - `Models.build_pilot_diff_context`
  - `Models.eval_pilot_derivatives`

## What Changed

### Script chain pilot replacement

- `calculate_derivatives.jl` Step-1 now maps key thermodynamic derivatives from pilot Jacobian path.
- CSV header and CLI remain unchanged.
- Preserved behaviors:
  - success/failure counters
  - failed grid points are counted but not written as data rows

Mapped fields (pilot path):

- `dP_dT` <- `pressure__dT_fm`
- `dP_dmu` <- `pressure__dmu_fm`
- `dEpsilon_dT` <- `energy__dT_fm`
- `dEpsilon_dmu` <- `energy__dmu_fm`
- `dn_dT` <- `rho_norm__dT_fm * ρ0`
- `dn_dmu` <- `rho_norm__dmu_fm * ρ0`

### Analysis chain pilot instrumentation

- `t190_sigma_chain_decomposition.jl` adds pilot derivative cache in `StatePoint` and appends pilot-only CSV columns.
- Core sigma/blocking physics flow is unchanged.

Appended detail columns:

- `pilot_dpressure_dT`
- `pilot_dentropy_dT`

Appended summary columns:

- `pilot_dpressure_dT_A`
- `pilot_dpressure_dT_B`
- `pilot_dentropy_dT_A`
- `pilot_dentropy_dT_B`

## Numerical Consistency Snapshot

Command used (single-point baseline vs pilot extraction):

```bash
julia --project=. -e 'include("src/models/Models.jl"); ...'
```

### Point A (`T=150 MeV`, `mu=0 MeV`, `xi=0`, `p_num=8`, `t_num=4`)

- `dP_dT`: abs `6.72e-10`, rel `4.03e-9`
- `dE_dT`: abs `1.57e-8`, rel `7.41e-9`
- `dn_dmu`: abs `3.27e-12`, rel `4.28e-10`
- near-zero channels (`dP_dmu`, `dE_dmu`, `dn_dT`) show larger relative ratios dominated by denominator effect (absolute errors remain small)

Timing (cold-ish single-point wall-time):

- baseline thermo-derivatives path: `68.83 s`
- pilot Jacobian path: `17.97 s`

### Point B (`T=150 MeV`, `mu=100 MeV`, `xi=0`, `p_num=8`, `t_num=4`)

- `dP_dT`: abs `1.70e-9`, rel `6.81e-9`
- `dP_dmu`: abs `2.04e-10`, rel `9.48e-9`
- `dE_dT`: abs `3.43e-8`, rel `9.39e-9`
- `dE_dmu`: abs `3.24e-10`, rel `8.79e-10`
- `dn_dT`: abs `1.25e-9`, rel `8.83e-9`
- `dn_dmu`: abs `8.80e-12`, rel `2.91e-10`

Timing:

- baseline thermo-derivatives path: `67.91 s`
- pilot Jacobian path: `18.33 s`

Interpretation:

- For non-near-zero channels, pilot and baseline are numerically consistent at ~`1e-8` relative scale in sampled points.
- For near-zero channels at `mu=0`, use absolute tolerance gating instead of pure relative tolerance.

## Runtime Cost (pilot chain evidence)

Script tiny-grid smoke:

```bash
julia --project=. scripts/pnjl/calculate_derivatives.jl --T_min=150 --T_max=150 --T_step=10 --mu_min=0 --mu_max=0 --mu_step=50 --p_num=8 --t_num=4 --output_dir=data/outputs/results/pnjl/derivatives_pilot_smoke
```

Observed:

- Step-1 success `1/1`, failure `0`, ~`74.3 s`
- Step-2 success `1/1`, failure `0`, ~`21.1 s`
- outputs generated:
  - `data/outputs/results/pnjl/derivatives_pilot_smoke/derivatives_xi0.0.csv`
  - `data/outputs/results/pnjl/derivatives_pilot_smoke/bulk_viscosity_xi0.0.csv`

## Verification Evidence

### Focused pilot adapter unit

```bash
julia --project=. -e 'ENV["UNIT_FILES"]="models/solver/test_solver_diff_pilot_adapters.jl"; include("tests/unit/runtests.jl")'
```

- Result: `Pass 17 / Total 17`

### Integration core (issue-required)

```bash
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
```

- Result summary: `Pass 130 / Total 130`
- Note: run log contains a known orchestrator warning/error text (`fallback detected while --fail-on-fallback=true`) but wrapper summary is passing.

### Validation smoke (issue-required)

```bash
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
```

- Result: `Pass 24 / Total 24`

### Unit smoke (supplemental)

```bash
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
```

- Result: `Pass 150 / Total 150`

### Docs consistency (api/doc sync)

```bash
julia --project=. scripts/dev/check_docs_consistency.jl
```

- Result: `OK`

## Maintainability Assessment

- New unified service entrypoints reduce ad-hoc derivative path duplication.
- Script and analysis chains now share one Jacobian contract surface (`target + param + by_name mapping`).
- Risk concentrated in pilot adapter naming/shape contract and therefore easier to test centrally.

## Phase-2 Decision (Actionable)

### Candidate callbacks to fold into `ProblemSpec`

Recommend minimal fold-in set:

- context construction normalization (`theta` alias normalization, `spec_override` forwarding)
- Jacobian shape/finite validation for `N x P`
- stable target registry for core thermo outputs (`pressure`, `energy`, `rho_norm`, `entropy`)

### Keep in outer diff layer (for now)

- script-/analysis-specific column mappings and naming conventions
- pilot-only convenience payload assembly (`by_name` flattening policy)
- chain-specific diagnostics and summary postprocess

### Migration template for follow-up chains

1. Build solver result with chain-local runtime knobs.
2. Create pilot diff context via `build_pilot_diff_context`.
3. Evaluate Jacobian via `eval_pilot_derivatives(targets, params)`.
4. Map through explicit key checks (no silent fallback).
5. Preserve original output contract; append pilot/provenance fields only.
6. Add focused unit contract test before chain integration.

## Outstanding Notes

- `calculate_derivatives.jl` still has structural duplication between Step-1/Step-2 loops by design choice for now; behavior and output contract are stable.
- For near-zero derivative channels, future regression governance should enforce mixed tolerance (`atol + rtol`) rules.
