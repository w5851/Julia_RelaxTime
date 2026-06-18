# PNJL PALC Spike

This directory is an isolated BifurcationKit experiment. It does not add
BifurcationKit to the repository root environment and does not change stable
`Models` entrypoints.

Instantiate once:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike -e 'using Pkg; Pkg.instantiate()'
```

Run the PALC branch only:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedmu_palc.jl
```

Run the spike with phase and `run_tmu_scan` references:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedmu_palc.jl --run-reference=true
```

Compare-only entrypoint with reference generation enabled by default:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/compare_fixedmu_palc.jl
```

The scripts stack the repository root project behind this isolated environment
at runtime, so `Models` and its existing dependencies still come from the main
project while BifurcationKit remains isolated here.

Useful options:

```sh
--T-MeV=120 --mu-min-MeV=240 --mu-max-MeV=380 --xi=0.0 --p-num=8 --t-num=4
--output-dir=data/outputs/results/analysis/palc_fixedmu_spike/dev
--run-reference=true
```

Each run writes:

- `palc_fixedmu_branch.csv`
- `tmu_scan_reference.csv` when references are enabled
- `phase_reference_summary.json`
- `comparison_summary.json`
- `comparison_report.md`

The reference phase grid uses the same coarse regression-scale window as the
phase pipeline baseline (`T=120:5:150 MeV`, `rho=0.1:0.1:3.0`) so the script can
usually select a spinodal hint automatically. A missing phase spinodal hint
leaves `continue_recommended=false`; use an explicit `T`/`mu` window around a
known spinodal to evaluate the 10 MeV fold agreement criterion.

## FixedAsymmetricRho

Run the asymmetric density PALC branch only:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedasymrho_palc.jl
```

Run with `run_trho_scan` reference enabled:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/compare_fixedasymrho_palc.jl --run-reference=true
```

Useful options:

```sh
--T-MeV=120 --rho-min=1.05 --rho-max=2.15 --rho-step=0.05
--asym-ud-ratio-target=0.876 --asym-s-target=0.0
--xi=0.0 --p-num=8 --t-num=4
--trho-reverse-rho=true
--output-dir=data/outputs/results/analysis/palc_fixedasymrho_spike/dev
```

Each asymmetric-density run writes:

- `palc_fixedasymrho_branch.csv`
- `trho_scan_reference.csv` when references are enabled
- `comparison_summary.json`
- `comparison_report.md`

The asymmetric-density branch jump metric compares the 8-dimensional solve
state in natural units: the stored `mu_*_MeV` columns are converted back to
fm^-1 before the state-norm jump is computed.

For production-like `trho_asymmetric` diagnostics, use
`--trho-reverse-rho=true` and inspect the `branch_comparison` section in
`comparison_summary.json`. A positive `reference_lower_pressure_count` means the
reference continuation selected roots with lower pressure than the PALC branch
at nearby rho values.

### Multi-branch FixedAsymmetricRho PALC

Use the multi-branch entrypoint when the question is not "can PALC follow one
branch?" but "can PALC preserve all relevant roots and then apply pressure
selection explicitly?"

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedasymrho_multibranch_palc.jl \
  --T-MeV=120 --rho-anchor=0.35 --rho-min=0.30 --rho-max=0.40 \
  --production-source-rho=1.0 --production-source-step=0.05 \
  --p-num=8 --t-num=4 \
  --output-dir=data/outputs/results/analysis/palc_fixedasymrho_spike/multibranch_t120_rho035
```

The anchor discovery stage currently records at least:

- a pressure-governed `solve_multi(...; evaluate_all_attempts=true)` root;
- a production-like high-to-low rho continuation root.

Distinct anchor roots are deduplicated by natural-unit state distance, then each
root starts its own PALC branch. The ground-state selection CSV compares the
branches at sampled rho values and selects the maximum-pressure constrained
branch. This keeps PALC as a continuation strategy while leaving branch
governance explicit.

Each multi-branch run writes:

- `palc_fixedasymrho_anchor_roots.csv`
- `palc_fixedasymrho_multibranch.csv`
- `palc_fixedasymrho_groundstate_selection.csv`
- `multibranch_summary.json`
- `multibranch_report.md`
