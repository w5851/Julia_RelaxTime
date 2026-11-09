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
- `path_solve_result.json`
- `path_branch_catalog.csv`
- `path_groundstate_selection.csv`
- `path_performance_summary.csv`
- `phase3_backend_comparison.csv`
- `phase3_formalization_review.json`
- `phase3_formalization_review.md`
- `multibranch_summary.json`
- `multibranch_report.md`

Run the aggregate formalization review over the two known multiroot targets:

```sh
julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/review_fixedasymrho_palc_formalization.jl
```

This review runs:

- `T=120 MeV, rho_anchor=0.35`
- `T=130 MeV, rho_anchor=0.80`

and compares `SeedContinuation`, independent pressure-governed `solve_multi`
scan, and multi-branch PALC. The current formalization criterion is diagnostic:
PALC may advance as an opt-in analysis backend when both scenarios expose
distinct branches and explicit ground-state selections, but BifurcationKit stays
out of the root `Project.toml` until broader regression, precompile, and
production opt-in evidence justify changing the stable environment.

### Matched performance pilot (2026-11-08)

The current task uses `.github/workflows/palc-matched-pilot.yml`. It runs the two
existing windows on Actions with 5 common density targets and 8 x 4 quadrature.
`run_matched_pilot.py` launches each method in a separate Julia process;
`run_matched_pilot.jl` records loading separately, then repeats the same fresh
workload six times by default (first call and five warm calls). BifurcationKit remains in
this isolated environment. PNJL numerical pilot execution is on Actions.
The workflow imports the shared dependency stack once before starting timed
subprocesses, so package precompilation cannot consume the numerical timeout.
The numerical workflow is dispatched explicitly for a chosen code revision;
ordinary PR checks remain automatic. Documentation updates do not rerun the pilot.

Methods are multi-anchor PALC, `SeedContinuation`, independent pressure-governed
`solve_multi`, the equilibrium helper used by the combined `trho_asymmetric`
script, a half-step SeedContinuation reference, the legacy nested SeedContinuation,
and a reverse-direction SeedContinuation. The first five and reverse-direction
methods use `fixedasymrho_seed_policy=:single_pass`. Every reported target is
checked at its actual rho with the same explicit residual and physicality gate.
The PALC total includes catalog/high-density anchor discovery, all continuation
branches, and single-seed correction to the common targets. Nearest/interpolated
branch points supply correction seeds only. All target crossings of a folded
branch are considered.

`summary.json` and `summary.md` report parity and timing. A cost ratio is emitted
only after all repetitions and same-target comparisons pass. `processes.json`
also records process wall time and runner compute time, including dependency
loading; Actions setup/install time remains visible in the workflow. First-call
timing includes remaining method JIT, and is not a package-cache-free cold start.
Legacy/new branch differences are reported separately from the new-policy parity
gate. `--require-parity` makes numerical CI fail on incomplete work or mismatches
among the new methods. This differs from the first-stage observational pilot.

The workflow accepts `profile=local|long|higher|transition`: five target points at `8 x 4`,
21 points at `8 x 4`, or 21 points at `24 x 8`. Long/higher windows are
`rho/rho0=0.25:0.01:0.45` and `0.70:0.01:0.90`. `repeat_processes=2` executes
two independent batches; `repetitions` controls identical calls within each process.
`check_candidate_search.jl` verifies the two retained counterexamples, seed
permutations, candidate visibility and default/nested compatibility.
The T120 `transition` profile covers `rho/rho0=0.6:0.1:3.0` with `24 x 8` nodes
and a half-step 49-point reference. It follows the homogeneous root through the
nonmonotone conjugate-chemical-potential region; it does not construct a mixed phase.
Use `verify_only=true` with `case=t120_rho035, profile=local` to rerun only the
candidate, rejection and compatibility checks after changes confined to those contracts.
Each accepted root includes chemical potentials, flavor densities and homogeneous
Helmholtz free energy. T=120 is treated as a first-order-transition temperature
region; branch retention is distinct from a global coexistence construction.

`FixedAsymmetricRho` telemetry separates requests, outer attempts and actual
inner multi-seed NLsolve work. Direct residual evaluations include AD or finite
difference Jacobian work; NLsolve f/g counts use the library's own convention
and are not added to those evaluations. BK callbacks include failed Newton
attempts and retries. Both paths use fixed tensor quadrature, so there is no
adaptive-quadrature workload in this pilot.

This establishes local solver consistency at fixed discretization. The dense
reference is a step-size check, and production parity covers the equilibrium
stage with successful downstream meson steps; neither runs meson densities nor
proves complete branch discovery or converged thermal integration.

`collect_matched_pilot.py --case-dir <downloaded-case> ... --candidate-checks <json>
--output <result.json>` retains timing/counts, comparison errors, one representative
path, branch corrections and differing legacy roots. It checks parity and independent
process results before replacing repeated state arrays with their verified differences.
The current qualified delivery and performance decision are documented in
`docs/analysis/pnjl/algorithmic_feasibility/palc_matched_20261108_analysis/README.md`.
