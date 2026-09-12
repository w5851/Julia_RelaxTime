# Charged GBU upstream μB sensitivity probe (2026-09-12)

## Purpose

This record tests whether the low-energy `K+/pi+` excess is sensitive to the
upstream baryon chemical potential.  It is a physical-background probe for a
future freeze-out-path study, not a pointwise fit and not a new production
baseline.

The observable remains the charged GBU partial yield on a fixed quark-only BQS
background.  The constraints are `FixedMuBConservedCharges`,
ρQ/ρB = 0.4, ρS = 0, no ΩM or ΣM feedback, and the existing charged
channel mapping (`pi_plus -> K12`, `K_plus -> K45`).

## Worktree and inputs

- worktree: `D:\w\charged-gbu-mub-sensitivity`
- branch: `codex/charged-gbu-mub-sensitivity`
- source HEAD: `b3922e5fc2e1f70150c812859c2464a8685a26f6`
- freeze-out profile: `config/physics/freezeout/default.toml`
- central background/result anchor:
  `data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/freezeout_20260907_v2/`
- selected energies: 3.0, 5.0, 7.7 GeV
- perturbation: fixed baseline `T`, μB = μB,fo + (-50, 0, +50) MeV
- channels evaluated: `pi_plus`, `K_plus`

At ΔμB = 0 the accepted v2 background JSON and production channel JSON are
reused as an explicit central anchor.  At nonzero ΔμB the BQS equilibrium is
solved again from the central seed; μQ, μS, quark masses and the derived meson
chemical potentials are not held fixed.

## Algorithm and numerical scope

The script
`scripts/analysis/relaxtime/scan_charged_gbu_mub_sensitivity.jl` reuses the
infinite-thermal kernel and independent gap/root bookkeeping from
`ChargedGBUResearchWorkflow`/`CausalGBUInfiniteYield`.  The probe settings are

```text
q nodes = 8, qmax = 8 fm^-1
profile mesh = 48, cut nodes = 24, omega nodes = 32, upper = 32 fm^-1
```

This is a sparse diagnostic: it does not certify the full q tail or replace the
accepted 10-point production calculation.  No phase folding, density clipping,
global sign change or Bose-domain substitution is used.

Command:

```powershell
julia --project=. scripts/analysis/relaxtime/scan_charged_gbu_mub_sensitivity.jl `
  --output data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_20260911_v3 `
  --sqrts-list 3.0,5.0,7.7 --mub-offsets -50,0,50 `
  --q-nodes 8 --qmax 8 --mesh 48 --cut-nodes 24 --omega-nodes 32 --upper 32
```

The lowest-energy refinement used `q-nodes=12`, `mesh=96`, `cut-nodes=48`,
`omega-nodes=64`, and the valid upper endpoint 32 fm^-1.  An earlier trial with
upper = 48 fm^-1 failed the existing split-window guard; that failed output was
retained separately and is not used below.

## Results

The accepted sparse diagnostic output is
`data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_20260911_v3/`.
All 9 probe points and 18 channel totals completed without a recorded failure.
The BQS residuals in `backgrounds.csv` remain below `7.7e-12`.

| √sNN (GeV) | T (MeV) | μB,fo (MeV) | R(μB−50) | R(μB) | R(μB+50) | dR/dμB (MeV⁻¹) |
|---:|---:|---:|---:|---:|---:|---:|
| 3.0 | 79.96 | 719.08 | 0.37858 | 0.50965 | 0.69005 | 0.003115 |
| 5.0 | 118.52 | 553.07 | 0.45900 | 0.53871 | 0.63929 | 0.001803 |
| 7.7 | 139.61 | 421.65 | 0.48506 | 0.54557 | 0.61944 | 0.001344 |

Here `R` is the sparse finite-q GBU `K+/pi+` ratio.  The central sparse
values differ from the accepted v2 production values (0.50586, 0.53634,
0.54373) by less than 0.8%, which is adequate for this mechanism probe but not
a new convergence claim.

At 3 GeV, the refined calculation gives

```text
R(muB-50) = 0.376790
R(muB)    = 0.507245
R(muB+50) = 0.686801
dR/dmuB   = 0.003100 MeV^-1
```

The refined values differ from the sparse run by about 0.5%, preserving the
positive slope and its scale.

The channel ledger shows that π+ changes very little while K+ carries almost
all of the response.  At 3 GeV, μK+ changes from about 124.0 to 172.1 MeV as
μB is raised by ±50 MeV; at 5 GeV it changes from 84.6 to 123.5 MeV, and at
7.7 GeV from 69.5 to 102.3 MeV.  This identifies the derived kaon chemical
potential/kaon sector as a high-leverage mechanism candidate, not as a proven
unique cause.

## Provenance

- result manifest:
  `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/muB_sensitivity_20260911_v3/manifest.json`
- manifest SHA-256:
  `440c1d649a862eb60ae8a6941b0d7e27b026f43c058c50a8a9b9bd547b432649`
- ratio CSV SHA-256:
  `bed1022e5e320e8029d694de51593072849f89c553ad537a0f371cacb1da3da1`
- derivative CSV SHA-256:
  `82b9ad191592e37a8ed52d7373212dd76c62db1b7c3882ea4734c3c798ef22ee`
- figure: `muB_sensitivity.png` and `muB_sensitivity.pdf` in the same output directory
- output manifest records `production_authorized=false`, `diagnostic_only=true`,
  `full_q_convergence=false`, and hashes for the current kernel, profile,
  config and reused v2 central files.

## Interpretation and next step

The positive and energy-dependent μB response supports the hypothesis that the
high-μB low-energy background can amplify the low-energy `K+/pi+` ratio,
mainly through the BQS-derived μK+ and kaon density.  It does not by itself
authorize moving the freeze-out curve or changing the production default.

The next physical-path test should vary `T` and μB jointly using a published
alternative freeze-out parameterization or its uncertainty range, re-solving
BQS at every point.  In parallel, the finite-q/q=0 route difference and the
full q-tail resolution should be retained as separate systematic checks.  A
pointwise μB adjustment to force agreement with data, phase rebranching,
negative-density clipping, and meson-feedback activation remain outside this
diagnostic track.
