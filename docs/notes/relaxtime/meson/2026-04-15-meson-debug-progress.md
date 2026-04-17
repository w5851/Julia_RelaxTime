# Meson Workflow Debug Progress (2026-04-15)

## Scope

- Goal: explain Julia vs Fortran meson-target mismatches with evidence.
- Focus: separate hard equation differences from solver/path/branch semantics.

## Confirmed findings

1. Upstream branch/path issue was real but is no longer dominant.
- After adding hard constraint `m_s > m_u` in `src/models/solver/runtime/ConstraintSolverCommon.jl`, the previous three severe equilibrium branch-jump points were addressed (`(0,200)`, `(0,240)` returned to expected branch neighborhood; `(600,180)` high-mass branch filtered).
- Legacy all8 mismatch count improved from about `37` to `20`.

2. Hard equation mismatch is largely excluded for current failing set.
- Script: `scripts/analysis/relaxtime/meson_mixed/compare_equation_semantics_hard_sources.jl`
- On current `20` failing rows:
  - mixed channels (`eta/eta_prime/sigma/sigma_prime`): Julia vs Fortran equation core is numerically equivalent (`~1e-16`).
  - nonmixed channels have one implementation-form difference in imaginary equation component:
    - Fortran `eq_meson.f90`: `fvec(2)=pi_im`
    - Julia `MesonMass.jl`: `Im[1-4KΠ] = -4K*pi_im`
  - For current failed nonmixed rows, `K` is not near zero (`~0.161-0.226`), so root sets are equivalent in practice for these points (scaling difference, not root-condition difference).

3. Current bottleneck is solver-side/root-selection-side behavior.
- Many failed points show `root_quality=:bad` and `residual=Inf` after attempts `newton|trust_region|newton|newton|newton`.
- This means both Newton and trust-region were attempted but still did not produce a good candidate.

4. Mixed sign binding alone does not improve failure count.
- Added mode: `mixed_branch_align=:strict_sign_binding` in `src/models/workflows/MesonMassWorkflow.jl`.
- Added validation coverage in `tests/validation/relaxtime/test_mixed_branch_identity_alignment_experiment.jl`.
- Legacy all8 comparison:
  - default (`identity_track_label_output`): `20 fail` (`mixed=13`, `nonmixed=7`)
  - strict sign binding: `20 fail` (`mixed=13`, `nonmixed=7`)
- Conclusion: label-swap suppression alone is not enough; root search quality remains limiting factor.

## Downstream equations (for quick reference)

- nonmixed (`pi/K/sigma_pi/sigma_K`): solve in `(M, Gamma)`
  - `Re[1 - 4KΠ(M,Gamma)] = 0`
  - `Im[1 - 4KΠ(M,Gamma)] = 0`
- mixed (`eta/eta_prime/sigma/sigma_prime`): solve branch eigen-equation in `(M, Gamma)`
  - `F_minus = (M00 + M88) - sqrt((M00 - M88)^2 + 4M08^2) = 0`
  - `F_plus  = (M00 + M88) + sqrt((M00 - M88)^2 + 4M08^2) = 0`

## Reuse-first plan (solver architecture)

Before introducing new meson-specific solver machinery, reuse existing `src/models/solver` patterns:

1. `src/models/solver/runtime/GenericRootEngine.jl`
- Existing reusable capabilities: primary/fallback methods, multiseed, continuation state, quality tagging, diagnostics.
- Already used by meson workflow via callback mode (`solve_root_with_policy`).

2. `src/models/solver/governance/CandidateGovernance.jl`
- Existing reusable capabilities: attempt-pool execution, candidate normalization, quality tags, selector contracts.
- Can be reused to make meson candidate selection explicit and auditable.

3. `src/models/solver/orchestrator/ProblemSpecOrchestrator.jl`
- Existing reusable capabilities: governed attempt planning, hard-constraint integration, fallback routing.
- Design direction: map meson root solve into the same governed selection semantics where feasible.

## Next technical direction

Prioritize solver robustness improvements (shared across mixed and nonmixed):

1. Two-stage root pipeline
- Stage A: cheap candidate generation (seed pool + bounded coarse scan).
- Stage B: local refinement (Newton/TR) with strict acceptance criteria.

2. Residual scaling/preconditioning
- Normalize `Re/Im` residual components to avoid one component dominating step acceptance.

3. Continuation + governed candidate selection
- Keep continuity seeds, but always evaluate a governed backup candidate set and select by explicit rule.

4. Diagnostics-first policy
- Persist per-point attempt traces (seed source, method, residual, selected reason) for deterministic triage.

## Status

- This note captures conclusions up to current session state and should be updated before further solver changes.

## Update: governance adapter pass (same day)

Implemented a zero-intrusion solver-core approach at workflow boundary:

- Added workflow-level meson governance probing and selection in `src/models/workflows/MesonMassWorkflow.jl`.
  - Keep existing `RootPolicy` execution.
  - When selected root is not `:good`, perform governed candidate probe over method/seed pool.
  - Reuse exported governance API (`execute_governance_selector`) with residual-min selector.
  - Add diagnostics fields:
    - `governance_candidate_count`
    - `governance_selection_reason`

- Added strict sign mode:
  - `mixed_branch_align=:strict_sign_binding` (fixed mixed-sign mapping, no label swap).

- Added/updated tests:
  - `tests/unit/models/test_meson_mass_workflow.jl`: diagnostics fields contract.
  - `tests/validation/relaxtime/test_mixed_branch_identity_alignment_experiment.jl`: strict-sign behavior.

### Effect on legacy all8

- Before governance adapter (after upstream hard-constraint fix): `20` fail.
- After governance adapter: `4` fail.
- Remaining fails are all mixed `eta` rows mapped to Julia `eta_prime` side in legacy comparator, with one clear sign issue (`~ -700 vs +700`) and three residual mixed-branch offsets.

### Interpretation

- This supports the hypothesis that most remaining divergence is in mixed-branch identity/selection semantics rather than hard equation mismatch.

## Update: mixed cleanup sprint (same day, later)

Completed two-step mixed-side cleanup in workflow boundary layer:

1) Semantic normalization guard
- For mixed channels, negative mass outputs are normalized to positive magnitude with positive width magnitude.
- This removes pure sign-mirror artifacts (e.g., `-700` vs `+700` mirror case).

2) Mixed lightweight solver enhancement
- Governance probe for mixed channels expanded with extra seeds around selected mass and threshold mass.
- Mixed candidate scoring changed from generic residual-only to label-bound branch score (`eta/sigma -> minus`, `eta_prime/sigma_prime -> plus`).
- For no-pass constraint fallback, selector now applies a mixed-specific metric with mild mass/gamma regularization.

### Current status after mixed sprint

- Legacy all8 reduced further from `2 fail` to `1 fail`.
- Remaining single failure:
  - `(muB=600, T=240, eta)` legacy target mapped to Julia `eta_prime` still ~10% high (`~623.49` vs `~562.46`).
  - Point diagnostics indicate all governed mixed candidates remain non-good (`governance_selection_reason=no_candidate_passed_constraints`), i.e. this is currently a hard mixed near-degeneracy/root-basin competition point under present residual/solver policy.

### Practical takeaway

- The major mismatch body has been eliminated without touching solver core.
- Residual one-point mismatch appears to require either:
  - a targeted mixed paired-root search objective (joint eta/eta' optimization), or
  - a small tolerance/expectation governance decision for this known near-degenerate point.

## Update: legacy baseline governance fix (final)

Per verification request, we checked whether the last remaining legacy mismatch point
`(muB=600, T=240, eta)` is a high-quality numerical root in Fortran output.

Verification evidence:

- Fortran root diagnostics row (`root_diag.dat`, id=3 at `mu_loop=600`, `T=240`) shows:
  - `newt1_converged = 0`
  - `newt1_resid ≈ 0.05126`
- Therefore this point is not a converged high-quality root, but a low-quality
  continuation-path value retained in legacy export.

Action taken:

- Marked that baseline row as low-quality excluded target:
  - file: `tests/validation/data/targets/relaxtime/legacy/meson/legacy_meson_scan_fortran_muB600_v1.csv`
  - status changed to `excluded_low_quality_nonconverged` for row `fortran_eta_240p0_xi_0p0_51`.
- Updated validator to skip excluded low-quality rows:
  - `tests/validation/relaxtime/test_legacy_fortran_meson_numeric_mu0_mu600_all8_targets.jl`
- Updated target data README governance note:
  - `tests/validation/data/targets/relaxtime/legacy/meson/README.md`
- Updated analysis script to skip excluded rows:
  - `scripts/analysis/relaxtime/meson_mixed/check_meson_target_backsub_residuals.jl`

Result after governance correction:

- Legacy all8 validation: pass (`254/254`).
- Literature targets: pass (`0` fails).
- Back-sub diagnostic script: both literature and legacy report `0` failed rows.

## Update: scripts/relaxtime Mott pipeline mixed integration check

### Objective

Verify that script-layer Mott pipeline reflects latest meson workflow fixes, and add optional mixed output/plot path to avoid validation-set-only bias.

### What was changed

1) `scripts/relaxtime/run_mott_phase_scan.jl`
- Added CLI/config switch `--include-mixed` (default `false`, non-breaking for existing nonmixed workflows).
- When enabled, compute/output mixed observables with same runtime settings as nonmixed:
  - masses: `M_eta`, `M_eta_prime`, `M_sigma`, `M_sigma_prime`
  - widths: `Gamma_eta`, `Gamma_eta_prime`, `Gamma_sigma`, `Gamma_sigma_prime`
  - residuals/qualities/methods for each mixed observable

2) `scripts/relaxtime/run_mott_phase_derived_csv.jl`
- Added optional derived column `M_eta_plus_M_eta_prime` when mixed columns exist.

3) `scripts/relaxtime/run_mott_phase_plot_modes.jl`
- Added `--include-mixed` switch.
- Plot observable list now extends with mixed observables only when flag is on.

### One-shot run (mixed enabled)

Command profile used (quick smoke-sized scan):

- `T: 120 -> 260 MeV, step 20 MeV`
- `xi list: profile default`
- `muB=0`
- `p_num=8, t_num=4, max_iter=25`
- output dir: `data/outputs/results/relaxtime/mott_phase/mixed_check`

Artifacts generated:

- scan CSV: `data/outputs/results/relaxtime/mott_phase/mixed_check/mott_phase_scan_mixed.csv`
- derived CSV: `data/outputs/results/relaxtime/mott_phase/mixed_check/mott_phase_scan_mixed_derived.csv`
- figures root: `data/outputs/results/relaxtime/mott_phase/mixed_check/figures`

Quick sanity summary:

- mixed columns all finite on this run (`40/40` finite per mixed observable).
- no NaN/Inf observed in mixed output fields.
- coarse continuity check (xi=0.0, step=20 MeV) shows largest jumps around `T≈200 MeV`, consistent with crossover/transition-sensitive region behavior rather than random script-level corruption.

### Conclusion

- Script-layer Mott pipeline is now connected to latest workflow-side meson fixes.
- Mixed path is optional and backward-compatible by default.
- Mixed outputs and plots are now available for ongoing robustness monitoring beyond fixed validation targets.

## Update: finer-grid mixed scan check (`T_step=5 MeV`)

### Run setup

- outdir: `data/outputs/results/relaxtime/mott_phase/mixed_check_fine`
- scan command: same as mixed-check run but `T_step=5`
- plotting command: same with `--include-mixed`

Generated artifacts:

- `data/outputs/results/relaxtime/mott_phase/mixed_check_fine/mott_phase_scan_mixed_fine.csv`
- `data/outputs/results/relaxtime/mott_phase/mixed_check_fine/mott_phase_scan_mixed_fine_derived.csv`
- `data/outputs/results/relaxtime/mott_phase/mixed_check_fine/figures/`

### Anomaly summary script

- Added: `scripts/analysis/relaxtime/meson_mixed/summarize_mixed_scan_anomalies.jl`
- Purpose: report per-xi temperature-step absolute jumps for mixed observables above threshold.

Observed with threshold `0.25`:

- Jump clusters are not random; they are concentrated in transition-sensitive temperature bands:
  - around `T≈185-210 MeV` (notably for `M_sigma`, `Gamma_eta`),
  - and around `T≈240-255 MeV` for higher positive xi slices (notably eta/eta' and sigma/sigma').
- No NaN/Inf spread pattern indicating script wiring faults; behavior looks like model/solver sensitivity in known difficult regions rather than pipeline breakage.

### Takeaway

- Script-layer mixed integration is functioning and exposes expected sensitivity regions under finer scan.
- Next robustness pass should focus on solver-stability diagnostics in those clustered T windows (instead of generic script changes).

## Update: jump root-cause analysis for scripts scan (same day, follow-up)

User observed large jumps in recent figures and requested comparison against historical
`reference_100_300_fine` outputs (where pi/K and sigma_pi/sigma_K looked smoother).

### Evidence and comparison

1) Config mismatch was significant in earlier ad-hoc checks:

- `reference_100_300_fine` uses `p_num=12, t_num=6, max_iter=40, T_step=1`.
- recent `mixed_check_fine` used `p_num=8, t_num=4, max_iter=25, T_step=5`.

2) After controlling configs, big discrepancy remained for script scans:

- Script nonmixed fine (before fix, models-default path) still showed large pi jumps
  (`M_pi` max jump ~3.58 fm^-1, many `>=0.15` events).
- Reference data in same T window, decimated to 5 MeV, had much smaller pattern
  (no giant 3+ fm^-1 flips; largest jumps mainly threshold-onset widths).

3) Backend and seed-path isolation test identified root cause:

- Running workflow point-by-point with `solver_backend=:models` and **no equilibrium continuation seed** reproduces large jump pattern.
- Running with `solver_backend=:models` + **equilibrium continuation seed (`x_state` carry-over)** removes giant flips and collapses stats to reference-like level.
- `solver_backend=:legacy` naturally follows continuation-like behavior and is stable under same sweep.

### Root cause conclusion

The jump issue is primarily **scan-loop equilibrium seeding policy**, not mixed-column wiring itself:

- `run_mott_phase_scan.jl` carried meson-level seeds but did not carry equilibrium seed across T for each xi.
- Under models-default solve path, this allows branch-basin switches around crossover-sensitive windows, producing artificial large discontinuities.

### Fix applied

File: `scripts/relaxtime/run_mott_phase_scan.jl`

- Added per-xi equilibrium continuation state in scan loop:
  - initialize `equilibrium_seed_state = nothing` per xi;
  - each point calls workflow with
    `seed_state=(equilibrium_seed_state === nothing ? Main.Models.HADRON_SEED_5 : equilibrium_seed_state)`;
  - after successful solve, update
    `equilibrium_seed_state = collect(res.equilibrium.x_state)`.

This keeps script-level behavior aligned with continuation expectation and suppresses nonphysical branch hopping in sequential scans.

### Post-fix quick verification

- Nonmixed fine run (`p_num=12,t_num=6,max_iter=40,T_step=5`) after fix:
  - `M_pi` jump stats return to reference-decimated level (no giant 3+ flips).
- Mixed fine run with same settings after fix:
  - pi/K jump stats also match reference-like profile;
  - remaining mixed-channel jumps are concentrated in known sensitive windows (mainly sigma/eta widths), indicating model/solver sensitivity rather than scan plumbing failure.

## Update: eta' residual edge jump final review and current baseline decision

### Final review scope

- User requested a final check on the remaining right-edge eta' jump (xi=-0.3) under the best current branch (`abctune2`).
- Goal: verify whether this can still be removed by micro-tuning, or is structural under current single-meson solve architecture.

### Final evidence

Target branch:

- `data/outputs/results/relaxtime/mott_phase/mixed_check_fine_eqcont_abctune2/mott_phase_scan_mixed_fine_eqcont_abctune2.csv`

At the problematic edge:

- `(xi=-0.3, T=255)`:
  - `M_eta_prime≈2.6649`, `Gamma_eta_prime≈1.5251`
  - `residual≈0.0501`, `root_quality=bad`
  - `governance_selection_reason=no_candidate_passed_constraints`
  - `second_pass_triggered=true`
- `(xi=-0.3, T=260)`:
  - `M_eta_prime≈2.9272`, `Gamma_eta_prime≈1.2804`
  - `residual≈8.8e-9`, `root_quality=good`

Interpretation:

- The edge jump is a transition from a non-qualified local compromise root (`T=255`) to a high-quality root (`T=260`), not a script wiring error.

### Micro-tuning exhaustion summary

Tried (incremental):

1. A+B+C tuning:
- continuity penalty, second-pass expansion, scan diagnostics columns.

2. eta'-targeted enhancements:
- second-pass eta' gamma fine grid (`0.8~2.8`),
- residual-adaptive continuity weight,
- continuity guard attempts.

3. joint-pair objective experiments (`jointpair`, `jointpair2`):
- coupled eta/eta' selection objective reduced some gamma jumps,
- but induced mass-vs-width tradeoff (Pareto), no strict global dominance over `abctune2`.

Pareto evidence at `(xi=-0.3, T=255)` candidate landscape:

- for continuity distance `<=0.25`, best residual still `~0.063` (far above good threshold),
- to reduce residual near `~0.03`, continuity distance rises to around `~1.0` scale.

Conclusion:

- Under current single-point mixed-root framework, this edge eta' jump is a structural residual in feasible set geometry.
- Micro-tuning can shift tradeoff, but cannot fully remove jump without sacrificing other criteria.

### Baseline decision (current)

- Keep `abctune2` as current most balanced operational baseline.
- Mark `(xi=-0.3, T≈255->260)` eta' edge jump as a known structural-sensitive point.
- For full removal, next stage should use a true coupled eta/eta' pair solver (not further local micro-tuning).

## Next research direction requested by user

User hypothesis to investigate:

- Phenomenon currently appears concentrated at `xi=-0.3`.
- Need to test whether this is primarily numerical (solver/path artifact) or potentially physical/non-accidental behavior.

Immediate operational probe requested:

- Extend to lower anisotropy (e.g. `xi=-0.4`) and re-check eta' vs T curve for analogous edge jump behavior.

Planned analysis framing:

1. Generate scan at `xi=-0.4` with same settings as `abctune2` profile.
2. Compare eta' jump metrics at high-T edge vs `xi=-0.3`.
3. Inspect diagnostics (`root_quality`, `governance_selection_reason`, second-pass fields).
4. If similar bad->good edge transition repeats, likely numerical-feasible-set issue;
   if pattern shifts systematically with xi, may indicate model-physical sensitivity requiring theory-side interpretation.

## Update: xi=-0.4 exploratory probe (for numerical vs physical interpretation)

### Probe setup

- Ran mixed scan with `xi_list = [-0.4, -0.3]`, `T=120:5:260`, `p_num=8`, `t_num=4`, `max_iter=25`, mixed enabled.
- Artifacts:
  - CSV: `data/outputs/results/relaxtime/mott_phase/mixed_xi_neg04_probe/mott_phase_scan_mixed_xi_neg04_probe.csv`
  - derived: `data/outputs/results/relaxtime/mott_phase/mixed_xi_neg04_probe/mott_phase_scan_mixed_xi_neg04_probe_derived.csv`
  - figures: `data/outputs/results/relaxtime/mott_phase/mixed_xi_neg04_probe/figures`

### Observations (eta' focused)

At threshold `|Δ|>=0.25` on eta' observables:

- `xi=-0.3`: 3 jump events (current branch)
  - includes one remaining high-T bad point around `T=255` (`root_quality=bad`, `no_candidate_passed_constraints`, second-pass triggered)
- `xi=-0.4`: 6 jump events
  - two early-T mass/gamma jumps near `T≈180-185` are still `root_quality=good`
  - additional high-T gamma jumps near `T≈250-260`, with mixed quality (`degraded/good`) and one joint objective selection event

### Interim interpretation

- Jump-like behavior is not unique to `xi=-0.3`; lowering to `xi=-0.4` does not remove it and can increase event count in this coarse metric.
- This weakens a “single-point accidental numeric glitch only at -0.3” interpretation.
- Current evidence is more consistent with xi-dependent sensitive regions (model + numerical path coupling), not a one-off artifact.

### Important caveat

- This probe was run on the current working branch where joint-pair mixed governance experiments are present.
- Therefore, direct quantitative comparison to historical `abctune2` baselines should be treated as directional, not final apples-to-apples, unless rerun in a clean `abctune2`-only algorithm state.

## Update: wide-xi + higher-T hypothesis check (window shift test)

User hypothesis:

- A jump window may exist for all xi and shift left/right with xi, making right-edge visibility stronger for sufficiently negative xi.

### Probe executed

- Scan range expanded:
  - `xi = [-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, 0.0]`
  - `T = 120:5:340 MeV`
- Output:
  - CSV: `data/outputs/results/relaxtime/mott_phase/mixed_xi_window_probe_wideT/mott_phase_scan_mixed_xi_window_probe_wideT.csv`
  - derived: `data/outputs/results/relaxtime/mott_phase/mixed_xi_window_probe_wideT/mott_phase_scan_mixed_xi_window_probe_wideT_derived.csv`
  - figures: `data/outputs/results/relaxtime/mott_phase/mixed_xi_window_probe_wideT/figures`
- Added analysis script:
  - `scripts/analysis/relaxtime/meson_mixed/analyze_eta_prime_window_shift.jl`

### Main findings (eta' jump threshold |Δ|>=0.25)

- Jump events are present across broad xi range, not only at `xi=-0.3`.
- High-T jump-center metric (Gamma_eta_prime weighted center, T>=220) is **non-monotonic** in xi under current branch state.
- Therefore, current evidence does **not** support a simple monotonic “as xi decreases, window shifts left” law.
- But it does support a broader statement: xi-dependent jump windows do appear across multiple xi slices.

### Interpretation status

- Current result favors “xi-sensitive mixed-channel landscape + numerical path coupling” over a single isolated accidental point.
- Monotonic shift hypothesis remains unproven in this run; likely requires:
  1) strict apples-to-apples rerun on frozen `abctune2` algorithm only,
  2) finer T grid around detected windows,
  3) optionally longer T upper bound for selected xi slices to avoid edge-truncation bias.

## Update: re-check on “window left-shift” and root-cause-level interpretation

User feedback after manual figure inspection:

- visually, a broad “left-shift tendency” of jump-window center with decreasing xi seems plausible.

Data-side re-check (same wide probe CSV) indicates:

- strict monotonic center metric remains non-monotonic globally;
- however, several onset-style indicators do show a directional tendency:
  - early jump-onset temperatures and non-good onset temperatures are generally lower for more negative xi in multiple sub-ranges,
  - with local reversals and overlaps (hence not strictly monotone law).

Working synthesis:

- “left-shift tendency” is a reasonable qualitative description,
- but should be stated as a trend with local exceptions, not an exact monotonic rule.

### Why eta' jumps while eta appears much smoother (root-cause-level view)

Observed asymmetry:

- eta channel remains mostly stable (very low non-good/no-candidate incidence, no large mass jumps under chosen threshold);
- eta' shows frequent non-good and no-candidate states and dominant jump signatures.

Mechanistic interpretation (current model + solver context):

1. Branch geometry asymmetry in mixed sector
- eta is tied to minus branch; eta' to plus branch.
- Near sensitive windows, plus branch objective surface is flatter/more competing, so eta' root qualification is harder.

2. Feasible-set conflict at edge points
- For key bad points (e.g. xi=-0.3, T=255), candidate landscape shows no simultaneously low-residual + high-continuity solution in local neighborhood.
- Thus solver must pick compromise points, producing visible jumps when next-T point returns to good root.

3. Width coupling amplifies eta' visibility
- eta' jump is strongest in Gamma channel; width enters mixed-score and continuity selection simultaneously.
- This makes eta' more sensitive to branch switching and path decisions than eta mass line.

4. xi acts as landscape shaper
- Lower xi changes equilibrium background and mixed matrix conditioning, moving/reshaping sensitive windows.
- This explains trend-like migration with local non-monotonic behavior.

Bottom line:

- eta' jump is not a plotting artifact and not merely random numerical noise.
- It is an emergent result of mixed-branch plus-channel sensitivity + solver-governance tradeoff under xi-dependent landscape deformation.
- eta smoothness vs eta' jump is therefore expected under current architecture.

## Update: direct answers on convergence status and continuity-failure root cause

User questions:

1) Are values inside jump window unconverged?
2) Where is the converged solution closest to continuity?
3) What is the upstream/root cause of continuity solve suddenly failing?

### 1) Convergence status in jump window (eta')

For the representative residual edge point (`xi=-0.3`, `T=255`, current baseline path):

- selected eta' has `root_quality=bad`, `residual≈0.0501`, `no_candidate_passed_constraints`, second-pass triggered.
- at next point (`T=260`) selected eta' is `root_quality=good`, `residual≈8.8e-9`.

Thus yes: inside this jump window, selected eta' value can be non-converged/under-qualified, and jump arises when next-T returns to a good root.

### 2) Closest converged solution to continuity (at T=255)

A dense local seed/method scan (`m0∈[2.2,3.3]`, `g0∈[0,3.0]`, both newton+trust-region, iterations=120):

- scanned candidates: `1426`
- converged candidates: `0`
- good candidates (`resid<=1e-6`): `0`

So at this point, **there is no converged candidate at all** in tested local basin. The nearest-to-continuity converged solution is therefore undefined in current tested domain.

### 3) Why continuity suddenly fails (upstream diagnosis)

Observed upstream changes from `T=255 -> 260` at `xi=-0.3`:

- `m_u`: `0.05413 -> 0.05044`
- `m_s`: `1.33591 -> 1.27261`
- `Phi`: `0.78378 -> 0.79250`
- anisotropic fields `A_u/A_s` and effective couplings `G_u/G_s` shift accordingly.

Implication:

- mixed plus-branch landscape deforms enough between adjacent T points that the continuity-near region remains high residual at T=255,
  while a good-root basin appears/opens at T=260.
- This is consistent with “feasible-set gap” at T=255 rather than simple seed insufficiency.

Additional local evidence:

- continuity-near neighborhood floor at T=255 remains relatively high in plus-score,
- lower-score candidates exist only at large continuity distance,
- matching the residual-vs-continuity conflict already identified in Pareto diagnostics.

### Residual-surface visualization support (human-review oriented)

To support direct visual audit of the two-equation residual system, added export+plot tooling:

- exporter: `scripts/analysis/relaxtime/meson_mixed/export_mixed_residual_surface.jl`
  - computes complex mixed residual `F(m, gamma)` on grid
  - outputs `Re(F)`, `Im(F)`, `||F||`, `log10||F||` surface CSV + mass/gamma slices
- plotter: `scripts/analysis/relaxtime/meson_mixed/plot_mixed_residual_surface.py`
  - outputs heatmap (`log10||F||`) and 1D slices for quick visual inspection

Generated example datasets/figures:

1. jump-window point (`T=255`, `xi=-0.3`, eta' plus)
- `data/outputs/results/relaxtime/mott_phase/residual_surface/edge255_xi-0p3_etaplus/residual_surface.csv`
- `data/outputs/results/relaxtime/mott_phase/residual_surface/edge255_xi-0p3_etaplus/figures/residual_surface_heatmap.png`

2. neighboring good point (`T=260`, `xi=-0.3`, eta' plus)
- `data/outputs/results/relaxtime/mott_phase/residual_surface/edge260_xi-0p3_etaplus/residual_surface.csv`
- `data/outputs/results/relaxtime/mott_phase/residual_surface/edge260_xi-0p3_etaplus/figures/residual_surface_heatmap.png`

These plots are intended for human side-by-side inspection of basin-structure changes across the edge transition.

### Iteration trajectory overlays and root-existence final probe

Added trajectory overlays (from continuity-like seeds):

- `T=255`: `data/outputs/results/relaxtime/mott_phase/residual_surface/edge255_xi-0p3_etaplus/figures_components/trajectory_overlay.png`
- `T=260`: `data/outputs/results/relaxtime/mott_phase/residual_surface/edge260_xi-0p3_etaplus/figures_components/trajectory_overlay.png`

Observation:

- at `T=255`, trust-region trajectory does not enter the visible low-residual basin;
- at `T=260`, trajectory rapidly enters and converges to the stable basin.

To answer “is this only step-size/precision too low?” directly, ran high-precision multi-seed probe:

- script (historical temporary probe, cleaned from repo): `tmp/final_root_existence_probe.jl`
- dense seeds around suspected basin (`m0: 2.60..3.00, g0: 1.80..2.80, step 0.02`, two methods)
- strict solver settings (`iterations=1200`, `ftol=xtol=1e-14`)

Results:

1) `T=255, xi=-0.3`
- solved attempts: `2142`
- converged: `170`
- but `residual<=1e-6`: `0`
- best residual only `~1.434e-2`

2) `T=260, xi=-0.3`
- solved attempts: `2142`
- converged: `1354`
- `residual<=1e-8`: `1353`
- dominant root cluster: `(m≈2.877647, g≈2.654725)`

Interpretation:

- This is not explained by a simple “coarse step / low precision missed root” hypothesis.
- At `T=255` in current formulation, no true high-quality eta' plus root is found in the examined basin even under aggressive precision and seed density.
- At `T=260`, the good basin is robust and wide-capture.

### Final enhanced visualization: direct-eval points vs interpolated surface

To address interpolation-illusion risk explicitly, added direct-evaluation overlay tools:

- exporter: `scripts/analysis/relaxtime/meson_mixed/export_direct_eval_scatter.jl`
  - computes residual norm directly from equation evaluations on dense `(m,gamma)` lattice
- overlay plotter: `scripts/analysis/relaxtime/meson_mixed/plot_direct_vs_surface_overlay.py`
  - overlays direct low-residual points (`<=1e-2`, `<=1e-3`, `<=1e-4`) on top of interpolated surface background

Generated overlays:

- `T=255`: `data/outputs/results/relaxtime/mott_phase/residual_surface/edge255_xi-0p3_etaplus/figures_components/direct_vs_surface_overlay.png`
- `T=260`: `data/outputs/results/relaxtime/mott_phase/residual_surface/edge260_xi-0p3_etaplus/figures_components/direct_vs_surface_overlay.png`

Dense direct-eval count summary:

- `T=255`: min direct residual `~1.474e-2`; counts `<=1e-2: 0`, `<=1e-3: 0`, `<=1e-4: 0`
- `T=260`: min direct residual `~4.416e-4`; counts `<=1e-2: 484`, `<=1e-3: 6`, `<=1e-4: 0`

Interpretation:

- At `T=255`, the visually deeper basin is not supported by direct low-residual equation evaluations in this grid domain.
- At `T=260`, direct low-residual points are present and clustered, matching robust solver convergence behavior.

## Update: minimal nudge restart scan evaluation (`mixed_check_fine_nudged`)

To answer whether the newly added minimal nudge restart can move trajectories into lower-residual basins,
we compared:

- baseline: `mixed_check_fine_eqcont_abctune2`
- nudged: `mixed_check_fine_nudged`

under the same grid (`T=120..260`, `step=5`, `xi=[-0.3,-0.15,0.0,0.15,0.3]`, `muB=0`).

### Data-level comparison summary (eta-prime channel)

- common compared rows: `145`
- max absolute difference:
  - `|ΔM_eta_prime| = 5.5801e-2`
  - `|ΔGamma_eta_prime| = 2.0193e-1`
  - `|Δresidual_eta_prime| = 9.2335e-3`
- residual change counts:
  - improved: `45`
  - worsened: `12`
  - unchanged: `88`

Root-quality distribution:

- baseline: `good=87, degraded=55, bad=3`
- nudged: `good=87, degraded=54, bad=4`

Interpretation:

- nudge restart changes are localized (many rows unchanged), but not uniformly beneficial;
- global quality does **not** improve (same `good` count, one extra `bad`).

### Key edge window check (xi=-0.3, T=245..260)

Critical rows around the original jump:

- `T=245`: residual `4.4139e-2 -> 4.9785e-2` (worse)
- `T=250`: residual `5.3389e-2 -> 5.8127e-2` (worse)
- `T=255`: residual `5.0134e-2 -> 5.9367e-2` (worse)
- `T=260`: remains high-quality tiny-residual root (`~8.82e-9`, unchanged quality)

So in the exact problematic boundary region, this minimal nudge did **not** pull the track into a better basin;
it increased pre-jump residuals and slightly enlarged `M_eta_prime` jump amplitude.

### Jump metric comparison (adjacent-step max)

For `xi=-0.3`:

- `max |ΔT M_eta_prime|`: `0.2623 -> 0.2901` (worse)
- `max |ΔT Gamma_eta_prime|`: `0.4020 -> 0.3114` (better)

The mass jump worsened while width jump improved; net continuity gain is mixed and does not resolve the original eta-prime mass jump concern.

### Side-effect scan (other channels)

Across all compared rows:

- nonmixed channels (`pi/K`) unchanged at machine precision;
- mixed channels other than eta-prime show only tiny numerical perturbations (`~1e-8` scale for masses, `~1e-8..1e-7` for widths).

This indicates the nudge logic is localized and does not destabilize broad pipeline outputs.

### Current conclusion

- Minimal nudge restart is a low-risk local intervention but **insufficient** for the eta-prime edge jump objective.
- It can improve some xi/T points, yet fails at the target boundary (`xi=-0.3`, `T~245..255`) and may worsen mass-jump continuity there.
- Next effective step likely requires stronger basin-selection policy (e.g., objective that explicitly balances residual and continuation distance for eta/eta' pair) rather than fixed one-shot positive nudge vectors.

## Update: ROI rescue experiment (Plan A) for T255 target basin

User target basin for eta-prime plus branch at `T=255, xi=-0.3`:

- preferred region: `(m, gamma) in ([2.82,2.92], [2.30,2.50])`
- expected local minimum scale from direct scan: residual `~1.47e-2` near `(2.862, 2.412)`

Implemented workflow-level ROI rescue (no solver-core changes):

- file: `src/models/workflows/MesonMassWorkflow.jl`
- trigger gate (mixed/no-pass path):
  - meson is `:eta_prime`
  - selected residual remains high (`>= 0.015`)
- action:
  - launch ROI lattice seeds over `[2.82,2.92]x[2.30,2.50]` (`5x5`, both TR/Newton)
  - build candidates as `seed_source=:roi_rescue`
  - accept replacement only when objective improves and continuity/residual guards pass
  - mark selection reason as `:roi_rescue_objective_min` when replacement is chosen

TDD evidence:

- updated unit tests: `tests/unit/models/test_meson_mass_workflow.jl`
- red phase observed first (missing diagnostics fields), then implementation, then green:
  - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_meson_mass_workflow.jl"; include("tests/unit/runtests.jl")'`
  - final: `25/25` pass

### Scan-level evaluation (same coarse grid as baseline)

Run:

- `scripts/relaxtime/run_mott_phase_scan.jl`
  - `--tstep 5 --p-num 8 --t-num 4 --max-iter 25 --include-mixed`
  - outdir: `data/outputs/results/relaxtime/mott_phase/mixed_check_fine_roiA`

Key result vs baseline `abctune2` (`xi=-0.3`, edge window):

- `T=245`: residual `0.0441 -> 0.0498` (worse)
- `T=250`: residual `0.0534 -> 0.0581` (worse, reason became `roi_rescue_objective_min`)
- `T=255`: residual `0.0501 -> 0.0594` (worse)
- `T=260`: stays good tiny residual (`~8.82e-9`)

Jump metrics (`xi=-0.3`):

- `max |ΔT M_eta_prime|`: `0.2623 -> 0.2901` (worse)
- `max |ΔT Gamma_eta_prime|`: `0.4020 -> 0.3114` (better)

### Interpretation

- In isolated eta-prime solve, ROI rescue can be attempted and can move candidates.
- In full mixed workflow path, eta/eta' pairing and downstream pair-objective selection can override local ROI benefit at the final reported track.
- Therefore current Plan A implementation is still insufficient to force the final T255 trajectory endpoint into the requested ROI basin under full scan governance.

### Practical next step (if continuing)

To hit the exact user goal, ROI awareness must be lifted to **pair selection level** (eta + eta_prime jointly), not only single-channel candidate rescue.

- candidate direction: add ROI-aware term into joint pair objective for eta-prime side (soft preference, not hard clamp), with guard that eta branch continuity does not collapse.

## Update: Plan B completed (pair-level ROI objective)

Implemented Plan B by lifting ROI preference into joint eta/eta-prime pair objective:

- file: `src/models/workflows/MesonMassWorkflow.jl`
- changes:
  - inject ROI seeds into joint candidate pool for `eta_prime` in target window (`xi<=-0.25`, `T=240..260 MeV`, residual high)
  - add ROI distance penalty term directly in pair objective (`pair` level, not single-channel post-select)
  - expose pair selection reason `:pair_min_roi` when ROI-aware mode is active

TDD check:

- test file: `tests/unit/models/test_meson_mass_workflow.jl`
- added continuation-context assertion (`245 -> 250 -> 255`) requiring T255 eta-prime enters ROI neighborhood with low residual
- run result: `29/29` pass

### Scan evidence (`mixed_check_fine_roiB`, same coarse settings)

Compared with baseline `abctune2` at `xi=-0.3`:

- `T=250`:
  - residual: `0.05339 -> 0.03125` (improved)
  - `(m,gamma): (2.7099,1.1231) -> (2.9725,2.4142)`
- `T=255`:
  - residual: `0.05013 -> 0.01435` (strongly improved; near expected local minimum scale)
  - `(m,gamma): (2.6649,1.5251) -> (2.9069,2.5115)`

Target basin hit status at `T=255`:

- mass in target box: yes (`2.82..2.92`)
- gamma in target box: marginally above upper edge (`2.5115` vs target max `2.50`, overshoot `~0.0115`)

Jump metrics (`xi=-0.3`):

- `max |ΔT M_eta_prime|`: `0.2623 -> 0.2466` (improved)
- `max |ΔT Gamma_eta_prime|`: `0.4020 -> 1.3056` (degraded; large width jump introduced around high-T edge)

### Interpretation

- Plan B successfully drives T255 endpoint into the intended low-residual neighborhood (core user objective substantially met on residual and mass).
- Remaining issue is gamma continuity/stability side effect (large width jump), so this is not yet a clean global replacement policy.
- Next refinement should keep Plan B as base and add a width-jump guard (or bounded gamma corridor) to prevent over-amplified `Gamma_eta_prime` transitions.
