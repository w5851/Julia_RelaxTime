# trho_asymmetric K+/pi+ Meson Density Scan v1

Verdict: `blocked`

This directory records the convergence-gate evidence for the
`FixedAsymmetricRho` density-constrained `K+ / pi+` meson-density path.
No production-grade formal data were generated.

## Scope

- Path strategy: `trho_asymmetric`
- Equilibrium source: `FixedAsymmetricRho`
- Meson chemical profile: `asymmetric_kplus_over_piplus_signed`
- Constraint targets: `rho_u / rho_d = 0.876`, `rho_s = 0 fm^-3`
- Grid attempted: `T = 120:10:220 MeV`, `rho/rho0 = 0.05:0.05:1.00`
- Regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`

## Evidence

The initial remote `convergence_low` run completed successfully as a workflow run, but
the physics/numerical domain guard failed for most non-stable regimes:

- rows: `880`
- path points: `220`
- all-regime-ok points: `13`
- status counts: `ok=245`, `failed=221`, `unsafe_bose_domain=414`
- stable failures: `14`
- strict BW failures: `207`
- phase-shift unsafe rows: `414`

The failures are not node-count convergence failures. They are caused by the
positive charged-kaon chemical potential in the selected `K+ / pi+` asymmetric
profile. Several stable rows already hit `mass <= mu_K`; the BW and phase-shift
energy integrals also include `omega <= mu_K` under the default
`strict_normal_domain` policy.

Follow-up remote convergence attempts on commit `4c6c1542` used the explicit
phase-shift-only policy `density_policy=x_min_cut`.

- `bose_x_min=1e-6`: BU/GBU phase-shift rows became runnable, but `mid -> high`
  K-related quantities still showed order-50% maximum relative differences.
- `bose_x_min=1e-2`: cutoff sensitivity improved, but `mid -> high` K-related
  maximum relative differences remained about 11%--13%, and high/custom
  comparisons exposed different upstream `FixedAsymmetricRho` equilibrium
  branches at the same `(T, rho_target)`.

The current blocker is therefore upstream equilibrium branch stability in the
combined scan path. Production remains blocked until the branch-selection /
continuity policy is fixed and a fresh convergence gate passes.

## Files

- `convergence/convergence_low/combined_meson_density_scan.csv`
- `convergence/convergence_low/README.md`
- `convergence/convergence_low/remote_run_manifest.json`
- `convergence/convergence_low/run.command.txt`
- `convergence/convergence_low/run.stdout.log`
- `convergence/convergence_low/run.stderr.log`
- `convergence/convergence_summary.json`
- `PRODUCTION_AUDIT.md`
- `manifest.json`

Figure-side convergence evidence:

`data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/`

## Next Decision

To continue toward production-grade output, the physics scope must be changed
or clarified. Candidate follow-ups are:

- rerun convergence after aligning `trho_asymmetric` scan order with the
  `TrhoScan` density-constrained continuity policy;
- switch to the negative charged channel profile, where `mu_K` is not positive;
- restrict the grid to a domain where all required Bose-domain constraints hold;
- explicitly downgrade to a diagnostic-only policy if production-grade
  convergence is not physically/numerically justified.

Those are material scope changes and were not applied in this blocked run.
