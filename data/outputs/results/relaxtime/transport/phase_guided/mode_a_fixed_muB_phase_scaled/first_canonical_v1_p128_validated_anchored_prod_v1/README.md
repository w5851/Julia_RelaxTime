# first_canonical_v1_p128_validated_anchored_prod_v1

Phase-guided transport production-grade case for `mode_a_fixed_muB_phase_scaled`.

## Verdict

`production-grade`

This result is the promoted p128 sharded GitHub Actions production artifact. The action artifacts were first treated as candidates, then promoted after a full-grid `p104 -> p128` convergence gate passed with no matched point above 1% relative difference.

## Scope

- mode: `mode_a_fixed_muB_phase_scaled`
- case name: `first_canonical_v1_p128_validated_anchored_prod_v1`
- summary: fixed `mu_B`, phase-reference `T/T_phase` bands, continuous `xi` scan inside each band
- grid: `mu_B = 0, 450, 900 MeV`, `alpha_T = 1.0, 1.1, 1.2`, `xi = -0.50:0.05:0.50`
- planned and completed points: `189`
- compute bulk viscosity (`zeta`): `true`

## Numerical Policy

- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `tau_p_nodes = 128`
- `tau_angle_nodes = 20`
- `tau_phi_nodes = 36`
- `tau_n_sigma_points = 24`
- `sigma_grid_n = 560`
- `channel_diagnostics = true`

The `validated_anchored` sigma-cache policy is part of the production numerical policy for this result. It reduces the threshold-subtraction / sigma-cache sensitivity diagnosed in the earlier `first_canonical_v1` result while preserving the physical `match_thermo` propagator-xi policy.

## Key Files

- `phase_guided_transport_scan.csv`: promoted transport scan output
- `sampling_plan.csv`: full phase-guided sampling plan
- `channel_diagnostics.csv`: per-channel scattering-rate diagnostics
- `failed_points.csv`: failure ledger; contains `0` data rows
- `effective_config.json`: repository-imported effective configuration snapshot
- `manifest.json`: production manifest with source runs, hashes, and convergence evidence
- `PRODUCTION_AUDIT.md`: human-readable production audit

Shared convergence evidence is stored in:

- `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence/`

Canonical figures are stored in:

- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_validated_anchored_prod_v1/`

## Convergence Summary

- adjacent comparison: `104/18/30/20/440 -> 128/20/36/24/560`
- matched points: `378` across mode A and mode B
- overall max relative difference: `0.009426486578504083`
- core observable max relative difference: `0.006176952435828459`
- points above 1% relative difference: `0`
- points above 0.5% relative difference: `12`

Worst p104-to-p128 point: `mode_a`, `muB=450 MeV`, `alpha_T=1.1`, `xi=0.35`, `tau_ubar/tau_dbar`, relative difference `0.9426486578504083%`.

## Supersession Note

This case supersedes `first_canonical_v1` for production-grade interpretation of local `xi` structures. The old result is retained as a historical comparison artifact, but it lacks explicit high-precision tau/sigma integration parameters, `validated_anchored` sigma-cache treatment, channel diagnostics, and a full-grid convergence gate. Local non-monotonic structures or precision-sensitive claims should therefore use this p128 production case instead of the old `first_canonical_v1` case.

The old result is not removed in this repository update. Keeping it avoids losing provenance and makes the numerical improvement auditable.
