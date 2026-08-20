# Author Review Decisions: Historical Figure Assets

Review date: 2026-08-16

Source inventory:

- docs/analysis/figure_asset_registry_v1/asset_registry.json
- docs/analysis/figure_asset_registry_v1/cleanup_candidates.csv

Registry snapshot used for this review:

- generated_at_utc: `2026-08-16T04:42:45.270159+00:00`
- source_git_commit: `bf249323641fd9976b303f8651fa43771e1c52e9`
- asset_count: `446`
- case_group_count: `82`
- review_group_count: `82`
- visible_untracked_excluded: `12`
- missing_tracked_assets: `0`
- asset_registry.json SHA-256: `B8C33851CEC2E44A372067FA8F29D368E9CA55FCDC2300EA1E72857B46DDA80D`
- cleanup_candidates.csv SHA-256: `F8D773188639328118A897A0C4B282F09A38AAAE993D994EC3EE3E5E700EEA64`

Status: author review complete. PR A remains a read-only inventory and dry-run. The decisions below are the approved input for a separate cleanup/retirement PR; no deletion, move, or rename is performed by this record.

## Guardrails

- Existing numerical CSV/JSON results, solver behavior, Maxwell/C2/reference semantics, and transport numerical semantics are unchanged.
- Untracked C1/C2/pilot assets remain untouched.
- A move must preserve the source result data, add a traceable README, update active references, and record source/destination hashes.
- A delete must use an explicit absolute-path allowlist, recheck tracked status and SHA-256 immediately before execution, and produce a relocation/deletion manifest.
- The original registry and cleanup_candidates.csv remain the PR A dry-run snapshot. They are not rewritten into executable delete/move instructions.

## Group 1

| Review group | Author decision |
| --- | --- |
| data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1 | Keep PNG as legacy because it is an early low-precision result; delete the PDF. No external citation was reported. |
| data/outputs/figures/pnjl/phase_diagrams | Delete the complete directory. |
| data/outputs/figures/relaxtime/cross_section/T100 | Keep as legacy. The plot is useful, but the data and plotting are early/low precision. |
| data/outputs/figures/relaxtime/cross_section/T200_muB0_xi0.0 | Keep the reviewed PNG assets as legacy. |
| data/outputs/figures/relaxtime/cross_section/T250 | Keep as legacy. |

The following T200 subdirectories have the same keep-as-legacy decision:

- process_ssbar_to_uubar
- sqrt_s_group1
- sqrt_s_group1/ylim_0_10
- sqrt_s_group2
- sqrt_s_group2/ylim_0_10
- sqrt_s_group3
- sqrt_s_group3/ylim_0_10

## Group 2

| Review group | Author decision |
| --- | --- |
| data/outputs/figures/relaxtime/fixed_temperature_xi_scan_muB0 | Delete the complete directory tree because newer results exist and the old results have substantial precision problems. |
| data/outputs/figures/relaxtime/gap_transport_by_xi_muB0p0 | Delete. |
| data/outputs/figures/relaxtime/gap_transport_by_xi_muB800 | Delete. |

## Group 3

| Review group | Author decision |
| --- | --- |
| data/outputs/figures/relaxtime/literature | Keep the comparison evidence but relocate the five image assets to docs/analysis/relaxtime/literature_comparison/ and add a context README. Leave validation target data untouched. |
| data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan | Keep. The manifest and source CSV were introduced together in commit 2f5ba06b on 2026-06-05; no later update was found. The author’s conditional legacy downgrade is therefore not triggered. This old manifest is not promoted to strict. |
| data/outputs/figures/relaxtime/meson_density/fig3_like_tmu_equal_flavor_folded_highres | Keep as legacy. |
| data/outputs/figures/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100 | Relocate the overlay/residual figure evidence to docs/analysis/relaxtime/meson_density/freezeout_kminus_piminus_mu_pi_100_analysis/ and add a README describing inputs, comparison meaning, metrics, and context. Keep source result CSVs and generator provenance. |
| data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1 | Keep. |
| data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low | Delete. |
| data/outputs/figures/relaxtime/meson_mass/path_scan/freezeout | Keep only the strict line-first v2 case listed below; delete the other three child variants. |
| data/outputs/figures/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0__plotv1__strict__line_first_v2 | Keep. |
| data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a | Keep as a current paper formal figure or formal candidate. |
| data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab | Keep as a current paper formal figure or formal candidate. |

The other freezeout child variants approved for deletion are:

- default_baseline_freezeout_xi0__plotv1__strict
- default_baseline_freezeout_xi0__plotv1__strict__line_first_v1
- default_baseline_freezeout_xi0__plotv1__strict__pilot_v2

## Group 4

| Review group | Author decision |
| --- | --- |
| data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_b | Keep as a current paper formal figure or formal candidate. |
| data/outputs/figures/relaxtime/paper_p1_mott_phase_isentropic/production_20260531 | Keep as a current paper formal figure or formal candidate. |
| data/outputs/figures/relaxtime/temperature_scan_muB0_xi0 | Delete. |
| data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_validated_anchored_prod_v1 | Keep and rename to first_canonical_v1_p128_xi005_validated_anchored_prod_v1. The result README and effective configuration confirm xi=-0.50:0.05:0.50. Rename the result and figure case paths together and update path-only metadata without changing numerical CSV values. |
| data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1 | Same xi005 rename and retention decision as mode A. |
| dense v1 xi001 mode A/B cases | Keep as legacy. They use the older RS-shared-energy transport-kernel energy semantics and remain needed for historical provenance and existing analysis references. |
| dense v2 xi001 mode A/B cases | Keep as current_candidate. They use onshell_kinematic_energy in the transport kernel and are the intended future current source, subject to promotion and manuscript eligibility review. |
| data/outputs/figures/relaxtime/validation | Relocate as transport error/diagnostic analysis under docs/analysis/relaxtime/ and add context documentation. |

The v1/v2 comparison is recorded inside both v2 result directories in COMPARISON_TO_V1.md. The two cases have the same 909-point grid and nearly unchanged upstream/tau quantities, but transport semantic drift reaches about 28% for eta, 25% for sigma, and 77% to 80% for zeta. Therefore v2 is not a format-only replacement. v1 is not deleted while v2 remains current_candidate and manuscript_eligible=false.

## Groups 5-9

The author approved deletion of the complete directory:

data/outputs/figures/relaxtime/xi_smoothness_sampling

This covers review groups S001 through S042. These are early low-precision random-xi parameter-point scans; the precision issue is resolved and the random sampling workflow is no longer needed. The individual S001-S042 groups therefore do not require separate follow-up review.

## Completeness review

The author supplied a decision for every tracked review group:

- retain as legacy;
- retain as current formal/candidate;
- retain as analysis evidence and relocate with context;
- retain with explicit variant rename;
- or delete using a later explicit allowlist.

The only conditional decisions were resolved by repository evidence during review:

- combined meson-density scan: no later update found;
- freeze-out K-/pi- plot review: confirmed as manual overlay/residual analysis;
- dense transport v1/v2: confirmed as a transport energy-semantics change, not duplicate rendering.

No unresolved author-choice gate remains for PR A. Actual deletion, relocation, and renaming remain separate PR B actions.
