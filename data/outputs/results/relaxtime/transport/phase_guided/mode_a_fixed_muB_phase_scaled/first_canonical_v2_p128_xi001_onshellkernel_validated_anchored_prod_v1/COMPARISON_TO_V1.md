# Comparison To Dense v1

This note records the historical relationship between this v2 result and the retained dense v1 result. It is provenance documentation only; no numerical result file is modified by this note.

## Cases

- v1: first_canonical_v1_p128_xi001_validated_anchored_prod_v1
- v2: first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1
- mode: mode_a_fixed_muB_phase_scaled
- grid: 909 rows, xi=-0.50:0.01:0.50

## Semantic Difference

The v1 case follows the older RS-shared-energy transport-kernel convention. The v2 case explicitly uses:

- transport kernel energy: onshell_kinematic_energy, E_kin=sqrt(p^2+M^2)
- distribution energy: rs_deformed, E_dist=sqrt(p^2+M^2+xi*(p*cos(theta))^2)
- propagator xi policy: match_thermo
- sigma cache policy: validated_anchored

This is a transport numerical-semantics change, not a figure restyle or interpolation change.

## Comparison Evidence

All 909 old/new rows are matched. The unchanged upstream and relaxation-time fields have a maximum relative difference of 1.66503885983704e-4, and the xi=0 transport comparison has a maximum relative difference of 3.11013314454199e-7.

The semantic transport drift maxima are:

| Observable | Maximum relative difference |
| --- | ---: |
| eta | 0.278495551584214 |
| sigma | 0.247471876832126 |
| zeta | 0.765990249103789 |
| eta_over_s | 0.278495551562842 |
| zeta_over_s | 0.765990249103789 |
| sigma_over_T | 0.247471876832126 |

Evidence files:

- convergence/production_import_summary.json
- convergence/old_vs_new_full_grid_comparison.csv
- effective_config.json
- manifest.json

## Promotion And Retention Decision

The v2 result is production-grade numerically but remains current_candidate with manuscript_eligible=false. It is the intended source for future current transport figures after the independent promotion gate.

The v1 result remains retained as legacy because existing analysis material references it and it is required to preserve the historical semantic comparison. It must not be silently replaced or deleted while those references remain.

The corresponding v1 and v2 figure directories are separate sibling cases. PR A does not delete, move, or rename either case.
