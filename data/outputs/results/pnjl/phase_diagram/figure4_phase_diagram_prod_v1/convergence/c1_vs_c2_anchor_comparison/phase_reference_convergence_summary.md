# PNJL Phase Reference Convergence Summary

- candidate: `c1_p24t8_full_grid`
- reference: `c2_p32t12_anchors`
- verdict: `needs_manual_review`

## Artifact Inventory

| side | artifact | rows | xi_count | naninf | converged_false |
|---|---:|---:|---:|---:|---:|
| candidate | boundary | 72 | 5 | 0 | None |
| candidate | cep | 5 | 5 | 0 | None |
| candidate | spinodals | 72 | 5 | 0 | None |
| candidate | crossover | 80 | 5 | 0 | 0 |
| reference | boundary | 72 | 5 | 0 | None |
| reference | cep | 5 | 5 | 0 | None |
| reference | spinodals | 72 | 5 | 0 | None |
| reference | crossover | 80 | 5 | 0 | 0 |

## Max Differences

| artifact | metric | n | max_abs | max_rel | max_abs_key | max_rel_key |
|---|---|---:|---:|---:|---|---|
| boundary | mu_transition_MeV | 72 | 0.26135915901812723 | 0.0007180345728965059 | 0.25|60.0 | 0.25|60.0 |
| boundary | rho_hadron | 72 | 0.005763311651990066 | 0.019988297936973612 | 0.5|60.0 | 0.25|60.0 |
| boundary | rho_quark | 72 | 0.039975326268347544 | 0.01601860207577932 | 0.0|60.0 | 0.0|60.0 |
| cep | T_CEP_MeV | 5 | 0.029296875 | 0.0002737850787132101 | 0.5 | 0.5 |
| cep | muB_CEP_MeV | 5 | 0.07596633469188419 | 7.26933482813994e-05 | 0.5 | 0.5 |
| cep | muq_CEP_MeV | 5 | 0.025322111563980343 | 7.269334828145378e-05 | 0.5 | 0.5 |
| crossover | T_crossover_MeV | 80 | 0.030414466061785106 | 0.00028423992503987113 | 0.5|16 | 0.5|16 |
| crossover | derivative | 80 | 1213.8497260563827 | 1.0691971821584236 | 0.5|16 | 0.5|16 |
| crossover | mu_MeV | 80 | 0.025322111563980343 | 7.269334828169857e-05 | 0.5|16 | 0.5|4 |
| crossover | rho | 80 | 0.1691914700899224 | 0.10131905834929714 | 0.25|16 | 0.25|16 |
| spinodals | mu_spinodal_hadron_MeV | 72 | 0.07164561553776139 | 0.00019038760331705967 | 0.5|60.0 | 0.5|60.0 |
| spinodals | mu_spinodal_quark_MeV | 72 | 1.4341450249249874 | 0.004188879452428898 | 0.0|60.0 | 0.0|60.0 |
| spinodals | rho_spinodal_hadron | 72 | 0.005876862020160578 | 0.010859701008807465 | 0.25|60.0 | 0.25|60.0 |
| spinodals | rho_spinodal_quark | 72 | 0.04928972563382872 | 0.027819350815090606 | 0.5|60.0 | 0.5|60.0 |
