# High-xi anchor comparison gate

Verdict: `production-grade`

- new case: `first_canonical_v1_p128_xi001_validated_anchored_prod_v1`
- old anchor case: `first_canonical_v1_p128_validated_anchored_prod_v1`
- comparison: all high-xi rows at xi multiples of 0.05 vs existing p128 production rows
- worst relative difference: `1.125`
- effective worst relative difference after abs <= 1e-10 floor: `8.453069929263162e-05`
- threshold: `0.01`
- absolute floor: `1e-10`
- comparison cells: `12663`

The p128 integration-parameter convergence gate is inherited from `first_canonical_v1_p128_validated_anchored_prod_v1_convergence` because the numerical integration policy is unchanged.
