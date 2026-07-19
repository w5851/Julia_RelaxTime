# Full-range PNJL phase-reference convergence plan

Status: `configured / diagnostic-only` until every gate below has been audited.

## Scope lock

- Case slug: `fullrange_phase_reference_adaptive_prod_v1`.
- Physics scope: PNJL first-order boundary, spinodals, CEP, and crossover for anisotropy `xi in [-0.5, 0.5]`.
- Formal positive-temperature range: `T=1--240 MeV`; this case does not represent strict `T=0`.
- Initial axes: `xi=-0.5:0.05:0.5`, `T=1:5:236` plus the explicit `240 MeV` endpoint, and `rho=0:0.05:4` in `rho/rho0`.
- Thermal policy: `rs_reduced_adaptive`; angular anisotropy enters the scalar phase thermal kernel only through the RS distribution energy. This policy is not asserted for magnetic or transport kernels.
- Execution surface: GitHub Actions workflow `.github/workflows/pnjl-dense-reference.yml`; no full production run is executed locally.
- Missing/invalid semantics: any missing required anchor, non-finite core observable, failed crossover row, unconverged grid record, or unresolved phase status blocks production promotion.
- Non-goals: no transport production, no overwrite of existing references, no canonical `boundary_dense.csv`/`cep_dense.csv` promotion, and no paper-text changes.

## Action convergence matrix

All dense runs use `crossover_only=false`, `crossover_n_mu=16`, `crossover_mu_max=450 MeV`, explicit `crossover_T_max=240 MeV`, `p_num/t_num=24/8` for provenance, and `iterations=80` unless stated otherwise.

| ID | tag | coverage | thermal controls | phase-grid controls | role |
| --- | --- | --- | --- | --- | --- |
| C0 | `fullrange_phase_conv_base_20260719` | full xi/T/rho | `rtol=1e-8`, `atol=1e-10`, `maxevals=1e7` | T/xi levels 2; default geometry gates | common candidate |
| C1 | `fullrange_phase_conv_integral_tight_20260719` | full xi/T/rho | `rtol=1e-10`, `atol=1e-12`, `maxevals=2e7` | identical to C0 | integral-tolerance gate |
| C2 | `fullrange_phase_conv_grid_tight_20260719` | full xi/T/rho | identical to C0 | T/xi levels 3; position/density/area/response gates halved; CEP width `0.05 MeV` | phase-grid gate |
| C3 | `fullrange_phase_conv_tensor_anchor_20260719` | `xi=-0.5,-0.25,0,0.25,0.5`, `T=60--240 MeV` | independent `tensor_gauss`, `p_num/t_num=48/16`, iterations 100 | C0 geometry gates, no adaptive xi | coupled phase-line oracle |
| O1 | nightly full validation at the locked source SHA | fixed states at `T=1--240 MeV`, `xi=-0.5,0,0.5` | independent split-domain 2-D tensor oracle | n/a | low-temperature kernel oracle |

After C0--C3 and O1 pass, a separate formal Action rerun will combine the C1 thermal controls and C2 phase-grid controls under tag `fullrange_phase_reference_adaptive_prod_v1`. Convergence artifacts remain under this result-side `convergence/` directory; the formal reference and figures use separate sibling directories.

## Predeclared acceptance gates

For C0 versus C1/C2 on the common initial xi/T anchors:

- first-order/crossover/CEP position differences: at most `0.10 MeV` in quark-chemical-potential or temperature coordinates;
- CEP baryon-chemical-potential difference: at most `0.30 MeV`;
- boundary/spinodal/crossover density difference: at most `0.01 rho0`;
- crossover response derivative relative difference: at most `5%`;
- no missing common anchor, non-numeric value, failed crossover row, or unconverged `phase_grid_convergence` record.

For C0 versus the C3 phase-line oracle at the five xi anchors and `T>=60 MeV`:

- first-order boundary difference: at most `0.30 MeV`;
- CEP and crossover position difference: at most `0.10 MeV` (`0.30 MeV` for `mu_B`);
- spinodal chemical-potential difference: at most `1.5 MeV` and must be reported as a diagnostic uncertainty;
- density difference: at most `0.02 rho0`.

O1 must pass its independent tensor-oracle self-convergence and adaptive-versus-oracle assertions. Any threshold failure yields `blocked` or `diagnostic-only`; it does not permit a formal rerun without an updated, author-reviewed convergence decision.
