# Thermal relative-velocity convention audit (2026-10-24)

## Scope

This is a pre-production convention audit. It does not replace the current
production profile, regenerate production CSVs, or dispatch GitHub Actions.
The implementation worktree is `codex/meson-freezeout-defaults` at base
`dffe4070eb5c0788e9c1a227cf40978b51ce3e`; the relative-velocity changes are
currently uncommitted.

## Formula and implementation

For an incoming thermal event with medium-frame energies
`e_i = sqrt(p_i^2 + m_i^2)` and `e_j = sqrt(p_j^2 + m_j^2)`, the two supported
conventions are

\[
 v_{\rm cm}=\frac{\sqrt{\lambda(s,m_i^2,m_j^2)}}{2E_i^*E_j^*},
 \qquad
 v_{\rm M}=\frac{\sqrt{\lambda(s,m_i^2,m_j^2)}}{2e_i e_j}.
\]

`cm` is the historical production convention. `medium_moller` evaluates the
second denominator before the event is compressed to its `s` node, so it cannot
be reconstructed by multiplying a final rate by one global factor.

The mode is part of `ThermalSamplingConfig` and therefore part of the thermal
cache identity. A cache built with one mode is rejected when requested with the
other mode. The default remains `cm` until the paper/production convention is
explicitly selected.

## Fixed algebra check

For `m_i=m_j=1`, `|p_i|=|p_j|=1`, and perpendicular incoming momenta,
`e_i=e_j=sqrt(2)` and `s=6`. The implementation returns

* `v_cm = 1.1547005383792515 = 2/sqrt(3)`;
* `v_M = 0.8660254037844384 = sqrt(3)/2`;
* `v_cm/v_M = 4/3`.

This check passed locally and on W. It is an algebraic convention difference,
not a quadrature or threshold-sampling error.

## W representative-point probe

The probe ran on host W through the `quote-remote` SSH alias, using an isolated
clone at `C:\\Users\\82306\\Desktop\\Julia\\Julia\\transport_probe_20261024`.
It used a small fixed-background fixture (`T=0.9 fm^-1`, `xi=0.27`,
`m=(0.73,0.73,1.51) fm^-1`, `mu=(0.2,0.2,0.5) fm^-1`), `p_nodes=4`,
`cos_nodes=3`, `phi_nodes=5`, `n_sigma_points=4`, and `sigma_cutoff=2.0`.
These are diagnostic settings, not production settings.

For the `us_to_us` channel:

| mode | averaged rate | estimated cache error | evaluations | cache nodes |
|---|---:|---:|---:|---:|
| `cm` | 0.08887235505499728 | 4.575774954069018e-5 | 115 | 83 |
| `medium_moller` | 0.05414652548183448 | 3.2314607354700505e-5 | 115 | 82 |

The Møller rate is `0.6092617378` of the historical rate in this diagnostic
point. A full small-grid relaxation-time probe with the existing unit-test
density fixture gave `tau_u=0.3959448873` (`cm`) versus `0.6960130529`
(`medium_moller`), with analogous increases for the other species. These
values establish sensitivity to the convention; they are not production or
paper benchmark values.

The first `cm` probe included W-side Julia dependency compilation. The two
timings (`7.09 s` for the first rate cache, `0.021 s` for the warm Møller cache;
`10.79 s` and `2.47 s` for the two small-grid full tau calls) are therefore not
a fair performance comparison.

## Local verification

The following focused tests passed after the change:

* `test_thermal_optimization_native.jl`: 51/51 assertions;
* `test_thermal_weighted_cache.jl`: 40/40 assertions;
* `test_thermal_weighted_sampling.jl`: 59/59 assertions;
* `test_average_scattering_rate.jl`: 31/31 assertions;
* `test_transport_sampling_profiles.jl`: 25/25 assertions.

The production sampling profile still resolves to `relative_velocity = :cm`
through the constructor default. No production output or Actions run has been
modified by this audit.

## Relation to the production-plan review

The review in Codex thread `01a16bbe-2c21-74a3-a7cb-2c43c0d52e05`
recommended the bounded transport-only route: keep the existing mean-field and
phase background, use the complete RS one-loop input consistently for the
transport-side effective coupling, and record the remaining background closure
as an approximation. That recommendation does not settle the incoming
relative-velocity convention. The manuscript interface audit is explicit that
the external legs use ordinary on-shell energies and defines

\[
v_{\rm rel}=\frac{\sqrt{(P_i\!\cdot P_j)^2-M_i^2M_j^2}}
                 {\varepsilon_i\varepsilon_j},
\]

which is `medium_moller` for the thermal event. Therefore the complete-RS-
coupling decision and the `cm`/Møller decision are separate production gates;
passing the former does not make the current historical `cm` results formula
aligned with the manuscript.

## Decision gate

Before a production rerun, select one of these explicit routes:

1. retain `cm` and state the thermal-average approximation in the manuscript;
2. promote `medium_moller`, add it explicitly to the production profile, and
   repeat the independent/joint convergence, Julia–Fortran alignment, and
   production pilot checks under the new cache identity.

Until that decision is made, existing production results remain tagged as the
historical `cm` convention and should not be compared to a Møller rerun as if
the difference were numerical error.
