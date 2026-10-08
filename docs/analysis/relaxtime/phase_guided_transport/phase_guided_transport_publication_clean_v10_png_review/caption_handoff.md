# publication_clean_v10 PNG review caption handoff

v10 changes figure layout and display transforms only. It does not add points,
replace values, reconnect phase branches, smooth data, call a solver, or
change the v5 provenance chain.

Recommended caption additions:

> Columns correspond to $\mu_B=0$, 450, and 900 MeV from left to right.
> Solid, dashed, and dash-dotted curves correspond to $\alpha_T=1.0$,
> 1.1, and 1.2, respectively, in both figures. The associated fixed
> temperatures are listed below. Figure 1 rows show $\tau_u$, $\tau_s$,
> $\tau_{\bar u}$, and $\tau_{\bar s}$ in fm; Figure 2 rows show
> $\eta/s$, $\zeta/s$, and $\sigma/T$. Open circles and squares mark the chirally restored and
> chirally broken endpoints of the first-order transition, respectively.
> At $\mu_B=900\;\mathrm{MeV}$ and $\alpha_T=1.0$, disconnected segments
> display the phase branches separately. In Figure 1 the complete
> $\mu_B=900\;\mathrm{MeV}$ column uses logarithmic y axes; other
> composite panels use linear y axes. Vertical ranges differ between panels.

The in-axes legend is a reviewed layout exception for this dense multi-panel
family. It is not a new repository-wide default: the default remains an
external legend unless a case-level contract records and reviews the exception.
The parameter key is in (a); the endpoint key is in (c), under the heading
"1st-order". The keys appear once each, with no key placed above the figure.

Typography status: the composites use 11 pt labels/titles/ticks and a 10.5 pt
legend. Ordinary capital/numeral glyphs exceed 2 mm at the native 6.75 in
width, but math scripts are about 1.69 mm. This is an explicitly registered
PNG-review exception, not a claim of full APS typography compliance. No PDF
or submission eligibility is granted by this stage.

Exact frozen temperature mapping (MeV):

| panel | alpha_T=1.0 | alpha_T=1.1 | alpha_T=1.2 |
| --- | ---: | ---: | ---: |
| muB0.0 | 200.088901 | 220.097792 | 240.106682 |
| muB450.0 | 182.763066 | 201.039372 | 219.315679 |
| muB900.0 | 125.737258 | 138.310984 | 150.884710 |

Panels using local log-y in v10:

- `mode_a` / `muB900.0` / `tau_u`: ratio=10.9391, positive data, local log-y.
- `mode_a` / `muB900.0` / `tau_s`: ratio=9.59061, positive data, local log-y.
- `mode_a` / `muB900.0` / `tau_ubar`: ratio=48.4833, positive data, local log-y.
- `mode_a` / `muB900.0` / `tau_dbar`: ratio=48.4833, positive data, local log-y.
- `mode_a` / `muB900.0` / `tau_sbar`: ratio=33.8425, positive data, local log-y.
- `mode_b` / `T120.0` / `eta_over_s`: ratio=146.168, positive data, local log-y.
- `mode_b` / `T120.0` / `zeta_over_s`: ratio=235.144, positive data, local log-y.
- `mode_b` / `T120.0` / `tau_u`: ratio=514.389, positive data, local log-y.
- `mode_b` / `T120.0` / `tau_d`: ratio=514.389, positive data, local log-y.
- `mode_b` / `T120.0` / `tau_s`: ratio=1003.4, positive data, local log-y.
- `mode_b` / `T120.0` / `tau_ubar`: ratio=11098, positive data, local log-y.
- `mode_b` / `T120.0` / `tau_dbar`: ratio=11098, positive data, local log-y.
- `mode_b` / `T120.0` / `tau_sbar`: ratio=2641.52, positive data, local log-y.
- `mode_b` / `T160.0` / `zeta_over_s`: ratio=56.4808, positive data, local log-y.
- `mode_b` / `T160.0` / `tau_u`: ratio=34.0882, positive data, local log-y.
- `mode_b` / `T160.0` / `tau_d`: ratio=34.0882, positive data, local log-y.
- `mode_b` / `T160.0` / `tau_s`: ratio=44.8428, positive data, local log-y.
- `mode_b` / `T160.0` / `tau_ubar`: ratio=337.093, positive data, local log-y.
- `mode_b` / `T160.0` / `tau_dbar`: ratio=337.093, positive data, local log-y.
- `mode_b` / `T160.0` / `tau_sbar`: ratio=222.436, positive data, local log-y.

The v10 PNG review layer remains manuscript_eligible=false and does not replace the v5/current pointer.
