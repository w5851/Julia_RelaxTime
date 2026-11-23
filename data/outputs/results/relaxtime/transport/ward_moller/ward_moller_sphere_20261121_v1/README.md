# Ward-matched A and Møller transport

## Source and adopted prescription

The complete production grid uses calculation commit
`2465c05d24d07c9c1c34949502e7243c4294bf15` and profile
`ward_moller_sphere_v1`. All 1867 original manuscript backgrounds are
preserved: 942 phase-referenced states and 925 fixed-temperature states.
The new calculation rebuilds all 31739 channel rates. Background masses,
Polyakov fields, densities, entropy and first-order branch choices are unchanged.

The terminal replay on W uses the frozen rates and the corrected
\(2n_{\bar u}\bar w_{s\bar u}\) term in the strange-quark rate. It does not
change any of the 31739 rates or any other species' relaxation times.
`rate_assembly_replay.json` binds the original run to the replayed data;
`rate_assembly_drift.json` records per-density maxima. The corresponding
direction and final-state blocking controls were also replayed. The earlier
independent integration refinements in `refinement_comparison.csv` and
`refinement_summary.json` retain the pre-correction terminal assembly.

The thermal one-line denominator is the variational RS energy, with the
vacuum term unchanged. The same A enters the coupling construction. The
incoming thermal flux uses medium-rest-frame Møller energies. Ordinary
on-shell dispersion and finite scattering support remain in the two-line
loop and external kinematics. The explicit linear two-line occupation is
averaged over the meson direction before constructing the propagator.

This restores the static soft-pion matching relation. It is not a derivation
of a finite-momentum Ward identity, a full nonlinear two-line loop, or an
average of the complete squared scattering amplitude.

## Retained evidence

- `numerical_manifest.json` binds the production SHA, configuration, Actions
  runs, returned remote identities and every retained payload.
- `execution_contracts.json` and `core_source_verification.json` retain the
  worker identities and the 256 identical remote/Actions core file hashes.
- `mode_a_grid.csv`, `mode_b_grid.csv` and `channel_rates.csv` are the complete
  native numerical results. `physics_audit.json` checks finiteness,
  symmetry, terminal replay, flavor ordering and conductivity decomposition.
- `spectra.csv` contains the zero-density pi/K poles on five frozen spectral
  backgrounds. `mott_temperatures.csv` uses direct threshold roots.
  The `threshold_inverse` diagnostic in the spectral scan can be infinite at
  the exact branch point of the native real-axis routine; it is not used to
  locate Mott roots or draw spectral curves.
- `frozen_mott_brackets.csv` confirms the pi/K intervals on the actual
  finite-density input rows. `dense_mott_roots.csv` is an earlier continuous
  diagnostic with separately solved backgrounds, not the source of Table III.
- `mechanism_*` stores eight complete local states, all process rates,
  energy-resolved elastic cross sections and weighted thermal areas.
- `direction_comparison.csv`, `direction_oracle.csv` and
  `refinement_summary.json` describe the new prescription's representative
  validation. `final_b0_poles.csv` and `final_blocking_comparison.csv` test
  retained physical approximations on the final Ward/sphere inputs.
- `midwindow_joint_*` records the four-point joint refinement around the
  138 MeV pion crossing. It verifies the tiny local maximum rather than
  changing the original production values.

## Results that changed

The smallest phase-referenced eta/s is 0.197, rather than the former 0.084.
All nine conductivity bands have a net decrease across the xi endpoints.
The zero-density intermediate bulk band has a small net decrease. At
160 MeV the intermediate/high-density shear ordering crosses at the
negative-xi edge, while both remain below the zero-density value.

The two pion intervals are (-0.04,-0.03) and (0.17,0.18). Elastic
d-ubar scattering supplies 59.1% and 58.6% of the inverse ubar rate
decrease. The ubar relaxation time grows by 21.5% and 13.0%; the total
conductivity changes by -0.76% and +0.05%. Joint refinement retains a
small maximum in the second window, with a +0.08% crossing change.
The u+d sector carries about 98%/97% of the current, and the antiquark
sector is strongly filtered by both occupation and relaxation weights.

In the phase-referenced bands at 450 MeV, the quark times cross near the
negative-xi edge. At 900 MeV, tau_s stays below tau_u throughout all three
bands. In the coolest and intermediate high-density bands, the strange-quark
contribution from light targets, `2 n_u w_us`, is already 1.35–1.60 and
1.26–1.33 times the complete light-quark collision rate, respectively.
The extra columns in `flavor_hierarchy.csv` close the relaxation-ordering
comparison using the same channel rates and target densities.

At xi=0.18 the pion pole is only about 0.00038 MeV below threshold in the
native real-axis diagnostic. The narrow cross-section peak nevertheless
has a smaller thermally weighted contribution: the 0–25 MeV d-ubar
collision-rate contribution falls from 0.09695 to 0.05052 fm^-1.
`physics_audit.json` records this local pole and energy-window comparison.

## Figures and review state

The five new manuscript figures are in
`data/outputs/figures/relaxtime/transport/phase_guided/ward_moller_manuscript_20261121_v1/`.
They use a single bundle manifest, native vertices, the original branch
gaps, 600 dpi color/gray PNGs and Times/STIX fonts. Figures 2–4 are
replaced; the new mechanism figure is Figure 5 and fixed T is Figure 6.
The previously accepted C2 phase figure is reused by the paper.

The author-requested figure corrections replace the PNGs in this same
bundle: collision-rate terminology, italic flavor labels, common rate
limits, expanded cross-section headroom, clear conductivity-change ticks,
inset shared legends, and two-decimal bulk ticks at fixed T=200 MeV.
The terminal replay updates Figs. 3–6 while preserving the layouts and
leaving the zero-density spectrum in Fig. 2 unchanged.

Numerical validation, machine plot validation, manuscript compilation and
author acceptance are separate. The new evidence and PNGs are in
`agent_review`/`png_review`; the old v15 author acceptance is historical.
