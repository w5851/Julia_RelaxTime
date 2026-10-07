# Prospective sampling and coding protocol

Recorded on 2026-10-05 before inspecting any sampled paper figures.

Question: In a bounded sample of close PNJL/NJL transport studies and neighboring momentum-anisotropic QCD transport studies, how often are legends outside the top of a whole figure, and which other design choices can inform the current 12-panel/9-panel transport figures?

Primary source: INSPIRE public literature API. arXiv/INSPIRE open full text for visual inspection; Crossref DOI metadata as a secondary bibliographic check if available. This is a targeted descriptive visual audit, not a systematic review of all QCD papers.

Search window: publication years 2015-2026, with records available on 2026-10-05. Original research in Physical Review D or C. The three query strings and raw responses are retained in this directory. Merge by INSPIRE ID and check DOI/arXiv duplicates.

Strata, fixed before figures are seen:

- Core: 8 most recently published original PNJL/NJL studies that calculate transport coefficients or relaxation times in quark matter. Model use and transport calculation must be confirmed in the abstract/full text, not merely cited. Exclude reviews, proceedings, theses, hydrodynamic/flow studies without coefficients, and thermodynamics-only papers.
- Adjacent: 8 most recently published original studies of momentum-anisotropic quark/QCD transport coefficients or relaxation times, using kinetic/RTA or a related effective model. Exclude magnetic-field-driven-only anisotropy, hydrodynamic flow only, jet transport only, and astrophysical rotation/tidal transport. A paper eligible for both strata enters core once.
- Rank within strata by journal publication year descending, then preprint date descending, then INSPIRE ID. If fewer than 8 are eligible, audit all; do not expand by visual appearance. If an open full text is unavailable, record the failed retrieval and use the next eligible paper under the same rule.

Count every numbered figure containing a two-dimensional quantitative line/point plot in selected papers. Record noneligible numbered figures and exclusions explicitly. A numbered figure is the unit; multi-panel figures are not split into independent samples. Record panel count and transport-specific status for sensitivity analyses. Record paper-level counts as well as figure-level counts; no independence claim or population confidence interval.

Mutually exclusive primary legend-location codes:

- inside: every series key is inside plotting axes (including a key appearing in only one panel and serving other panels).
- top_external: at least one key is above the whole plotting area; record whether other inside keys coexist and whether it is global/shared or axis-local.
- other_external: an external key beside, below, or between plotting axes, without a top_external key.
- direct_or_caption: no series key, but labels on curves or caption/text provide the mapping.
- no_key_needed: no grouped series requiring a key.
- uncertain: insufficient resolution or ambiguous key/border; resolve with a larger render or leave uncertain.

Do not count a condition label (e.g. fixed temperature), panel label, title, or axis offset as a legend. A key inside the top corner is inside, not top_external. Record global two-row header legends comparable to v12 separately from any top_external key.

For each eligible figure record DOI/arXiv, PDF version/source/hash, printed figure number, PDF page, panel count, key location and sharing, visible border/tick treatment, grid, color/dash/marker encoding, parameter annotation, axis scale, detail view, and observations. Do not infer source font sizes, automatic placement, or exact final journal millimeters from PDF appearance alone. Do not infer grayscale usability from color alone.

Keep selection decisions, PDF figure coverage, coding, and arithmetic reproducible. Report sampling, author-group, publication-version, and single-reviewer limitations. Separate observed practice, official journal requirements, and proposed local SOP changes. Preserve the frozen v12 files and SOP inputs during this research-only audit.

## Metadata-screening clarification (before PDF inspection)

The broad anisotropic query returned 880 records and only the newest 250 were retrieved as a pilot. It is not the sampling frame. The final anisotropic query adds `date 2015->2026 and (j Phys.Rev.D or j Phys.Rev.C)` and returns all 148 hits. The final union with the complete PNJL (39) and NJL (165) queries contains 317 unique INSPIRE records before journal/year/topic screening.

To keep the core stratum close to the manuscript, exclude studies dominated by imposed magnetic fields, rigid rotation, cold color superconductivity/neutron-star matter, or chiral-charge equilibration only. Hot quark matter with chiral imbalance remains relevant when shear/bulk transport are computed. In the adjacent stratum, retain field-response studies only when they explicitly calculate an additional momentum-distribution anisotropy, and retain heavy-quark drag/diffusion as neighboring (not identical) observables. Exclude pure attractor/evolution studies whose objective is recovering known coefficients rather than studying a transport-coefficient dependence. These are topical decisions made from metadata, not figure style.

Selected core INSPIRE IDs: 3105682, 2842485, 2767151, 2902771, 1828803, 1717205, 1793578, 1694808.
Selected adjacent IDs: 3079053, 2973399, 2668341, 1968667, 1905653, 1846668, 1836409, 1778152.
