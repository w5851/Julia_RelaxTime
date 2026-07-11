# Models API Taxonomy

Read this reference only when reorganizing `docs/api/models/` topics or migrating historical API pages.

## Migration policy

- Absorb still-valid parameter conventions, contracts, and cautions from historical `docs/api/pnjl/` or `docs/api/relaxtime/` pages into the new topic.
- Reduce absorbed historical pages to migration notes or redirects; do not make new topics depend on them as their primary explanation.
- Treat references to the removed `src/pnjl/PNJL.jl` path as compatibility history, not a current implementation entrypoint.

## Topic grouping

- Organize model families and variants under an explicit family/variant category rather than flattening every topic directly under `models/`.
- Treat magnetic-field capabilities as model variants alongside NJL, PNJL, and RPNJL families.
- Group susceptibility, cumulants, and derivative capabilities as derived quantities rather than workflow peers of phase, scans, and solvers.
- Treat transport as a derived capability conceptually, but do not force a final `docs/api/models/...` location while it remains coupled to relaxtime domain and workflow pages.

## Responsibility-core filenames

- Use `Algorithms.md` when the topic is dominated by numerical or physical criteria.
- Use `CoreConcepts.md` when the topic is dominated by responsibilities, data flow, and module collaboration.

The semantic layer is stable; filenames may follow the topic's dominant reader need.
