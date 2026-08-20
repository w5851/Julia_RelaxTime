# CEP And Maxwell Evidence

This directory is the physical namespace for the Issue #130 CEP/Maxwell evidence line. The subgroups preserve the distinction between pilot, feasibility, shadow, endpoint-local, and Maxwell-contract evidence:

- `narrow_pilot/`: v1 and v2 local CEP pilot packages.
- `stagec/`: Stage-C offline, certificate, extrema-guard, and tolerance replay packages.
- `cascade_shadow/`: production-cascade shadow and the earlier full endpoint-local shadow.
- `endpoint_local/`: endpoint-local contract feasibility and later full-shadow packages.
- `maxwell_contracts/`: endpoint candidate, endpoint-limit, and tolerance-contract packages.

Each versioned package keeps its own manifest, source provenance, hashes, tables, figures, and verdict. The directory grouping does not merge evidence or promote diagnostic candidates. `phase_reference_*` packages remain at the PNJL root because they are decision/freeze layers that summarize these inputs rather than belong to one CEP/Maxwell run.
