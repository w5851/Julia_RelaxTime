## Evidence boundary

- v5 tables are copied byte-for-byte into `tables/v5_*` and remain the baseline.
- Endpoint expansion summary is copied from the solver-free replay artifact; the full response curves remain in the external artifact referenced by the replay manifest.
- The plotting script does not call the PNJL solver and does not write a reference.
- The figure is a visual diagnostic. It does not convert unresolved v5 geometry into a certificate and does not infer missing xi slices.
