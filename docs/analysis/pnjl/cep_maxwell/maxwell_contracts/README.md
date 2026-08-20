# Maxwell Contract Evidence

This group contains three solver-free Maxwell endpoint/tolerance contract packages:

- `pnjl_maxwell_endpoint_candidate_feasibility_v1/`: strict three-intersection candidate replay.
- `pnjl_maxwell_endpoint_limit_contract_v1/`: endpoint-limit refinement certificate.
- `pnjl_maxwell_tolerance_contract_feasibility_v1/`: tolerance-frontier feasibility replay.

These packages have separate runs, inputs, and verdict boundaries. `feasible_candidate` and endpoint-limited candidate states are diagnostic evidence, not production or phase-reference promotion. The move changes the repository namespace and current generator paths only; package-internal generation commands and external artifact paths remain generation-time provenance.
