# Endpoint-Local CEP Evidence

This group contains three endpoint-local diagnostic packages:

- `pnjl_cep_endpoint_local_contract_feasibility_v2/`: solver-free contract feasibility candidate.
- `pnjl_cep_endpoint_local_production_shadow_v4/`: full 24-anchor shadow from the earlier calculation line.
- `pnjl_cep_endpoint_local_production_shadow_v4_20260813/`: the later full-shadow evidence package.

The packages preserve separate calculation/postprocess provenance and verdicts. `feasible_candidate` and `full_hybrid_candidate` are diagnostic states, not reference promotion. Complete rho-mu curve data remains in the referenced Actions/local artifacts; this move changes only the repository namespace and current references.

The `phase_reference_current_state_freeze_v1` package retains generation-time paths and hashes as a frozen snapshot and is not rewritten by this namespace move.
