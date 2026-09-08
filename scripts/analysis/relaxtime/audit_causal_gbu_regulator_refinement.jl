"""Refine the independent loop nodes without changing any closure tolerance."""
module CausalGBURegulatorRefinement
include("audit_causal_gbu_regulator_closure.jl")
const _CLOSURE_BASE=joinpath(CausalGBURegulatorClosure.R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
function main()
    ENV["GBU_CLOSURE_OUTPUT"]=get(ENV,"GBU_CLOSURE_OUTPUT",joinpath(_CLOSURE_BASE,"regulator_routing_closure_refined_20260905"))
    CausalGBURegulatorClosure.main(configurations=[(q,512) for q in (0.,1.,3.,6.)],np=512,nx=256)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
