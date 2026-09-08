"""One shifted infinite-target bracket, using retained backgrounds only."""
module CausalGBUInfiniteMottRefinement
include("causal_gbu_infinite_thermal.jl")
const I=CausalGBUInfiniteThermal
const R=I.R
using CSV,JSON3
function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    previous=joinpath(base,"direction_b_infinite_mott_20260907")
    input=joinpath(base,"method_v1_mott")
    output=get(ENV,"GBU_INFINITY_MOTT_REFINED_OUTPUT",joinpath(base,"direction_b_infinite_mott_refined_20260907"))
    paths=[joinpath(previous,"normal_mott.csv"),joinpath(previous,"manifest.json"),joinpath(input,"backgrounds.csv")]
    input_hashes=Dict(p=>R.hashfile(p) for p in paths)
    rows=NamedTuple.(collect(CSV.File(first(paths))))
    hashes=R.start_output(output);probes=NamedTuple[]
    for T in (195.6640625,195.78125)
        bg=R.saved_background(input,T,24)
        k=I.kernel(bg,:K_plus,1.;cut_nodes=192,split_inv_fm=36.)
        lo=I.polarization(k,k.threshold;nodes=64);hi=I.polarization(k,k.threshold;nodes=128)
        push!(probes,(T_MeV=T,threshold_inverse=real(1-4bg.coupling[:K_plus]*hi),node_change=abs(hi-lo)))
    end
    l,h=probes
    index=only(findall(r->r.channel=="K_plus" && r.q_inv_fm==1.,rows))
    old=rows[index]
    passed=l.threshold_inverse<0<h.threshold_inverse && h.T_MeV-l.T_MeV<=0.25 &&
        max(l.node_change,h.node_change)<1e-6 && old.before_count==1 && old.after_count==0 &&
        abs(old.phase_drop_over_pi-1)<0.005
    rows[index]=merge(old,(T_low_MeV=l.T_MeV,T_high_MeV=h.T_MeV,
        low_threshold_inverse=l.threshold_inverse,high_threshold_inverse=h.threshold_inverse,
        node_change_inv_fm2=max(old.node_change_inv_fm2,l.node_change,h.node_change),passed=passed))
    CSV.write(joinpath(output,"normal_mott.csv"),rows);CSV.write(joinpath(output,"new_bracket.csv"),probes)
    all(R.hashfile(p)==h for (p,h) in input_hashes) || error("retained Mott evidence drift")
    R.finish_output(output,hashes,Dict("status"=>all(r.passed for r in rows) ? "normal_mott_infinite_accepted" : "mott_refinement_failed",
        "rows"=>length(rows),"all_passed"=>all(r.passed for r in rows),"input_hashes_absolute"=>input_hashes,
        "reused_successful_rows"=>7,"reused_count_phase_rows"=>8,"new_threshold_probes"=>2,
        "solver_called"=>false,"thermal_target"=>"infinity"))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
