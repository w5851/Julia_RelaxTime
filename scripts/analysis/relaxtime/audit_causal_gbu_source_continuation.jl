"""Retained-BQS source checks across individual quark onset; no new solves."""
module CausalGBUSourceContinuationAudit
include("causal_gbu_source_continuation.jl")
const C=CausalGBUSourceContinuation
const R=C.R
using CSV

function main(;nodes=128,channels=R.CHANNELS,qs=(0.,1.))
    nodes>=32 && all(c->c in R.CHANNELS,channels) && all(q->isfinite(q) && q>=0,qs) ||
        throw(ArgumentError("invalid source continuation coverage"))
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=joinpath(base,"method_v1_mott")
    bgpath=joinpath(input,"backgrounds.csv");bh=R.hashfile(bgpath)
    output=get(ENV,"GBU_SOURCE_CONTINUATION_OUTPUT",joinpath(base,"direction_b_source_continuation_20260907"))
    R.method_contract();hashes=R.start_output(output)
    rows=NamedTuple[]
    for ch in channels,q in qs
        bg=R.saved_background(input,205.625,24)
        i,j=R.charged_rpa_spec(ch).pair
        a,u,b,v=bg.m[i],bg.mu[i],bg.m[j],bg.mu[j]
        rs=[C.source_response(q,a,u,b,v,bg.T,bg.vacuum,10.,[0.,0.8im,2+0.8im];
            Phi=bg.Phi,PhiBar=bg.PhiBar,nodes=n,source_steps=[0.004,0.002,0.001]) for n in (div(nodes,2),nodes)]
        g=R.build_spectral_bubble(q,a,u,b,v,bg.T;Phi=bg.Phi,PhiBar=bg.PhiBar,
            vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=10.,momentum_nodes=2nodes,angle_nodes=nodes)
        dynamic=[R.spectral_bubble(g,real(z);eta_inv_fm=imag(z)).value for z in (0.8im,2+0.8im)]
        node=maximum(abs.(rs[2].pi_inv_fm2-rs[1].pi_inv_fm2))
        parity=maximum(abs.(dynamic-rs[2].pi_inv_fm2[2:3]))
        curvature=abs(last(rs[2].source_curvatures_inv_fm2)+real(first(rs[2].pi_inv_fm2)))
        step=abs(diff(rs[2].source_curvatures_inv_fm2)[end])
        curvnode=maximum(abs.(rs[2].source_curvatures_inv_fm2-rs[1].source_curvatures_inv_fm2))
        passed=max(node,parity,curvature,step,curvnode)<1e-6
        push!(rows,(channel=String(ch),q_inv_fm=q,T_MeV=205.625,muB_MeV=240.,
            reference_gap_lower_bound_inv_fm=rs[2].reference_gap_lower_bound_inv_fm,
            old_pair_proof_applicable=rs[2].below_individual_onset,
            source_node_change_inv_fm2=node,dynamic_parity_error_inv_fm2=parity,
            curvature_error_inv_fm2=curvature,source_step_change_inv_fm2=step,
            curvature_node_change_inv_fm2=curvnode,numerics_passed=passed,production_authorized=false))
        CSV.write(joinpath(output,"source_onset.csv"),rows)
        println("[source-onset] $(ch) q=$(q) passed=$(passed) max=$(max(node,parity,curvature,step,curvnode))");flush(stdout)
    end
    R.hashfile(bgpath)==bh || error("input drift")
    R.finish_output(output,hashes,Dict("status"=>all(r.numerics_passed for r in rows) ? "source_onset_probes_passed" : "source_onset_probe_failure",
        "input_directory"=>input,"input_hashes"=>Dict("backgrounds.csv"=>bh),"rows"=>length(rows),
        "numerics_passed"=>all(r.numerics_passed for r in rows),"nodes"=>[div(nodes,2),nodes],
        "channels"=>String.(channels),"q_inv_fm"=>collect(qs),
        "source_steps_inv_fm"=>[0.004,0.002,0.001],"threshold_inv_fm2"=>1e-6,
        "solver_called"=>false,"meson_density_computed"=>false,"full_stationarity_certified"=>false,
        "reference"=>"fixed_rank_four_Dirac_sea_not_zero_temperature_Fermi_medium"))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
