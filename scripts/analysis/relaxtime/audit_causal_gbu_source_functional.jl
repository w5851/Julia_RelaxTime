"""B1 source-Hessian audit on a retained background; no solve or density."""
module CausalGBUSourceFunctionalAudit
include("causal_gbu_source_functional.jl")
const F=CausalGBUSourceFunctional
const R=F.R
using CSV,JSON3

function main(;nodes=128)
    nodes>=16 && iseven(nodes) || throw(ArgumentError("even nodes>=16 required"))
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_SOURCE_OUTPUT",joinpath(base,"direction_b_source_v2_20260906"))
    hashes=R.start_output(output)
    rows=NamedTuple[]; statics=NamedTuple[]
    zs=ComplexF64[0.8im,2+0.8im,6+1.1im]
    tolerance=1e-6
    for c in (:pi_plus,:K_plus),q in (0.,1.4,3.2),occ in ("fermi","pnjl"),
        tc in (bg.vacuum,10.),ch in (:P,:S)
        i,j=R.charged_rpa_spec(c).pair
        a,u,b,v=bg.m[i],bg.mu[i],bg.m[j],bg.mu[j]
        phi,bar=occ=="fermi" ? (1.,1.) : (bg.Phi,bg.PhiBar)
        args=(q,a,u,b,v,bg.T,bg.vacuum,tc)
        coarse=F.source_response(args...,zs;Phi=phi,PhiBar=bar,channel=ch,nodes=div(nodes,2))
        fine=F.source_response(args...,zs;Phi=phi,PhiBar=bar,channel=ch,nodes=nodes)
        grid(n,swapped=false)=R.build_spectral_bubble(q,swapped ? b : a,swapped ? v : u,
            swapped ? a : b,swapped ? u : v,bg.T;Phi=phi,PhiBar=bar,channel=ch,
            vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=tc,
            momentum_nodes=n,angle_nodes=n)
        g,gsmall,greverse=grid(nodes),grid(div(nodes,2)),grid(nodes,true)
        ep,en,ef,ec=0.,0.,0.,0.
        for (k,z) in enumerate(zs)
            cur=R.spectral_bubble(g,real(z);eta_inv_fm=imag(z))
            old=R.spectral_bubble(gsmall,real(z);eta_inv_fm=imag(z))
            rev=R.spectral_bubble(greverse,-real(z);eta_inv_fm=imag(z))
            ep=max(ep,abs(cur.value-fine.pi_inv_fm2[k]))
            en=max(en,abs(cur.value-old.value),abs(fine.pi_inv_fm2[k]-coarse.pi_inv_fm2[k]))
            ef=max(ef,abs(cur.value-conj(rev.value)))
            ec=max(ec,cur.contact_identity_residual)
        end
        push!(rows,(channel=String(c),spin_channel=String(ch),occupation=occ,q_inv_fm=q,
            thermal_cutoff_inv_fm=tc,nodes=nodes,source_core_error_inv_fm2=ep,
            node_change_inv_fm2=en,flavor_reflection_error_inv_fm2=ef,
            contact_identity_error_inv_fm2=ec,passed=max(ep,en,ef,ec)<tolerance,
            solver_called=false,production_authorized=false))
        # Independent eigenvalue trace-log curvature, not a phase derivative.
        if ch===:P && q<=1.4
            steps=[0.004,0.002,0.001]
            static=F.source_response(args...,[0.];Phi=phi,PhiBar=bar,nodes=64,source_steps=steps)
            small=F.source_response(args...,[0.];Phi=phi,PhiBar=bar,nodes=32)
            cs=static.source_curvatures_inv_fm2
            error=abs(cs[end]+real(static.pi_inv_fm2[1]))
            stepchange=abs(cs[end]-cs[end-1])
            nodechange=abs(static.pi_inv_fm2[1]-small.pi_inv_fm2[1])
            push!(statics,(channel=String(c),occupation=occ,q_inv_fm=q,
                thermal_cutoff_inv_fm=tc,pi_static_inv_fm2=real(static.pi_inv_fm2[1]),
                curvature_coarse_inv_fm2=cs[1],curvature_middle_inv_fm2=cs[2],
                curvature_fine_inv_fm2=cs[3],hessian_error_inv_fm2=error,
                source_step_change_inv_fm2=stepchange,node_change_inv_fm2=nodechange,
                passed=max(error,stepchange,nodechange)<tolerance,
                full_stationarity_certified=false,production_authorized=false))
        end
        CSV.write(joinpath(output,"dynamic_hessian.csv"),rows)
        isempty(statics) || CSV.write(joinpath(output,"static_source_curvature.csv"),statics)
        println("[source-functional] $(length(rows))/48 $(c) $(ch) q=$(q) $(occ) Lth=$(tc) passed=$(rows[end].passed)")
    end
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    R.finish_output(output,hashes,Dict("status"=>"B1_fixed_projector_functional_diagnostic",
        "background"=>bg,"dynamic_rows"=>length(rows),"static_rows"=>length(statics),
        "dynamic_passed"=>all(r.passed for r in rows),"static_passed"=>all(r.passed for r in statics),
        "absolute_tolerance_inv_fm2"=>tolerance,"external_frequency_probes"=>string.(zs),
        "static_source_steps_inv_fm"=>[0.004,0.002,0.001],"static_nodes"=>[32,64],
        "functional"=>"fixed_projector_vacuum_plus_medium_minus_same_source_reference",
        "physical_response_selected"=>false,"positive_physical_trace_certified"=>false,
        "full_charged_action_stationarity_certified"=>false,"continuous_UHP_count_certified"=>false,
        "GBU_observable_derived"=>false,"meson_density_computed"=>false,"solver_called"=>false,
        "limitations"=>["Fixed-background fermion Hessian only; no PNJL or KMT full-action Hessian",
            "Signed subtraction determinant is not a positive physical ensemble",
            "PNJL holonomy is complexified for unequal Phi and PhiBar",
            "No stability, Mott, Levinson, or density acceptance follows from this comparison"]))
    println("[source-functional] dynamic $(count(r->r.passed,rows))/$(length(rows)); static $(count(r->r.passed,statics))/$(length(statics))")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
