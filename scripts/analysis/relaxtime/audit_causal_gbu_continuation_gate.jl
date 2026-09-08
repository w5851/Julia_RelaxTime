"""Conditional normal Mott checks; full-spectrum/production gates remain separate."""
module CausalGBUContinuationGateAudit
include("causal_gbu_continuation_gate.jl")
const C=CausalGBUContinuationGate
const A=C.A
const R=C.R
using CSV,JSON3

function kernel(bg,ch,q,thermal)
    i,j=R.charged_rpa_spec(ch).pair
    g=R.build_spectral_bubble(q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=thermal,momentum_nodes=16,angle_nodes=8)
    edges=A.cut_panels(g)
    rho=x->x isa BigFloat ? C.wide_cut(bg,ch,q,x,thermal;nodes=64) : A.direct_cut(bg,ch,q,x,thermal;nodes=64).imaginary
    return (;g,edges,rho,shift=bg.mu[i]-bg.mu[j],
        threshold=hypot(q,bg.m[i]+bg.m[j]))
end

function direct_value(k,z)
    values=ComplexF64[]
    fallback=false
    for n in (64,128)
        v=try
            C.mapped_cauchy(k.rho,k.edges,z;nodes=n)
        catch err
            err isa ArgumentError && occursin("unrepresentable",sprint(showerror,err)) || rethrow()
            fallback=true
            p=A.continuous_cauchy(k.rho,k.edges,z;nodes=n==32 ? 8 : 12,atol=n==32 ? 1e-7 : 1e-8)
            p.converged || error("adaptive fallback not converged")
            p.value
        end
        push!(values,v)
    end
    node_change=abs(values[1]-values[2])
    return (value=last(values),node_change=node_change,passed=node_change<1e-6,
        endpoint_fallback_used=fallback)
end

function normal_gap(bg,ch,q,thermal;loop_nodes=(128,256))
    length(loop_nodes)==2 && 8<=first(loop_nodes)<last(loop_nodes) || throw(ArgumentError("two increasing loop orders required"))
    i,j=R.charged_rpa_spec(ch).pair
    a,b,u,v=bg.m[i],bg.m[j],bg.mu[i],bg.mu[j]
    lo=max(0.,hypot(q,a-b)-u+v);hi=hypot(q,a+b)-u+v
    checks=map(loop_nodes) do n
        g=R.build_spectral_bubble(q,a,u,b,v,bg.T;Phi=bg.Phi,PhiBar=bg.PhiBar,
            vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=thermal,
            momentum_nodes=n,angle_nodes=div(n,2))
        f(w)=1-4bg.coupling[ch]*R.spectral_bubble(g,w;eta_inv_fm=0.).value
        r=R.certify_gap_roots((w,_)->f(w),q,[(lo,hi)];physical_sheet=true,real_axis=true,omega_nodes=128)
        slopes=[real(-4bg.coupling[ch]*R.spectral_bubble(g,x.omega_inv_fm;eta_inv_fm=0.).derivative) for x in r.roots]
        return (r=r,slopes=slopes)
    end
    same=checks[1].r.count==checks[2].r.count
    change=same ? maximum((abs(a.omega_inv_fm-b.omega_inv_fm) for (a,b) in zip(checks[1].r.roots,checks[2].r.roots));init=0.) : Inf
    k=kernel(bg,ch,q,thermal)
    phases=Float64[];num=true
    for h in (1e-9,1e-11)
        p=direct_value(k,complex(k.threshold+h))
        push!(phases,-angle(1-4bg.coupling[ch]*p.value)/pi)
        num &= p.passed
    end
    return (count=checks[2].r.count,root_change=change,phase=last(phases),
        direct_phase_numerics_passed=num,gap_root_checks_passed=all(x.r.passed for x in checks),
        phase_limit_change=abs(phases[1]-phases[2]),
        passed=all(x.r.passed for x in checks) && same && change<1e-7 &&
            all(<(0),checks[2].slopes) && num && abs(phases[1]-phases[2])<0.0005)
end

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=joinpath(base,"method_v1_mott")
    bgpath=joinpath(input,"backgrounds.csv")
    bh=R.hashfile(bgpath)
    saved=collect(CSV.File(bgpath))
    ts=sort!(unique([r.T_MeV for r in saved if r.p_nodes==24 && r.muB_MeV==240.]))
    output=get(ENV,"GBU_CONTINUATION_OUTPUT",joinpath(base,"direction_b_continuation_20260907"))
    R.method_contract()
    hashes=R.start_output(output)
    probes,mott,finiteq,failures=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    cache=Dict{Tuple{Symbol,Float64,Float64},Any}()
    function probe(ch,q,T)
        key=(ch,q,T)
        get!(cache,key) do
            bg=R.saved_background(input,T,24)
            k=kernel(bg,ch,q,10.)
            p=direct_value(k,complex(k.threshold))
            dom=C.source_domain(bg,ch)
            row=merge((channel=String(ch),q_inv_fm=q,T_MeV=T,
                threshold_inverse=real(1-4bg.coupling[ch]*p.value),
                threshold_pi_node_change_inv_fm2=p.node_change,numerics_passed=p.passed,
                endpoint_fallback_used=p.endpoint_fallback_used,production_authorized=false),dom)
            push!(probes,row)
            return (bg=bg,row=row)
        end
    end
    function checkpoint()
        for (name,rows) in (("threshold_probes",probes),("normal_mott",mott),
                           ("finite_q_endpoints",finiteq),("failures",failures))
            isempty(rows) || CSV.write(joinpath(output,name*".csv"),rows)
        end
        flush(stdout)
    end
    for ch in R.CHANNELS,q in (0.,1.)
        try
            left,right=1,length(ts)
            low,high=probe(ch,q,ts[left]),probe(ch,q,ts[right])
            low.row.threshold_inverse<0<high.row.threshold_inverse || error("no retained threshold bracket")
            while right-left>1
                mid=div(left+right,2)
                p=probe(ch,q,ts[mid])
                if p.row.threshold_inverse<0
                    left,low=mid,p
                else
                    right,high=mid,p
                end
            end
            preT=maximum(filter(T->T<=ts[left]-0.5,ts))
            postT=minimum(filter(T->T>=ts[right]+0.5,ts))
            prebg,postbg=R.saved_background(input,preT,24),R.saved_background(input,postT,24)
            pre,post=normal_gap(prebg,ch,q,10.),normal_gap(postbg,ch,q,10.)
            domain=C.source_domain(low.bg,ch).previous_pair_source_proof_applicable &&
                C.source_domain(high.bg,ch).previous_pair_source_proof_applicable &&
                C.source_domain(prebg,ch).previous_pair_source_proof_applicable &&
                C.source_domain(postbg,ch).previous_pair_source_proof_applicable
            passed=low.row.numerics_passed && high.row.numerics_passed &&
                ts[right]-ts[left]<=0.25 && pre.passed && post.passed &&
                pre.count-post.count==1 && abs(pre.phase-post.phase-1)<0.005
            push!(mott,(channel=String(ch),q_inv_fm=q,thermal_cutoff_inv_fm=10.,
                T_low_MeV=ts[left],T_high_MeV=ts[right],T_pre_MeV=preT,T_post_MeV=postT,
                before_count=pre.count,after_count=post.count,phase_drop=pre.phase-post.phase,
                before_root_node_change=pre.root_change,after_root_node_change=post.root_change,
                before_phase_node_passed=pre.direct_phase_numerics_passed,after_phase_node_passed=post.direct_phase_numerics_passed,
                before_gap_checks_passed=pre.gap_root_checks_passed,after_gap_checks_passed=post.gap_root_checks_passed,
                before_phase_limit_change=pre.phase_limit_change,after_phase_limit_change=post.phase_limit_change,
                normal_mott_numerics_passed=passed,source_pair_proof_covers_endpoints=domain,
                full_levinson_accepted=false,production_authorized=false))
            println("[continuation] $(ch) q=$(q) T=[$(ts[left]),$(ts[right])] normal=$(passed) source_scope=$(domain)")
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(stage="normal_mott",channel=String(ch),q_inv_fm=q,reason=sprint(showerror,err)))
            println(stderr,"[continuation-failed] $(last(failures))")
        end
        checkpoint()
    end
    bg=R.saved_background(input,170.,24)
    for ch in R.CHANNELS,q in (0.01,1.4,3.2),thermal in (10.,24.)
        try
            k=kernel(bg,ch,q,thermal)
            S=last(k.edges)
            # At any fixed 0<q<2L the energy interval shrinks to zero at S.
            # This is not the q=0 discontinuous endpoint theorem.
            ps=[C.mapped_integral(x->k.rho(x)/(x-S)/pi,k.edges;nodes=n) for n in (64,128)]
            value=last(ps)
            ys=[k.rho(S-h) for h in (1e-3,1e-5,1e-7)]
            push!(finiteq,(channel=String(ch),q_inv_fm=q,thermal_cutoff_inv_fm=thermal,
                rho_1e3=ys[1],rho_1e5=ys[2],rho_1e7=ys[3],
                endpoint_inverse=real(1-4bg.coupling[ch]*value),
                pi_node_change_inv_fm2=abs(ps[1]-ps[2]),
                numerics_passed=abs(ps[1]-ps[2])<1e-6,
                q0_log_divergence_theorem_applies=false,full_UHP_certified=false,
                production_authorized=false))
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(stage="finite_q_edge",channel=String(ch),q_inv_fm=q,reason="Lth=$(thermal): "*sprint(showerror,err)))
        end
        checkpoint()
    end
    R.hashfile(bgpath)==bh || error("retained background drift")
    complete=isempty(failures) && length(mott)==8 && length(finiteq)==24
    R.finish_output(output,hashes,Dict("status"=>complete ? "conditional_normal_mott_and_finite_q_edges" : "evaluation_failed",
        "input_directory"=>input,"input_hashes"=>Dict("backgrounds.csv"=>bh),
        "normal_mott_rows"=>length(mott),"finite_q_rows"=>length(finiteq),
        "all_normal_mott_numerics_passed"=>length(mott)==8 && all(r.normal_mott_numerics_passed for r in mott),
        "source_pair_proof_covers_all_events"=>length(mott)==8 && all(r.source_pair_proof_covers_endpoints for r in mott),
        "finite_q_numerics_passed"=>length(finiteq)==24 && all(r.numerics_passed for r in finiteq),
        "step3_density_acceptance"=>"not_run_pending_endpoint_counting_and_source_scope",
        "solver_called"=>false,"meson_density_computed"=>false,"full_UHP_certified"=>false,
        "full_levinson_accepted"=>false,
        "limitations"=>["Retained p24 BQS backgrounds; no new solves or background interpolation",
            "Normal low-energy Mott checks do not include exterior regulator roots",
            "Negative finite-q tail is retained; positive endpoint inverse does not prove all exterior zeros absent",
            "Source-domain failure is a proof-coverage gap, not a proof that finite-density PNJL is invalid"] ))
    complete || error("continuation diagnostics incomplete; evidence retained")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
