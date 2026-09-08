"""Explicit quark-only backgrounds, Mott brackets, numerical axes and sparse freezeout.

Each invocation selects one stage. No production defaults or baselines are changed.
"""
module CausalGBUValidation
isdefined(Main,:Models) || Base.include(Main,normpath(joinpath(@__DIR__,"..","..","..","src","models","Models.jl")))
include("causal_gbu_research_utils.jl")
using .CausalGBUResearch
using CSV, JSON3
using Main.Models
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_coupling
const R = CausalGBUResearch
const HBARC = Main.Constants_PNJL.ħc_MeV_fm

function background(model,T_MeV,muB_MeV,pn;seed=nothing)
    kwargs = (p_num=pn,t_num=8,residual_norm_max=1e-7,iterations=300,xi=0.0)
    seed===nothing || (kwargs=merge(kwargs,(seed_guess=seed,)))
    t = time_ns()
    eq = Models.solve(model,Models.FixedMuBConservedCharges(muB_MeV/HBARC,0.4,0.0),T_MeV/HBARC;kwargs...)
    eq.converged || error("background failed T=$(T_MeV),muB=$(muB_MeV),residual=$(eq.residual_norm)")
    state = Models.meanfield_state(eq.x_state)
    kernel = build_full_kmt_interaction(state.phi;G=model.params.G_fm2,K=model.params.K_fm5)
    coupling = Dict(c=>charged_rpa_coupling(R.charged_rpa_spec(c),kernel) for c in R.CHANNELS)
    bg = (m=(u=Float64(eq.masses[1]),d=Float64(eq.masses[2]),s=Float64(eq.masses[3])),
        mu=(u=Float64(eq.mu_vec[1]),d=Float64(eq.mu_vec[2]),s=Float64(eq.mu_vec[3])),
        T=T_MeV/HBARC,Phi=Float64(state.Phi),PhiBar=Float64(state.PhiBar),coupling=coupling,
        vacuum=Float64(model.params.Λ_inv_fm),label="T$(T_MeV)_muB$(muB_MeV)_p$(pn)",
        T_MeV=Float64(T_MeV),muB_MeV=Float64(muB_MeV),residual=Float64(eq.residual_norm),solver_called=true)
    row = (T_MeV=Float64(T_MeV),muB_MeV=Float64(muB_MeV),p_nodes=pn,angle_nodes=8,
        m_u=bg.m.u,m_d=bg.m.d,m_s=bg.m.s,mu_u=bg.mu.u,mu_d=bg.mu.d,mu_s=bg.mu.s,
        phi_u=Float64(state.phi[1]),phi_d=Float64(state.phi[2]),phi_s=Float64(state.phi[3]),
        Phi=bg.Phi,PhiBar=bg.PhiBar,K12=coupling[:pi_plus],K45=coupling[:K_plus],
        residual=bg.residual,omega=eq.omega,elapsed_s=(time_ns()-t)/1e9,production_authorized=false)
    return (bg=bg,row=row,seed=Float64.(eq.solution))
end

function scan_mott(getbg,output)
    rows,events = NamedTuple[],NamedTuple[]
    # Offset is explicit and compared after refinement, not silently interpreted as exact threshold.
    function probe(T,channel,q,mesh)
        bg = getbg(T,240.0,24)
        s = R.Settings(mesh=mesh)
        b = R.bubble_at(bg,channel,q,s)
        return (value=real(b.inverse(b.threshold-1e-7)),bubble=b,bg=bg,settings=s)
    end
    for channel in R.CHANNELS, q in (0.0,0.5,1.0), mesh in (128,256)
        a,b = 170.0,230.0
        pa,pb = probe(a,channel,q,mesh),probe(b,channel,q,mesh)
        if !(pa.value<0<pb.value)
            push!(events,(channel=String(channel),q_inv_fm=q,mesh=mesh,T_lo_MeV=a,T_hi_MeV=b,
                before_count=-1,after_count=-1,phase_drop=NaN,threshold_value_lo=pa.value,
                threshold_value_hi=pb.value,passed=false,status="no_threshold_sign_bracket",production_authorized=false))
            continue
        end
        while b-a > 0.125
            mid=(a+b)/2
            pm=probe(mid,channel,q,mesh)
            if pm.value<0
                a,pa=mid,pm
            else
                b,pb=mid,pm
            end
        end
        endpoints = []
        # Audit separated points as well as the tight threshold bracket: near-Mott
        # roots can lie closer to the gap boundary than a fixed scan margin.
        for T in (a-0.5,b+0.5)
            p=probe(T,channel,q,mesh)
            g=R.gap_audit(p.bubble,p.settings)
            phase=-angle(p.bubble.inverse(p.bubble.threshold+1e-7))/pi
            push!(endpoints,(gap=g,phase=phase))
            push!(rows,(channel=String(channel),q_inv_fm=q,mesh=mesh,T_MeV=T,
                roots=g.count,contour_count=g.contour_count,passed=g.passed,phase_threshold_over_pi=phase,
                closest_gap_to_unitary=isempty(g.roots) ? NaN : minimum(r.distance_to_gap_upper for r in g.roots),
                threshold_margin_1e7=p.value,threshold_margin_1e9=real(p.bubble.inverse(p.bubble.threshold-1e-9)),
                production_authorized=false))
        end
        pre,post=endpoints
        drop=pre.phase-post.phase
        passed=pre.gap.passed && post.gap.passed && pre.gap.count-post.gap.count==1 && abs(drop-1)<0.005
        push!(events,(channel=String(channel),q_inv_fm=q,mesh=mesh,T_lo_MeV=a,T_hi_MeV=b,
            before_count=pre.gap.count,after_count=post.gap.count,phase_drop=drop,
            threshold_value_lo=pa.value,threshold_value_hi=pb.value,passed=passed,
            status=passed ? "conditional_thermal_mott_bracket" : "endpoint_count_unresolved",production_authorized=false))
        CSV.write(joinpath(output,"thermal_mott_brackets.csv"),events)
        CSV.write(joinpath(output,"thermal_mott_endpoints.csv"),rows)
        println("[gbu-mott] $(channel) q=$(q) mesh=$(mesh) T=[$a,$b] pass=$(passed)")
    end
    qrows=NamedTuple[]
    bg=getbg(170.0,240.0,24)
    for channel in R.CHANNELS, q in 0.0:0.25:4.0
        b=R.bubble_at(bg,channel,q,R.Settings())
        g=R.gap_audit(b,R.Settings())
        push!(qrows,(channel=String(channel),q_inv_fm=q,count=g.count,contour_count=g.contour_count,
            passed=g.passed,k0_root=isempty(g.roots) ? NaN : first(g.roots).omega_inv_fm,
            gap_to_unitary=isempty(g.roots) ? NaN : minimum(r.distance_to_gap_upper for r in g.roots),
            threshold_inverse=real(b.inverse(b.threshold-1e-7)),
            threshold_phase_over_pi=-angle(b.inverse(b.threshold+1e-7))/pi,
            event_interpretation="momentum_dissolution_regulator_dependent",production_authorized=false))
    end
    CSV.write(joinpath(output,"q_continuation.csv"),qrows)
    return Dict("thermal_passed"=>all(r.passed for r in events),"q_count_passed"=>all(r.passed for r in qrows),
        "thermal_scope"=>"muB=240,background p24,threshold offset 1e-7; mesh128/256; not full spectrum certification")
end

function numerical_axes(bg,output)
    # One-axis changes at a time; model sensitivity is not included in numerical acceptance.
    base=(mesh=128,nw=1200,nq=12,qmax=8.0,thermal=20.0,np=128,nx=64)
    variants=[("base",(;)),("mesh256",(mesh=256,)),("omega2400",(nw=2400,)),
        ("q16",(nq=16,)),("qmax10",(qmax=10.0,)),("upper56",(upper=56.0,)),
        ("lower1e5",(lower=1e-5,)),("thermal24",(thermal=24.0,)),
        ("model_all_hard",(thermal=bg.vacuum,))]
    rows,shells=NamedTuple[],NamedTuple[]
    for (name,change) in variants, channel in R.CHANNELS
        s=R.Settings(;merge(base,change)...)
        d=R.density(bg,channel,s)
        append!(shells,[merge((variant=name,),r) for r in d.rows])
        push!(rows,(variant=name,channel=String(channel),density_inv_fm3=d.density,passed=d.passed,
            mesh=s.mesh,omega_nodes=s.nw,q_nodes=s.nq,qmax=s.qmax,thermal_max=s.thermal,
            omega_min=s.lower,omega_max=s.upper,momentum_nodes=s.np,angle_nodes=s.nx,
            interpretation=name=="model_all_hard" ? "model_sensitivity_not_upstream_matched" : "numerical_axis",
            production_authorized=false))
        CSV.write(joinpath(output,"axis_densities.csv"),rows)
        CSV.write(joinpath(output,"axis_shells.csv"),shells)
        println("[gbu-axis] $(name) $(channel) n=$(d.density) pass=$(d.passed)")
    end
    return Dict("all_conditional_gates"=>all(r.passed for r in rows if r.interpretation=="numerical_axis"),
        "background"=>bg,"spectral_energy_nodes"=>64)
end

function numerical_refinement(bg,output)
    rows,shells=NamedTuple[],NamedTuple[]
    for (name,change) in (("energy128",(ne=128,)),("q24",(nq=24,)),
                           ("q_tail_8_10",(qmax=10.0,ne=128,))), channel in R.CHANNELS
        s=R.Settings(;change...)
        lower=name=="q_tail_8_10" ? 8.0 : 0.0
        q,weights=R.gauleg(lower,s.qmax,s.nq)
        parts=[R.shell(bg,channel,v,s) for v in q]
        total=sum(weights[i]*parts[i].shell_inv_fm2 for i in eachindex(q))
        # A separated tail is signed; only the full partial yield has a positivity gate.
        passed=all(r.passed for r in parts) && isfinite(total) && (lower>0 || total>=0)
        append!(shells,[merge((variant=name,),r) for r in parts])
        push!(rows,(variant=name,channel=String(channel),density_inv_fm3=total,passed=passed,
            mesh=s.mesh,energy_nodes=s.ne,omega_nodes=s.nw,q_nodes=s.nq,qmin=lower,qmax=s.qmax,
            interpretation=lower>0 ? "signed_tail_band" : "numerical_axis",production_authorized=false))
        CSV.write(joinpath(output,"refined_axes.csv"),rows)
        CSV.write(joinpath(output,"refined_shells.csv"),shells)
        println("[gbu-refine] $(name) $(channel) n=$(total) pass=$(passed)")
    end
    return Dict("all_conditional_gates"=>all(r.passed for r in rows),"background"=>bg,
        "reference_base"=>"Settings defaults equal original method_v1_axes/base; ne=64 was implicit")
end

function sparse_failure_refinement(input,output)
    files=("manifest.json","backgrounds.csv","sparse_densities.csv","sparse_shells.csv")
    hashes=Dict(f=>R.hashfile(joinpath(input,f)) for f in files)
    failed=filter(r->!r.passed,collect(CSV.File(joinpath(input,"sparse_shells.csv"))))
    densities=collect(CSV.File(joinpath(input,"sparse_densities.csv")))
    rows=NamedTuple[]
    backgrounds=Dict{String,Any}()
    for old in failed
        parent=only(filter(r->r.sqrt_s_NN_GeV==old.sqrt_s_NN_GeV &&
            r.background_p==old.background_p && r.channel==old.channel,densities))
        bg=R.saved_background(input,parent.T_MeV,old.background_p)
        backgrounds[bg.label]=bg
        base=(mesh=256,nw=2400,nq=16)
        for (variant,change) in (("contour_only",(;)),("energy128",(ne=128,)),
                                ("mesh512",(mesh=512,)),("omega4800",(nw=4800,)))
            s=R.Settings(;merge(base,change)...)
            r=R.shell(bg,Symbol(old.channel),old.q_inv_fm,s)
            push!(rows,merge((sqrt_s_NN_GeV=old.sqrt_s_NN_GeV,background_p=old.background_p,
                variant=variant,mesh=s.mesh,energy_nodes=s.ne,omega_nodes=s.nw,
                previous_shell=old.shell_inv_fm2,previous_contour_passed=old.contour_passed),r))
            CSV.write(joinpath(output,"failed_shell_refinement.csv"),rows)
            println("[gbu-sparse-refine] $(old.channel) p=$(old.background_p) $(variant) pass=$(r.passed)")
        end
    end
    all(R.hashfile(joinpath(input,f))==h for (f,h) in hashes) || error("sparse inputs changed")
    return Dict("input_directory"=>abspath(input),"input_hashes"=>hashes,"backgrounds"=>backgrounds,
        "previous_failed_shells"=>length(failed),"all_targeted_gates"=>!isempty(rows) && all(r.passed for r in rows),
        "scope"=>"failed shells only, not a new complete sparse acceptance or global spectrum count",
        "complete_curve_authorized"=>false)
end

function saved_sparse_refinement(input,output)
    files=("manifest.json","backgrounds.csv","sparse_densities.csv","sparse_shells.csv")
    hashes=Dict(f=>R.hashfile(joinpath(input,f)) for f in files)
    prior=filter(r->r.background_p==48,collect(CSV.File(joinpath(input,"sparse_densities.csv"))))
    length(prior)==12 || error("expected three four-channel backgrounds at p48")
    s=R.Settings(mesh=512,ne=128,nw=4800,nq=24,lower=1e-5,upper=56.0,thermal=24.0)
    rows,shells,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    target=R.method_contract()["gates"]["relative_density_target"]
    for old in prior
        bg=R.saved_background(input,old.T_MeV,old.background_p)
        d=try
            R.density(bg,Symbol(old.channel),s)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(sqrt_s_NN_GeV=old.sqrt_s_NN_GeV,channel=String(old.channel),
                stage="saved_spectrum",reason=sprint(showerror,err)))
            CSV.write(joinpath(output,"refined_failures.csv"),failures)
            continue
        end
        change=d.density/old.density_inv_fm3-1
        append!(shells,[merge((sqrt_s_NN_GeV=old.sqrt_s_NN_GeV,background_p=48,),r) for r in d.rows])
        push!(rows,(sqrt_s_NN_GeV=old.sqrt_s_NN_GeV,T_MeV=old.T_MeV,muB_MeV=old.muB_MeV,
            background_p=48,channel=String(old.channel),density_inv_fm3=d.density,
            previous_density=old.density_inv_fm3,previous_gate_passed=old.passed,
            relative_change=change,numerical_change_passed=isfinite(change) && abs(change)<target,
            passed=d.passed,failed_shells=count(r->!r.passed,d.rows),production_authorized=false))
        CSV.write(joinpath(output,"refined_densities.csv"),rows)
        CSV.write(joinpath(output,"refined_shells.csv"),shells)
        println("[gbu-saved-sparse] s=$(old.sqrt_s_NN_GeV) $(old.channel) n=$(d.density) change=$(change) pass=$(d.passed)")
    end
    all(R.hashfile(joinpath(input,f))==h for (f,h) in hashes) || error("sparse inputs changed")
    complete=length(rows)==12 && isempty(failures)
    return Dict("settings"=>s,"input_directory"=>abspath(input),"input_hashes"=>hashes,
        "all_refined_conditional_gates"=>complete && all(r.passed for r in rows),
        "density_change_target_passed"=>complete && all(r.numerical_change_passed for r in rows),
        "failed_evaluations"=>length(failures),"complete_curve_authorized"=>false,
        "scope"=>"joint refinement on retained p48 backgrounds; no full-spectrum or regulator independence claim")
end

function momentum_events(getbg,output)
    bg=getbg(170.0,240.0,24)
    rows=NamedTuple[]
    for channel in R.CHANNELS, mesh in (128,256)
        s=R.Settings(mesh=mesh)
        probe(q)=begin
            b=R.bubble_at(bg,channel,q,s)
            (value=real(b.inverse(b.threshold-1e-7)),bubble=b)
        end
        lo,hi=0.0,6.0
        pl,ph=probe(lo),probe(hi)
        pl.value<0<ph.value || error("no momentum threshold bracket for $(channel)")
        while hi-lo>0.001
            mid=(lo+hi)/2
            p=probe(mid)
            if p.value<0
                lo,pl=mid,p
            else
                hi,ph=mid,p
            end
        end
        before,after=probe(lo-0.1),probe(hi+0.1)
        gb,ga=R.gap_audit(before.bubble,s),R.gap_audit(after.bubble,s)
        phase(p)=-angle(p.bubble.inverse(p.bubble.threshold+1e-7))/pi
        drop=phase(before)-phase(after)
        push!(rows,(channel=String(channel),mesh=mesh,q_lo=lo,q_hi=hi,
            before_count=gb.count,after_count=ga.count,before_contour=gb.contour_count,after_contour=ga.contour_count,
            phase_drop=drop,passed=gb.passed && ga.passed && gb.count-ga.count==1 && abs(drop-1)<0.005,
            interpretation="regulator_dependent_momentum_dissolution",production_authorized=false))
        CSV.write(joinpath(output,"momentum_brackets.csv"),rows)
        println("[gbu-q-event] $(channel) mesh=$(mesh) q=[$lo,$hi]")
    end
    b=R.bubble_at(bg,:K_minus,3.5,R.Settings())
    g=R.gap_audit(b,R.Settings())
    return Dict("momentum_events_passed"=>all(r.passed for r in rows),
        "previous_q3p5_resolution"=>g,"scope"=>"conditional kinematic gaps, not all cutoff-created gaps")
end

function sparse_freezeout(getbg,output,s)
    points=Models.build_freezeout_scan_points([3.0,7.7,200.0];profile=Models.load_freezeout_profile(profile="default"),traversal=:sqrts_descending)
    rows,shells=NamedTuple[],NamedTuple[]
    failures=NamedTuple[]
    for pt in points, pn in (24,48)
        bg=try
            getbg(pt.T_MeV,pt.muB_MeV,pn)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,background_p=pn,channel="all",
                stage="background",reason=sprint(showerror,err)))
            CSV.write(joinpath(output,"sparse_failures.csv"),failures)
            continue
        end
        for channel in R.CHANNELS
            d=try
                R.density(bg,channel,s)
            catch err
                err isa InterruptException && rethrow()
                push!(failures,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,background_p=pn,channel=String(channel),
                    stage="spectrum",reason=sprint(showerror,err)))
                CSV.write(joinpath(output,"sparse_failures.csv"),failures)
                continue
            end
            append!(shells,[merge((sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,background_p=pn,),r) for r in d.rows])
            push!(rows,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,T_MeV=pt.T_MeV,muB_MeV=pt.muB_MeV,
                background_p=pn,channel=String(channel),density_inv_fm3=d.density,passed=d.passed,
                failed_shells=count(r->!r.passed,d.rows),production_authorized=false))
            CSV.write(joinpath(output,"sparse_densities.csv"),rows)
            CSV.write(joinpath(output,"sparse_shells.csv"),shells)
            println("[gbu-freezeout] s=$(pt.sqrt_s_NN_GeV) p=$(pn) $(channel) n=$(d.density) pass=$(d.passed)")
        end
    end
    return Dict("conditional_gates_passed"=>length(rows)==24 && isempty(failures) && all(r.passed for r in rows),
        "failed_evaluations"=>length(failures),
        "complete_curve_authorized"=>false,"experimental_comparison_performed"=>false)
end

function main()
    R.method_contract()
    stage=get(ENV,"GBU_RESEARCH_STAGE","mott")
    stage in ("backgrounds","mott","q_refine","axes","axes_refine","freezeout","sparse_refine","saved_sparse") || error("unknown stage")
    sparse=stage=="freezeout" ? R.sparse_settings() : nothing
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    output=get(ENV,"GBU_RESEARCH_OUTPUT",joinpath(base,"method_v1_$(stage)"))
    hashes=R.start_output(output)
    model=Models.create_model(:PNJL)
    cache=Dict{Tuple{Float64,Float64,Int},Any}()
    bgs=NamedTuple[]
    function getbg(T,muB,pn)
        key=(Float64(T),Float64(muB),pn)
        if !haskey(cache,key)
            candidate=isempty(cache) ? nothing : argmin(k->abs(k[1]-T)+abs(k[2]-muB),collect(keys(cache)))
            seed=candidate===nothing ? nothing : cache[candidate].seed
            r=background(model,T,muB,pn;seed=seed)
            cache[key]=r
            push!(bgs,r.row)
            CSV.write(joinpath(output,"backgrounds.csv"),bgs)
            println("[gbu-background] T=$(T) muB=$(muB) p=$(pn) residual=$(r.bg.residual)")
        end
        return cache[key].bg
    end
    record=Dict{String,Any}("stage"=>stage,"solver_called"=>!(stage in ("axes","axes_refine","sparse_refine","saved_sparse")),"status"=>"started",
        "sparse_settings"=>sparse)
    try
        result = if stage=="backgrounds"
            for pn in (8,24,48,96)
                getbg(170.0,240.0,pn)
            end
            Dict("note"=>"background quadrature comparison only")
        elseif stage=="mott"
            scan_mott(getbg,output)
        elseif stage=="axes"
            numerical_axes(R.frozen_background(joinpath(base,"negative_density_phase_fig2_like")),output)
        elseif stage=="axes_refine"
            numerical_refinement(R.frozen_background(joinpath(base,"negative_density_phase_fig2_like")),output)
        elseif stage=="sparse_refine"
            sparse_failure_refinement(get(ENV,"GBU_RESEARCH_INPUT",joinpath(base,"method_v1_freezeout_v3")),output)
        elseif stage=="saved_sparse"
            saved_sparse_refinement(get(ENV,"GBU_RESEARCH_INPUT",joinpath(base,"method_v1_freezeout_v3")),output)
        elseif stage=="q_refine"
            momentum_events(getbg,output)
        else
            sparse_freezeout(getbg,output,sparse)
        end
        merge!(record,result)
        record["status"]="diagnostic_stage_completed"
    catch err
        err isa InterruptException && rethrow()
        record["status"]="failed_retained"
        record["failure"]=sprint(showerror,err,catch_backtrace())
        println(stderr,record["failure"])
    end
    record["background_solves"]=length(cache)
    R.finish_output(output,hashes,record)
    record["status"]=="failed_retained" && error("stage failed; diagnostics preserved at $(output)")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
