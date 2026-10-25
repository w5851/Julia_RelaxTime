#!/usr/bin/env julia
# Thermodynamic-only qualification of the reanchored mode-A production plan.
using CSV, JSON, SHA
include(joinpath(@__DIR__,"..","..","relaxtime","run_phase_guided_transport_scan.jl"))
include(joinpath(@__DIR__,"ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon: start_run, finish_run, write_rows, write_json

function fixed_temperature_boundary(T_MeV,muB_MeV,lower,upper,numerics;seed_pair=nothing)
    phase=GapTransportScanPhaseEquilibrium
    function solve(x)
        seeds=seed_pair===nothing ? NamedTuple() :
            (hadron_seed=seed_pair.high_mass_hadron.x_state,quark_seed=seed_pair.low_mass_quark.x_state)
        result=phase.solve_two_branch_candidates(T_MeV/ħc_MeV_fm,muB_MeV/(3ħc_MeV_fm),x,numerics;seeds...)
        seed_pair=result
        return result
    end
    a,b=Float64(lower),Float64(upper)
    left,right=solve(a),solve(b)
    left.distinct && right.distinct && left.delta_omega_high_minus_low*right.delta_omega_high_minus_low<=0 ||
        error("fixed-T first-order boundary has no two-branch sign bracket")
    for _ in 1:24
        x=(a+b)/2;mid=solve(x)
        mid.distinct || error("fixed-T phase boundary lost distinct branches")
        if left.delta_omega_high_minus_low*mid.delta_omega_high_minus_low<=0
            b=x;right=mid
        else
            a=x;left=mid
        end
    end
    return (xi=(a+b)/2,lower=a,upper=b)
end

function reuse_phase_anchors(source,inputdir,outdir,profile,opts)
    manifest_path=joinpath(source,"manifest.json")
    manifest=JSON.parsefile(manifest_path)
    manifest["parameters"]["profile"]==JSON.parse(JSON.json(Models.TransportSamplingProfiles.effective_config(profile))) ||
        error("anchor sampling profile differs from the current profile")
    # The scan below changes only how already certified anchors are consumed.
    # Require every other snapshotted source/config byte to match before reuse.
    current_script=relpath(@__FILE__,PROJECT_ROOT)
    for entry in manifest["source_files"]
        relative=entry["path"]
        frozen=joinpath(source,"source",relative)
        bytes2hex(sha256(read(frozen)))==entry["sha256"] || error("anchor snapshot mismatch: $relative")
        normpath(relative)==normpath(current_script) && continue
        bytes2hex(sha256(read(joinpath(PROJECT_ROOT,relative))))==entry["sha256"] ||
            error("anchor source changed: $relative; requalify before reusing anchors")
    end
    for entry in manifest["input_files"]
        bytes2hex(sha256(read(joinpath(inputdir,basename(entry["path"])))))==entry["sha256"] ||
            error("anchor input identity changed")
    end
    anchors=[(;pairs(row)...) for row in CSV.File(joinpath(source,"anchors.csv"))]
    Set(a.muB_MeV for a in anchors)==Set(opts.muB_values) || error("incomplete anchor set")
    length(anchors)==length(opts.muB_values) && all(a->a.certified &&
        isfinite(a.refinement_delta_MeV) && abs(a.refinement_delta_MeV)<=0.1,anchors) ||
        error("uncertified or duplicate anchor")
    Point=PhaseGuidedPlan.PhaseGuidedPoint
    points=Point[]
    for row in CSV.File(joinpath(source,"sampling_plan.csv"))
        values=map(zip(fieldnames(Point),fieldtypes(Point))) do (name,type)
            value=getproperty(row,name)
            type===Symbol ? Symbol(value) : convert(type,value)
        end
        push!(points,Point(values...))
    end
    length(points)==910 && length(unique((p.T_MeV,p.muB_MeV,p.xi) for p in points))==910 ||
        error("incomplete or duplicate mode-A plan")
    Set(p.alpha_T for p in points)==Set(opts.alpha_T_values) || error("anchor plan alpha_T differs")
    all(p->p.anchor_convergence_certified && any(a->a.muB_MeV==p.muB_MeV &&
        a.T_MeV==p.T_phase_base_MeV,anchors),points) || error("plan and anchor table disagree")
    copied=NamedTuple[]
    for name in ("anchors.csv","sampling_plan.csv","manifest.json")
        destination=joinpath(outdir,"anchor_reuse",name)
        mkpath(dirname(destination));cp(joinpath(source,name),destination)
        push!(copied,(name=name,sha256=bytes2hex(sha256(read(destination)))))
    end
    write_json(joinpath(outdir,"anchor_reuse.json"),(source=abspath(source),files=copied,
        all_other_source_bytes_identical=true,excluded_current_source=current_script))
    return anchors,points
end

function trace_phase_direction(points,numerics;reverse_order=false,label="",trace_path=nothing)
    phase=GapTransportScanPhaseEquilibrium
    results=Vector{Any}(undef,length(points))
    seed=reverse_order ? Models.HADRON_SEED_5 : Models.QUARK_SEED_5
    indices=reverse_order ? reverse(eachindex(points)) : eachindex(points)
    rows=NamedTuple[]
    for (count,i) in enumerate(indices)
        point=points[i]
        eq=phase.solve_phase_branch(point.T_MeV/ħc_MeV_fm,point.muB_MeV/(3ħc_MeV_fm),point.xi,seed,numerics)
        eq===nothing && error("phase continuation failed at $(point.scan_group) xi=$(point.xi)")
        results[i]=eq;seed=eq.x_state
        push!(rows,(T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,xi=point.xi,
            phi_u=seed[1],phi_d=seed[2],phi_s=seed[3],Phi=seed[4],Phibar=seed[5],
            m_u=eq.masses[1],m_d=eq.masses[2],m_s=eq.masses[3],omega=eq.omega,
            residual_norm=eq.residual_norm))
        if count%25==0 || count==length(points)
            trace_path===nothing || write_rows(trace_path,rows)
            println("TRACE ",label," reverse=",reverse_order," ",count,"/",length(points));flush(stdout)
        end
    end
    return results
end

function combine_phase_candidates(candidate_a,candidate_b)
    high,low=candidate_a.masses[1]>=candidate_b.masses[1] ? (candidate_a,candidate_b) : (candidate_b,candidate_a)
    delta=high.omega-low.omega
    distinct=GapTransportScanPhaseEquilibrium._relative_delta(high.masses[1],low.masses[1])>1e-3
    return (distinct=distinct,high_mass_hadron=high,low_mass_quark=low,
        delta_omega_high_minus_low=delta,stable=delta<=0 ? high : low,
        stable_phase=delta<=0 ? :hadron : :quark)
end

function phase_revalidation(inputdir,outdir;anchor_source=nothing)
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    args=["--mode","fixed-muB-phase-scaled","--sampling-profile",profile.source_path,
        "--muB-list","0,450,900","--alphaT-list","1,1.1,1.2",
        "--xi-list",join(collect(-50:50)./100,','),"--outdir",outdir,"--dry-run"]
    opts=PhaseGuidedCLI.parse_args(args)
    reference=_load_runtime_phase_reference(opts)
    base=(p_num=profile.background.p_nodes,t_num=profile.background.cos_nodes,
        model=Models.TransportSamplingProfiles.profile_model(profile),residual_norm_max=1e-8)
    refined=(p_num=128,t_num=32,
        model=Models.TransportSamplingProfiles.profile_model(profile;thermal_p_max_inv_fm=30.0),residual_norm_max=1e-8)
    sources=filter(p->endswith(p,".jl"),readdir(joinpath(PROJECT_ROOT,"scripts","relaxtime");join=true))
    start_run(outdir,inputdir,@__FILE__,(purpose="phase anchors and stable branches under the new background",
        profile=Models.TransportSamplingProfiles.effective_config(profile),muB_MeV=opts.muB_values,
        alpha_T=opts.alpha_T_values,xi=opts.xi_values,anchor_tolerance_MeV=0.1,background_relative_target=0.01,
        reference=(thermal_p_max_inv_fm=30.0,p_nodes=128,cos_nodes=32),
        branch_strategy="opposite-end continuation with independent cold-seed midpoint audits",
        branch_iterations_per_attempt=GapTransportScanPhaseEquilibrium.gap_iteration_limit(base),
        anchor_source=anchor_source,
        scope="mode-A plan only; not a global phase diagram"),extra_sources=sources)
    phase=GapTransportScanPhaseEquilibrium
    anchors=NamedTuple[];points=PhaseGuidedPlan.PhaseGuidedPoint[]
    if anchor_source!==nothing
        anchors,points=reuse_phase_anchors(anchor_source,inputdir,outdir,profile,opts)
        write_rows(joinpath(outdir,"anchors.csv"),anchors)
        println("REUSED ",length(anchors)," certified anchors and ",length(points)," plan points");flush(stdout)
    else
      for mu in opts.muB_values
        println("ANCHOR START muB=",mu);flush(stdout)
        one=PhaseGuidedCLI.parse_args(vcat(args,["--muB-list",string(mu)]))
        measured=@timed PhaseGuidedPlan.build_plan(one;phase_reference=reference)
        plan=measured.value
        append!(points,plan.points)
        firstpoint=first(plan.points)
        oldkind,oldT=PhaseGuidedPlan._phase_reference_for_mode_a(mu,0.0;phase_reference=reference)
        push!(anchors,(muB_MeV=mu,kind=String(firstpoint.phase_reference_kind),
            method=String(firstpoint.phase_anchor_method),old_reference_T_MeV=oldT,
            T_MeV=firstpoint.T_phase_base_MeV,shift_MeV=firstpoint.T_phase_base_MeV-oldT,
            refinement_delta_MeV=firstpoint.anchor_convergence_delta_MeV,
            certified=firstpoint.anchor_convergence_certified,side_delta_xi=firstpoint.coexistence_delta_xi,
            seconds=measured.time,gc_seconds=measured.gctime,allocated_bytes=measured.bytes))
        write_rows(joinpath(outdir,"anchors.csv"),anchors)
        println("ANCHOR DONE muB=",mu," T=",firstpoint.T_phase_base_MeV,
            " refinement=",firstpoint.anchor_convergence_delta_MeV);flush(stdout)
      end
    end
    plan=PhaseGuidedPlan.PhaseGuidedPlan(points,length(points))
    PhaseGuidedAssets.write_plan_csv(joinpath(outdir,"sampling_plan.csv"),plan)
    rows=NamedTuple[];failures=NamedTuple[];boundaries=NamedTuple[];audits=NamedTuple[]
    for group in unique(p.scan_group for p in points)
      group_points=sort(filter(p->p.scan_group==group,points);by=p->p.xi)
      println("GROUP START ",group);flush(stdout)
      try
        traced=map((("candidate",base),("refined",refined))) do (label,numerics)
            forward=trace_phase_direction(group_points,numerics;label="$group $label",
                trace_path=joinpath(outdir,"traces","$(group)_$(label)_forward.csv"))
            backward=trace_phase_direction(group_points,numerics;reverse_order=true,label="$group $label",
                trace_path=joinpath(outdir,"traces","$(group)_$(label)_backward.csv"))
            candidates=combine_phase_candidates.(forward,backward)
            i=cld(length(group_points),2);point=group_points[i]
            measured=@timed phase.solve_two_branch_candidates(point.T_MeV/ħc_MeV_fm,
                point.muB_MeV/(3ħc_MeV_fm),point.xi,numerics;branch_solver=phase.solve_models_equilibrium)
            cold=measured.value
            mass_delta=maximum(abs.(candidates[i].stable.masses./cold.stable.masses.-1))
            omega_delta=abs(candidates[i].stable.omega-cold.stable.omega)
            passed=mass_delta<=1e-6 && omega_delta<=1e-8
            push!(audits,(scan_group=group,level=label,T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,
                xi=point.xi,max_mass_relative_change=mass_delta,omega_absolute_change=omega_delta,
                seconds=measured.time,passed=passed))
            write_rows(joinpath(outdir,"cold_seed_audits.csv"),audits)
            passed || error("cold-seed phase audit failed for $group $label")
            return candidates
        end
        previous=nothing
        for (i,point) in enumerate(group_points)
            candidates,check=traced[1][i],traced[2][i]
            ea,eb=candidates.stable,check.stable
            sa=Models.model_thermo(base.model,ea.x_state,ea.mu_vec,point.T_MeV/ħc_MeV_fm;
                xi=point.xi,p_num=base.p_num,t_num=base.t_num)[3]
            sb=Models.model_thermo(refined.model,eb.x_state,eb.mu_vec,point.T_MeV/ħc_MeV_fm;
                xi=point.xi,p_num=refined.p_num,t_num=refined.t_num)[3]
            mass_delta=maximum(abs.(ea.masses./eb.masses.-1))
            entropy_delta=sa/sb-1
            same_phase=!(candidates.distinct && check.distinct) || candidates.stable_phase===check.stable_phase
            passed=same_phase && isfinite(entropy_delta) && mass_delta<=0.01 && abs(entropy_delta)<=0.01
            push!(rows,(T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,alpha_T=point.alpha_T,xi=point.xi,
                distinct_branches=candidates.distinct,refined_distinct_branches=check.distinct,
                delta_omega=candidates.delta_omega_high_minus_low,refined_delta_omega=check.delta_omega_high_minus_low,
                stable_phase=String(candidates.stable_phase),refined_stable_phase=String(check.stable_phase),
                m_u=ea.masses[1],m_s=ea.masses[3],entropy=sa,
                max_mass_relative_change=mass_delta,entropy_relative_change=entropy_delta,passed=passed))
            if candidates.distinct
                prev=previous
                if prev!==nothing && prev.delta_omega*candidates.delta_omega_high_minus_low<0
                    a=fixed_temperature_boundary(point.T_MeV,point.muB_MeV,prev.xi,point.xi,base;seed_pair=candidates)
                    b=fixed_temperature_boundary(point.T_MeV,point.muB_MeV,prev.xi,point.xi,refined;seed_pair=check)
                    push!(boundaries,(T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,alpha_T=point.alpha_T,
                        xi=a.xi,refined_xi=b.xi,refinement_delta_xi=b.xi-a.xi,
                        bracket_lower_xi=prev.xi,bracket_upper_xi=point.xi))
                    write_rows(joinpath(outdir,"phase_boundaries.csv"),boundaries)
                end
                previous=(xi=point.xi,delta_omega=candidates.delta_omega_high_minus_low)
            end
        end
      catch err
            err isa InterruptException && rethrow()
            push!(failures,(scan_group=group,error=sprint(showerror,err)))
            write_rows(joinpath(outdir,"failures.csv"),failures)
      end
            write_rows(joinpath(outdir,"background.csv"),rows)
            write_json(joinpath(outdir,"progress.json"),(completed=length(rows),total=length(points),failures=length(failures)))
            println("BACKGROUND ",length(rows),"/",length(points)," failures=",length(failures));flush(stdout)
    end
    passed=isempty(failures) && length(rows)==length(points) && all(r->r.passed,rows) && all(a->a.certified,anchors)
    finish_run(outdir;points=length(points),completed_points=length(rows),failures=length(failures),
        phase_plan_passed=passed,first_order_boundaries=length(boundaries),
        cold_seed_audits=length(audits),cold_seed_audits_passed=all(a->a.passed,audits))
    passed || error("new-background phase plan failed; inspect retained evidence")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (2,3) || error("Usage: run_transport_phase_revalidation.jl INPUT_DIR NEW_OUTPUT_DIR [ANCHOR_SOURCE]")
    phase_revalidation(abspath(ARGS[1]),abspath(ARGS[2]);anchor_source=length(ARGS)==3 ? abspath(ARGS[3]) : nothing)
end
