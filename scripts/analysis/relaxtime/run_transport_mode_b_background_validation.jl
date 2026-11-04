#!/usr/bin/env julia
# Qualify the fixed-temperature mode-B grid without treating mode-A anchors as input.
include(joinpath(@__DIR__,"run_transport_phase_revalidation.jl"))

function mode_b_audit_indices(traced,n)
    selected=Int[1,cld(n,4),cld(n,2),cld(3n,4),n]
    for candidates in traced, i in 2:n
        left,right=candidates[i-1],candidates[i]
        if left.distinct && right.distinct &&
            left.delta_omega_high_minus_low*right.delta_omega_high_minus_low<=0
            append!(selected,(i-1,i))
        end
    end
    return sort!(unique!(selected))
end

function mode_b_background_validation(inputdir,outdir)
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    args=["--mode","fixed-T-sparse-muB","--sampling-profile",profile.source_path,
        "--T-list","120,160,200","--muB-list","0,450,900",
        "--xi-list",join(collect(-50:50)./100,','),"--outdir",outdir,"--dry-run"]
    opts=PhaseGuidedCLI.parse_args(args)
    reference=_load_runtime_phase_reference(opts)
    plan=PhaseGuidedPlan.build_plan(opts;phase_reference=reference)
    length(plan.points)==909 || error("incomplete mode-B plan")
    base=(p_num=profile.background.p_nodes,t_num=profile.background.cos_nodes,
        model=Models.TransportSamplingProfiles.profile_model(profile),residual_norm_max=1e-8)
    refined=(p_num=128,t_num=32,
        model=Models.TransportSamplingProfiles.profile_model(profile;thermal_p_max_inv_fm=30.0),residual_norm_max=1e-8)
    sources=filter(p->endswith(p,".jl"),readdir(joinpath(PROJECT_ROOT,"scripts","relaxtime");join=true))
    push!(sources,joinpath(@__DIR__,"run_transport_phase_revalidation.jl"))
    start_run(outdir,inputdir,@__FILE__,(purpose="mode-B stable backgrounds on the requested fixed-temperature grid",
        mode="B",profile=Models.TransportSamplingProfiles.effective_config(profile),
        T_MeV=opts.T_values,muB_MeV=opts.muB_values,xi=opts.xi_values,
        units=(T_MeV="MeV",muB_MeV="MeV",xi="dimensionless",masses="fm^-1",entropy="fm^-3",omega="fm^-4"),
        background_relative_target=0.01,reference=(thermal_p_max_inv_fm=30.0,p_nodes=128,cos_nodes=32),
        branch_strategy="opposite-end continuation; independent cold seeds at five evenly spaced positions and sign-bracket endpoints",
        cold_seed_mass_relative_tolerance=1e-6,cold_seed_omega_absolute_tolerance=1e-8,
        scope="fixed mode-B grid backgrounds only; transport and global phase diagram require separate qualification"),
        extra_sources=sources)
    PhaseGuidedAssets.write_plan_csv(joinpath(outdir,"sampling_plan.csv"),plan)
    phase=GapTransportScanPhaseEquilibrium
    rows=NamedTuple[];failures=NamedTuple[];boundaries=NamedTuple[];audits=NamedTuple[]
    for group in unique(p.scan_group for p in plan.points)
        points=sort(filter(p->p.scan_group==group,plan.points);by=p->p.xi)
        println("GROUP START ",group);flush(stdout)
        try
            traced=map((("candidate",base),("refined",refined))) do (label,numerics)
                forward=trace_phase_direction(points,numerics;label="$group $label",
                    trace_path=joinpath(outdir,"traces","$(group)_$(label)_forward.csv"))
                backward=trace_phase_direction(points,numerics;reverse_order=true,label="$group $label",
                    trace_path=joinpath(outdir,"traces","$(group)_$(label)_backward.csv"))
                combine_phase_candidates.(forward,backward)
            end
            selected=mode_b_audit_indices(traced,length(points))
            for (level,(label,numerics)) in enumerate((("candidate",base),("refined",refined))), i in selected
                point=points[i]
                measured=@timed phase.solve_two_branch_candidates(point.T_MeV/ħc_MeV_fm,
                    point.muB_MeV/(3ħc_MeV_fm),point.xi,numerics;branch_solver=phase.solve_models_equilibrium)
                cold=measured.value
                mass_delta=maximum(abs.(traced[level][i].stable.masses./cold.stable.masses.-1))
                omega_delta=abs(traced[level][i].stable.omega-cold.stable.omega)
                passed=mass_delta<=1e-6 && omega_delta<=1e-8
                push!(audits,(scan_group=group,level=label,T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,
                    xi=point.xi,max_mass_relative_change=mass_delta,omega_absolute_change=omega_delta,
                    seconds=measured.time,passed=passed))
                write_rows(joinpath(outdir,"cold_seed_audits.csv"),audits)
                passed || error("cold-seed phase audit failed for $group $label xi=$(point.xi)")
            end
            for (i,point) in enumerate(points)
                candidate,check=traced[1][i],traced[2][i]
                ea,eb=candidate.stable,check.stable
                sa=Models.model_thermo(base.model,ea.x_state,ea.mu_vec,point.T_MeV/ħc_MeV_fm;
                    xi=point.xi,p_num=base.p_num,t_num=base.t_num)[3]
                sb=Models.model_thermo(refined.model,eb.x_state,eb.mu_vec,point.T_MeV/ħc_MeV_fm;
                    xi=point.xi,p_num=refined.p_num,t_num=refined.t_num)[3]
                mass_delta=maximum(abs.(ea.masses./eb.masses.-1))
                entropy_delta=sa/sb-1
                same_phase=!(candidate.distinct && check.distinct) || candidate.stable_phase===check.stable_phase
                resolved=!(candidate.distinct && candidate.delta_omega_high_minus_low==0) &&
                    !(check.distinct && check.delta_omega_high_minus_low==0)
                passed=same_phase && resolved && isfinite(entropy_delta) && sa>0 && sb>0 &&
                    mass_delta<=0.01 && abs(entropy_delta)<=0.01
                push!(rows,(scan_group=group,T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,xi=point.xi,
                    distinct_branches=candidate.distinct,refined_distinct_branches=check.distinct,
                    delta_omega=candidate.delta_omega_high_minus_low,refined_delta_omega=check.delta_omega_high_minus_low,
                    stable_phase=String(candidate.stable_phase),refined_stable_phase=String(check.stable_phase),
                    coexistence_resolved=resolved,m_u=ea.masses[1],m_d=ea.masses[2],m_s=ea.masses[3],entropy=sa,
                    max_mass_relative_change=mass_delta,entropy_relative_change=entropy_delta,passed=passed))
                i==1 && continue
                previous=traced[1][i-1]
                if candidate.distinct && previous.distinct &&
                    previous.delta_omega_high_minus_low*candidate.delta_omega_high_minus_low<0
                    a=fixed_temperature_boundary(point.T_MeV,point.muB_MeV,points[i-1].xi,point.xi,base;seed_pair=candidate)
                    b=fixed_temperature_boundary(point.T_MeV,point.muB_MeV,points[i-1].xi,point.xi,refined;seed_pair=check)
                    push!(boundaries,(scan_group=group,T_MeV=point.T_MeV,muB_MeV=point.muB_MeV,
                        xi=a.xi,refined_xi=b.xi,refinement_delta_xi=b.xi-a.xi,
                        bracket_lower_xi=points[i-1].xi,bracket_upper_xi=point.xi))
                    write_rows(joinpath(outdir,"phase_boundaries.csv"),boundaries)
                end
            end
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(scan_group=group,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
        end
        write_rows(joinpath(outdir,"background.csv"),rows)
        write_json(joinpath(outdir,"progress.json"),(completed=length(rows),total=length(plan.points),failures=length(failures)))
        println("BACKGROUND ",length(rows),"/",length(plan.points)," failures=",length(failures));flush(stdout)
    end
    passed=isempty(failures) && length(rows)==909 && all(r->r.passed,rows) && length(audits)>=90 && all(r->r.passed,audits)
    finish_run(outdir;mode="B",points=length(plan.points),completed_points=length(rows),failures=length(failures),
        phase_plan_passed=passed,first_order_boundaries=length(boundaries),cold_seed_audits=length(audits),
        cold_seed_audits_passed=!isempty(audits) && all(r->r.passed,audits),transport_qualified=false)
    passed || error("mode-B background grid failed; inspect retained evidence before preparing transport batches")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: run_transport_mode_b_background_validation.jl INPUT_DIR NEW_OUTPUT_DIR")
    mode_b_background_validation(abspath(ARGS[1]),abspath(ARGS[2]))
end
