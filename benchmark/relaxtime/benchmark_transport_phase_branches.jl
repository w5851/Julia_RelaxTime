#!/usr/bin/env julia
# Compare explicit local phase seeds with the general solver's perturbation pool.
include(joinpath(@__DIR__,"..","..","scripts","analysis","relaxtime","run_transport_phase_revalidation.jl"))

function benchmark_phase_branches(inputdir,outdir)
    phase=GapTransportScanPhaseEquilibrium
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    sources=filter(p->endswith(p,".jl"),readdir(joinpath(PROJECT_ROOT,"scripts","relaxtime");join=true))
    start_run(outdir,inputdir,@__FILE__,(purpose="local branch seeds versus general perturbed-seed solver",
        mass_relative_tolerance=1e-6,omega_absolute_tolerance=1e-8,
        residual_norm_max=1e-8,includes_jit_in_first_case=true),extra_sources=sources)
    rows=NamedTuple[]
    for (cutoff,p,c) in ((25.0,96,24),(30.0,128,32))
        numerics=(model=Models.TransportSamplingProfiles.profile_model(profile;thermal_p_max_inv_fm=cutoff),
            p_num=p,t_num=c,residual_norm_max=1e-8)
        for (T,mu,xi) in ((200.08890140807532,0.0,0.0),(182.76306553605758,450.0,0.0),
            (125.73725802686799,900.0,-0.003),(125.73725802686799,900.0,0.003))
            println("BRANCH START T=",T," muB=",mu," xi=",xi," p=",p);flush(stdout)
            GC.gc()
            local_result=@timed phase.solve_two_branch_candidates(T/ħc_MeV_fm,mu/(3ħc_MeV_fm),xi,numerics)
            println("  local seconds=",local_result.time);flush(stdout)
            GC.gc()
            original=@timed phase.solve_two_branch_candidates(T/ħc_MeV_fm,mu/(3ħc_MeV_fm),xi,numerics;
                branch_solver=phase.solve_models_equilibrium)
            a,b=local_result.value,original.value
            mass_delta=maximum(abs.(a.stable.masses./b.stable.masses.-1))
            omega_delta=abs(a.stable.omega-b.stable.omega)
            branches_match=a.distinct==b.distinct && (!a.distinct || a.stable_phase==b.stable_phase)
            passed=branches_match && mass_delta<=1e-6 && omega_delta<=1e-8
            push!(rows,(T_MeV=T,muB_MeV=mu,xi=xi,p_nodes=p,cos_nodes=c,thermal_p_max_inv_fm=cutoff,
                local_seconds=local_result.time,original_seconds=original.time,
                local_gc_seconds=local_result.gctime,original_gc_seconds=original.gctime,
                local_bytes=local_result.bytes,original_bytes=original.bytes,
                mass_relative_change=mass_delta,omega_absolute_change=omega_delta,
                local_distinct=a.distinct,original_distinct=b.distinct,
                local_phase=String(a.stable_phase),original_phase=String(b.stable_phase),passed=passed))
            write_rows(joinpath(outdir,"comparisons.csv"),rows)
            println("BRANCH DONE original seconds=",original.time," mass_delta=",mass_delta," passed=",passed);flush(stdout)
            passed || error("local phase branch solve differs from the original solver")
        end
    end
    finish_run(outdir;points=length(rows),branch_solver_passed=all(r->r.passed,rows))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=3 || error("Usage: benchmark_transport_phase_branches.jl INPUT_DIR BENCHMARK_OUTPUT [PHASE_OUTPUT]")
    benchmark_phase_branches(abspath(ARGS[1]),abspath(ARGS[2]))
    length(ARGS)==3 && phase_revalidation(abspath(ARGS[1]),abspath(ARGS[3]))
end
