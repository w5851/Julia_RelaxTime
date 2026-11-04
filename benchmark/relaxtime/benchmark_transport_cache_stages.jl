#!/usr/bin/env julia
# Attribute native cache construction cost before choosing production orders.
include(joinpath(@__DIR__,"..","..","scripts","analysis","relaxtime","run_transport_joint_validation.jl"))

function cache_stage_measurement(process,state,config)
    ts=AS.ThermalWeightedSampling
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    started=time_ns()
    prepared=@timed AS._prepare_thermal_cache_problem(process,qp,tp,kc;
        p_nodes=config.p,angle_nodes=config.c,phi_nodes=config.phi,n_sigma_points=config.t)
    problem=prepared.value;sampling=problem.sampling;context=problem.context
    measured=@timed ts.thermal_measure(problem.plan,context.bounds;max_events=sampling.max_events)
    measure=measured.value;s0=first(context.bounds)
    featured=@timed ts.scattering_feature_seeds(process,state,measure,s0)
    direct_nanoseconds=Ref(UInt64(0));calls=Ref(0)
    direct=s->begin
        before=time_ns()
        value=CS.total_cross_section(process,s,context.quark_params,context.thermo_params,
            context.K_coeffs;n_points=context.n_sigma_points)
        direct_nanoseconds[]+=time_ns()-before;calls[]+=1
        value
    end
    built=@timed ts.build_thermal_weighted_table(direct,measure,s0;
        rtol=sampling.rtol,atol=sampling.atol,weighted_window=sampling.weighted_window,
        seed_order=sampling.seed_order,thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2],
        extra_seeds=featured.value.seeds,max_evaluations=sampling.max_evaluations,
        max_nodes=sampling.max_nodes,max_rounds=sampling.max_rounds)
    table=ts.validated_table(built.value)
    return (rate=built.value.rate,estimated_relative_error=built.value.estimated_error/built.value.rate,
        events=measure.events,unique_events=length(measure.nodes),nodes=length(table.nodes),
        evaluations=calls[],prepare_seconds=prepared.time,measure_seconds=measured.time,
        feature_seconds=featured.time,adaptive_seconds=built.time,
        direct_sigma_seconds=Float64(direct_nanoseconds[])*1e-9,
        total_seconds=Float64(time_ns()-started)*1e-9,
        allocated_bytes=prepared.bytes+measured.bytes+featured.bytes+built.bytes)
end

function cache_stage_benchmark(inputdir,outdir)
    rows=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    row=only(r for r in rows if (r.muB_MeV,r.alpha_T,r.xi)==(0.0,1.0,-0.22))
    candidate=joint_config("candidate")
    channels=(:uu_to_uu,:uubar_to_uubar,:us_to_us)
    analysis=joinpath(ROOT,"scripts","analysis","relaxtime")
    start_run(outdir,inputdir,@__FILE__,(purpose="native construction cost attribution and radial order pilot",
        point_id=String(row.point_id),channels=String.(channels),orders=[48,64,96],
        reference_config=candidate,cache_target=1e-3,
        timing="serial stage timings; direct-sigma timer is included in adaptive time",
        acceptance="radial pilot only; not a complete-point production benchmark or final convergence qualification"),
        extra_sources=[joinpath(analysis,f) for f in ("run_transport_joint_validation.jl",
            "run_transport_background_tail_validation.jl","run_transport_background_validation.jl",
            "TerminalTransportQuadrature.jl","ThermalWeightedSampling.jl","ThermalAverageCandidates.jl",
            "ThresholdCacheCandidates.jl","StreamedThermalValidation.jl")])
    background=joint_background(row,candidate);state=background.state
    warmup=@elapsed cache_stage_measurement(first(channels),state,merge(candidate,(p=4,c=3,phi=5,t=4,)))
    write_json(joinpath(outdir,"warmup.json"),(seconds=warmup,includes_JIT=true))
    output=NamedTuple[]
    for process in channels,p in (48,64,96)
        GC.gc();Main.RelaxTime.PolarizationCache.reset_cache!()
        result=cache_stage_measurement(process,state,merge(candidate,(p=p,)))
        push!(output,merge((point_id=String(row.point_id),channel=String(process),p_nodes=p,
            angle_nodes=candidate.c,phi_nodes=candidate.phi),result))
        write_rows(joinpath(outdir,"timings.csv"),output)
        println("STAGES ",process," p=",p," rate=",result.rate," seconds=",result.total_seconds,
            " measure=",result.measure_seconds," adaptive=",result.adaptive_seconds,
            " sigma=",result.direct_sigma_seconds);flush(stdout)
    end
    finish_run(outdir;measurements=length(output),warmup_seconds=warmup)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: benchmark_transport_cache_stages.jl INPUT_DIR NEW_OUTPUT_DIR")
    cache_stage_benchmark(abspath(ARGS[1]),abspath(ARGS[2]))
end
