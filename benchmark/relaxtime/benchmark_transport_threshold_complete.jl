#!/usr/bin/env julia
include(joinpath(@__DIR__,"..","..","scripts","analysis","relaxtime","run_transport_joint_validation.jl"))

function benchmark_cpu_seconds()
    if Sys.iswindows()
        handle=ccall((:GetCurrentProcess,"kernel32"),Ptr{Cvoid},())
        created,exited,kernel,user=(Ref{UInt64}(0) for _ in 1:4)
        ok=ccall((:GetProcessTimes,"kernel32"),Int32,
            (Ptr{Cvoid},Ref{UInt64},Ref{UInt64},Ref{UInt64},Ref{UInt64}),handle,created,exited,kernel,user)
        ok!=0 || error("GetProcessTimes failed")
        return Float64(kernel[]+user[])*1e-7
    end
    return ccall(:clock,Clong,())/1e6
end

function legacy_benchmark_config()
    return (name="legacy",bg_limit=10.0,bg_p=24,bg_c=8,a_p=16,a_c=4,legacy_A=true,
        p=128,c=20,phi=36,t=24,tr_p=24,tr_c=8,tr_limit=8.0)
end

function benchmark_native_point(row,config)
    total0=time_ns();cpu0=benchmark_cpu_seconds()
    bg=@timed joint_background(row,config)
    state=bg.value.state;qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    caches=Dict{Symbol,AS.ThermalWeightedCache}()
    options=(p_nodes=config.p,angle_nodes=config.c,phi_nodes=config.phi,n_sigma_points=config.t,
        sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm)
    builds=NamedTuple[]
    for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
        GC.gc();Main.RelaxTime.PolarizationCache.reset_cache!()
        measured=@timed AS.build_thermal_weighted_cache(process,qp,tp,kc;options...)
        cache=measured.value;caches[process]=cache
        push!(builds,(channel=String(process),seconds=measured.time,allocated_bytes=measured.bytes,
            evaluations=cache.statistics.evaluations,events=cache.statistics.events,nodes=cache.statistics.nodes))
        println("BUILD thermal_weighted ",process," seconds=",measured.time);flush(stdout)
    end
    kwargs=merge(options,(sigma_cache_policy=:thermal_weighted,cs_caches=caches,propagator_quark_params=qp))
    output=@timed joint_transport(row,bg.value,config;tau_kwargs=kwargs,check_density=false)
    elapsed=(time_ns()-total0)/1e9
    cpu_elapsed=benchmark_cpu_seconds()-cpu0
    # A second complete downstream call demonstrates actual cache reuse.
    reused=@timed joint_transport(row,bg.value,config;tau_kwargs=kwargs,check_density=false)
    output.value.result.tau==reused.value.result.tau || error("cache reuse changed native tau")
    return (background_seconds=bg.time,construction_seconds=sum(r.seconds for r in builds),
        projection_seconds=0.0,transport_seconds=output.time,total_seconds=elapsed,
        cpu_seconds=cpu_elapsed,reuse_seconds=reused.time,result=output.value.result,builds=builds)
end

function benchmark_legacy_point(row,config;grid_n=560)
    total0=time_ns();cpu0=benchmark_cpu_seconds();bg=@timed joint_background(row,config)
    state=bg.value.state;qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    rates=Dict(method=>Dict{Symbol,Float64}() for method in (:none,:grouped))
    times=Dict(method=>0.0 for method in (:none,:grouped));builds=NamedTuple[]
    projection_cpu=Dict(method=>0.0 for method in (:none,:grouped))
    for (i,process) in enumerate(Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU)
        GC.gc();Main.RelaxTime.PolarizationCache.reset_cache!()
        built=@timed AS.build_w0cdf_pchip_cache(process,qp,tp,kc;N=grid_n,
            p_cutoff=Main.Constants_PNJL.Λ_inv_fm,n_sigma_points=config.t,
            threshold_subtraction=true,sigma_cache_policy=:validated_anchored)
        push!(builds,(channel=String(process),seconds=built.time,allocated_bytes=built.bytes,
            evaluations=length(built.value.s_vals),events=0,nodes=length(built.value.s_vals)))
        options=(p_nodes=config.p,angle_nodes=config.c,phi_nodes=config.phi,n_sigma_points=config.t,
            cs_cache=built.value,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
            threshold_subtraction=true,interpolation_mode=:linear,sigma_cache_policy=:validated_anchored)
        for method in (isodd(i) ? (:none,:grouped) : (:grouped,:none))
            cpu_start=benchmark_cpu_seconds()
            measured=@timed AS.average_scattering_rate(process,qp,tp,kc;options...,angular_reduction=method)
            projection_cpu[method]+=benchmark_cpu_seconds()-cpu_start
            rates[method][process]=measured.value;times[method]+=measured.time
        end
        isapprox(rates[:none][process],rates[:grouped][process];rtol=1e-10) ||
            error("legacy angular reduction changed a rate beyond roundoff")
        println("BUILD legacy ",process," seconds=",built.time);flush(stdout)
    end
    result=@timed joint_transport(row,bg.value,config;
        tau_kwargs=(existing_rates=NamedTuple(rates[:none]),),check_density=false)
    construction=sum(r.seconds for r in builds)
    wall=(time_ns()-total0)/1e9
    total_cpu=benchmark_cpu_seconds()-cpu0
    overhead=max(0.0,wall-bg.time-construction-sum(values(times))-result.time)
    return [(method=String(method),background_seconds=bg.time,construction_seconds=construction,
        projection_seconds=times[method],transport_seconds=result.time,
        total_seconds=bg.time+construction+times[method]+result.time+overhead,
        cpu_seconds=total_cpu-sum(values(projection_cpu))+projection_cpu[method],
        reuse_seconds=times[method]+result.time,result=result.value.result,builds=builds) for method in (:none,:grouped)]
end

function complete_benchmark(inputdir,outdir;repeats=2,candidate_p::Union{Nothing,Int}=nothing)
    repeats>0 || error("repeats must be positive")
    points=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    rows=[only(r for r in points if (r.muB_MeV,r.alpha_T,r.xi)==key)
        for key in ((0.0,1.0,-0.22),(900.0,1.0,-0.04))]
    candidate=joint_config("candidate");legacy=legacy_benchmark_config()
    candidate_p===nothing || candidate_p>0 || error("candidate momentum order must be positive")
    candidate_p===nothing || (candidate=merge(candidate,(p=candidate_p,)))
    analysis=joinpath(ROOT,"scripts","analysis","relaxtime")
    start_run(outdir,inputdir,@__FILE__,(purpose="serial complete-point timing including background, cache, all 17 rates, bulk and terminal transport",
        points=[String(r.point_id) for r in rows],repeats=repeats,candidate=candidate,legacy=legacy,
        threads=Threads.nthreads(),cpu=Sys.CPU_NAME,total_memory_bytes=Sys.total_memory(),free_memory_bytes=Sys.free_memory(),
        warmup="both complete paths with small outer grids; excluded from steady timings",
        legacy_timing="shared construction plus independently timed original/grouped projections; reported as component sums",
        memory="serial channels with explicit GC included in total; no parallel numerical jobs",
        accuracy="legacy is a cost baseline, not a qualified numerical reference"),
        extra_sources=[joinpath(analysis,f) for f in ("run_transport_joint_validation.jl",
            "run_transport_background_tail_validation.jl","run_transport_background_validation.jl",
            "TerminalTransportQuadrature.jl","ThermalWeightedSampling.jl","ThermalAverageCandidates.jl",
            "ThresholdCacheCandidates.jl","StreamedThermalValidation.jl")])
    println("WARMUP START");flush(stdout)
    warmup=@elapsed begin
        benchmark_legacy_point(first(rows),merge(legacy,(p=6,c=3,phi=5,t=4,));grid_n=32)
        benchmark_native_point(first(rows),merge(candidate,(p=6,c=3,phi=5,t=4,)))
    end
    write_json(joinpath(outdir,"warmup.json"),(seconds=warmup,includes_JIT=true))
    println("WARMUP COMPLETE seconds=",warmup);flush(stdout)
    timings=NamedTuple[];channels=NamedTuple[];observables=NamedTuple[]
    for row in rows,repeat in 1:repeats
        for kind in (isodd(repeat) ? (:legacy,:candidate) : (:candidate,:legacy))
            GC.gc()
            println("START ",row.point_id," repeat=",repeat," method=",kind);flush(stdout)
            measured=kind===:legacy ? benchmark_legacy_point(row,legacy) :
                [merge((method="thermal_weighted",),benchmark_native_point(row,candidate))]
            for result in measured
                key=(point_id=String(row.point_id),repeat=repeat,method=result.method)
                push!(timings,merge(key,(background_seconds=result.background_seconds,
                    construction_seconds=result.construction_seconds,projection_seconds=result.projection_seconds,
                    transport_seconds=result.transport_seconds,total_seconds=result.total_seconds,cpu_seconds=result.cpu_seconds,
                    reuse_seconds=result.reuse_seconds)))
                append!(channels,[merge(key,r) for r in result.builds])
                push!(observables,merge(key,joint_observables(result.result)))
                println("TIMING ",key," seconds=",result.total_seconds);flush(stdout)
            end
            write_rows(joinpath(outdir,"timings.csv"),timings)
            write_rows(joinpath(outdir,"channels.csv"),channels)
            write_rows(joinpath(outdir,"observables.csv"),observables)
        end
    end
    finish_run(outdir;points=length(rows),measurements=length(timings),warmup_seconds=warmup)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=4 || error("Usage: benchmark_transport_threshold_complete.jl INPUT_DIR NEW_OUTPUT_DIR [REPEATS CANDIDATE_P]")
    complete_benchmark(abspath(ARGS[1]),abspath(ARGS[2]);repeats=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 2,
        candidate_p=length(ARGS)>=4 ? parse(Int,ARGS[4]) : nothing)
end
