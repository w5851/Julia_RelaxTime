#!/usr/bin/env julia
# Include the original inaccurate cache as a cost baseline, not as an
# accuracy-qualified alternative. Prior matched-error runs remain unchanged.
include(joinpath(@__DIR__,"benchmark_threshold_cache.jl"))

function build_baseline_strategy(case,method)
    if method!==:stock
        return merge(build_strategy(case,method),(native_cache=nothing,))
    end
    state=case.state
    cache=AS.build_w0cdf_pchip_cache(case.process,state.quark_params,state.thermo_params,state.K_coeffs;
        N=560,p_cutoff=Main.Constants_PNJL.Λ_inv_fm,n_sigma_points=24,
        threshold_subtraction=true,asym_window=0.6,asym_fit_min_points=8,
        asym_extra_points=10,sigma_cache_policy=:validated_anchored)
    query=s->evaluate_cache(cache,state,s)
    return (query=query,nodes=length(cache.s_vals),new_evaluations=0,
        fallback_calls=Ref(0),J=local_integral(cache,state;order=12),native_cache=cache)
end

function original_consumer(case,cache;p_nodes=128,angle_nodes=20,phi_nodes=36)
    state=case.state
    return AS.average_scattering_rate(case.process,state.quark_params,state.thermo_params,state.K_coeffs;
        p_nodes=p_nodes,angle_nodes=angle_nodes,phi_nodes=phi_nodes,n_sigma_points=24,
        cs_cache=cache,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,interpolation_mode=:linear,
        threshold_subtraction=true,sigma_cache_policy=:validated_anchored,
        asym_window=0.6,asym_fit_min_points=8,asym_extra_points=10)
end

function run_baseline_cost(ablationdir,outdir,refdir)
    start_run(outdir,joinpath(ablationdir,"inputs"),@__FILE__,(
        purpose="original inaccurate production cache versus two accurate candidates: runtime decomposition",
        xi=0.37,channel="uubar_to_uubar",target_relative_error=1e-4,
        methods=["stock","dense","adaptive_x_sigma"],measured_repetitions=3,
        trial_order="balanced cyclic order; cache reset before each construction",
        original_consumer="timed before callable parity measurement on stock trials",
        reference_run=refdir,cpu=Sys.CPU_NAME,threads=Threads.nthreads(),
        excluded_costs="Julia startup/JIT, gap solve, other channels, transport integration and plot rendering"),
        extra_sources=[joinpath(@__DIR__,"benchmark_threshold_cache.jl"),
            joinpath(THRESHOLD_CACHE_BENCH_ROOT,"scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"),
            joinpath(THRESHOLD_CACHE_BENCH_ROOT,"scripts","analysis","relaxtime","run_threshold_cache_candidates.jl")])
    case=load_upgrade_case(ablationdir,0.37,:uubar_to_uubar)
    kept=retained_reference(refdir,outdir)
    kept["J"]==case.reference || error("Local integral reference differs")
    reference=kept["rate768"]
    write_json(joinpath(outdir,"reference.json"),(
        J=case.reference,rate=reference,original_rate=case.source_rate,
        retained_file=joinpath(refdir,"reference.json"),
        retained_sha256=bytes2hex(sha256(read(joinpath(refdir,"reference.json"))))))
    methods=(:stock,:dense,:adaptive_x_sigma)
    process_cpu_seconds()
    for method in methods
        println("WARM ",method);flush(stdout)
        built=build_baseline_strategy(case,method)
        functional_rate(built.query,case.process,case.state;p_nodes=4,angle_nodes=2,phi_nodes=2)
        if method===:stock
            original_consumer(case,built.native_cache;p_nodes=4,angle_nodes=2,phi_nodes=2)
        end
    end
    result=NamedTuple[]
    for trial in 1:3, offset in 0:2
        method=methods[mod1(trial+offset,3)]
        println("TIMED ",method," ",trial);flush(stdout)
        Main.RelaxTime.PolarizationCache.reset_cache!()
        GC.gc()
        start_cpu=process_cpu_seconds()
        construction=@timed build_baseline_strategy(case,method)
        build_cpu=process_cpu_seconds()-start_cpu
        built=construction.value
        native_seconds=missing; native_cpu=missing; native_value=missing
        if method===:stock
            start_cpu=process_cpu_seconds()
            native=@timed original_consumer(case,built.native_cache)
            native_seconds=native.time;native_cpu=process_cpu_seconds()-start_cpu;native_value=native.value
            isapprox(native.value,case.source_rate;rtol=1e-10,atol=0) || error("Original consumer no longer reproduces retained rate")
        end
        start_cpu=process_cpu_seconds()
        rate=@timed functional_rate(built.query,case.process,case.state)
        rate_cpu=process_cpu_seconds()-start_cpu
        if method===:stock
            isapprox(rate.value,native_value;rtol=1e-10,atol=0) || error("Callable consumer parity failed")
        end
        j_error=built.J/case.reference-1
        rate_error=rate.value/reference-1
        meets=abs(j_error)<=1e-4 && abs(rate_error)<=1e-4
        method===:stock || meets || error("Accurate candidate failed its target")
        push!(result,(method=String(method),trial=trial,nodes=built.nodes,
            new_sigma_evaluations=built.new_evaluations,
            build_seconds=construction.time,build_cpu_seconds=build_cpu,
            callable_rate_seconds=rate.time,callable_rate_cpu_seconds=rate_cpu,
            original_rate_seconds=native_seconds,original_rate_cpu_seconds=native_cpu,
            total_seconds=construction.time+(method===:stock ? native_seconds : rate.time),
            total_cpu_seconds=build_cpu+(method===:stock ? native_cpu : rate_cpu),
            callable_total_seconds=construction.time+rate.time,
            callable_total_cpu_seconds=build_cpu+rate_cpu,
            build_bytes=construction.bytes,rate=rate.value,J_relative_error=j_error,
            rate_relative_error=rate_error,meets_target=meets))
        write_rows(joinpath(outdir,"timings.csv"),result)
    end
    finish_run(outdir;measurements=length(result),original_retained_rate_reproduced=true,
        accurate_candidates_meet_target=all(r.meets_target for r in result if r.method!="stock"),
        stock_meets_target=all(r.meets_target for r in result if r.method=="stock"))
    println("COMPLETE ",outdir)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==3 || error("Usage: benchmark_threshold_baseline.jl ABLATION_DIR NEW_OUTPUT_DIR REFERENCE_RUN_DIR")
    run_baseline_cost(abspath(ARGS[1]),abspath(ARGS[2]),abspath(ARGS[3]))
end
