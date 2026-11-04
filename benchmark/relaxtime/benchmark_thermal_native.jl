#!/usr/bin/env julia
include(joinpath(@__DIR__, "benchmark_threshold_baseline.jl"))

function native_benchmark(inputdir,outdir)
    start_run(outdir,inputdir,@__FILE__,(purpose="actual optimized production consumer versus original independent quadrature",
        cases=["muB0 alpha1 xi0.37 uubar","muB450 alpha1 xi0.26 uubar","muB0 alpha1 xi0.37 us"],
        repetitions=3,p_nodes=128,angle_nodes=20,phi_nodes=36,rtol=1e-11),
        extra_sources=[joinpath(@__DIR__,"benchmark_threshold_baseline.jl"),
            joinpath(@__DIR__,"benchmark_threshold_cache.jl"),
            joinpath(ROOT,"scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"),
            joinpath(ROOT,"scripts","analysis","relaxtime","run_threshold_cache_candidates.jl")])
    rows=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    measurements=NamedTuple[]
    for (mu,xi,process) in ((0.0,0.37,:uubar_to_uubar),(450.0,0.26,:uubar_to_uubar),(0.0,0.37,:us_to_us))
        row=only(r for r in rows if r.muB_MeV==mu && r.alpha_T==1.0 && r.xi==xi)
        state=retained_state(row);qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
        cache=AS.build_w0cdf_pchip_cache(process,qp,tp,kc;N=560,p_cutoff=Main.Constants_PNJL.Λ_inv_fm,
            n_sigma_points=24,threshold_subtraction=true,sigma_cache_policy=:validated_anchored)
        query=s->evaluate_cache(cache,state,s)
        case=(state=state,process=process)
        original_consumer(case,cache;p_nodes=4,angle_nodes=2,phi_nodes=2)
        functional_rate(query,process,state;p_nodes=4,angle_nodes=2,phi_nodes=2)
        reference=functional_rate(query,process,state)
        for trial in 1:3, method in (isodd(trial) ? (:old,:native) : (:native,:old))
            GC.gc();cpu0=process_cpu_seconds()
            t=method===:old ? @timed(functional_rate(query,process,state)) : @timed(original_consumer(case,cache))
            cpu=process_cpu_seconds()-cpu0
            error=abs(t.value/reference-1)
            error<=1e-11 || Base.error("native numerical parity failed")
            push!(measurements,(muB_MeV=mu,xi=xi,channel=String(process),trial=trial,method=String(method),
                wall_seconds=t.time,cpu_seconds=cpu,bytes=t.bytes,value=t.value,reference=reference,relative_error=error))
            write_rows(joinpath(outdir,"timings.csv"),measurements)
            println(mu," ",xi," ",process," ",method," wall=",t.time," cpu=",cpu," error=",error);flush(stdout)
        end
    end
    finish_run(outdir;measurements=length(measurements),numerical_parity=true)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: benchmark_thermal_native.jl INPUT_DIR NEW_OUTPUT_DIR")
    native_benchmark(abspath(ARGS[1]),abspath(ARGS[2]))
end
