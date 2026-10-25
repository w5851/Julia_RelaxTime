#!/usr/bin/env julia
include(joinpath(@__DIR__, "benchmark_threshold_cache.jl"))
include(joinpath(THRESHOLD_CACHE_BENCH_ROOT, "scripts", "analysis", "relaxtime", "ThermalAverageCandidates.jl"))
using .ThermalAverageCandidates

function main_precomputation(ablationdir, outdir)
    start_run(outdir, joinpath(ablationdir, "inputs"), @__FILE__, (
        purpose="same-node thermal quadrature optimization; original cache held fixed",
        cases=["muB0 alpha1 xi0.37 uubar", "muB0 alpha1 xi-0.20 uubar", "muB0 alpha1 xi0 us"],
        repetitions=3, p_nodes=128, angle_nodes=20, phi_nodes=36,
        accuracy_target=1e-11, threads=Threads.nthreads(), cpu=Sys.CPU_NAME),
        extra_sources=[joinpath(@__DIR__, "benchmark_threshold_cache.jl"),
            joinpath(ROOT, "scripts", "analysis", "relaxtime", "ThresholdCacheCandidates.jl"),
            joinpath(ROOT, "scripts", "analysis", "relaxtime", "ThermalAverageCandidates.jl"),
            joinpath(ROOT, "scripts", "analysis", "relaxtime", "run_threshold_cache_candidates.jl")])
    measurements = NamedTuple[]
    for (xi, process) in ((0.37, :uubar_to_uubar), (-0.2, :uubar_to_uubar), (0.0, :us_to_us))
        state = retained_state(retained_row(joinpath(ablationdir, "inputs"), xi))
        qp, tp, kc = state.quark_params, state.thermo_params, state.K_coeffs
        stock = AS.build_w0cdf_pchip_cache(process, qp, tp, kc; N=560,
            p_cutoff=Main.Constants_PNJL.Λ_inv_fm, n_sigma_points=24,
            threshold_subtraction=true, sigma_cache_policy=:validated_anchored)
        query = s -> evaluate_cache(stock, state, s)
        reference = functional_rate(query, process, state)
        for velocity in (:legacy, :equal_mass)
            prepared_rate(query, process, state; velocity=velocity, p_nodes=4, angle_nodes=2, phi_nodes=2)
        end
        shared_rates((query, query), (process, process), state; p_nodes=4, angle_nodes=2, phi_nodes=2)
        println("CASE ", xi, " ", process, " reference=", reference); flush(stdout)
        methods = (:original, :precomputed, :equal_mass, :shared_pair)
        for trial in 1:3, offset in 0:3
            method = methods[mod1(trial+offset, 4)]
            GC.gc()
            cpu0 = process_cpu_seconds()
            timed = if method === :original
                @timed functional_rate(query, process, state)
            elseif method === :shared_pair
                @timed shared_rates((query, query), (process, process), state)
            else
                @timed prepared_rate(query, process, state;
                    velocity=method === :precomputed ? :legacy : :equal_mass)
            end
            cpu = process_cpu_seconds() - cpu0
            values = method === :shared_pair ? timed.value : (timed.value,)
            err = maximum(abs(v/reference-1) for v in values)
            err <= 1e-11 || error("Numerical parity failed: $method $err")
            push!(measurements, (xi=xi, process=String(process), trial=trial,
                method=String(method), wall_seconds=timed.time, cpu_seconds=cpu,
                allocated_bytes=timed.bytes, process_count=length(values),
                reference=reference, value=first(values), relative_error=err))
            write_rows(joinpath(outdir, "timings.csv"), measurements)
            println(method, " wall=", timed.time, " cpu=", cpu, " error=", err); flush(stdout)
        end
    end
    finish_run(outdir; measurements=length(measurements), numerical_parity=true)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("Usage: benchmark_thermal_precomputation.jl ABLATION_DIR NEW_OUTPUT_DIR")
    main_precomputation(abspath(ARGS[1]), abspath(ARGS[2]))
end
