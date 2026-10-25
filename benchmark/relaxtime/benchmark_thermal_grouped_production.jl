#!/usr/bin/env julia
include(joinpath(@__DIR__,"benchmark_grouped_thermal_projection.jl"))

function production_group_benchmark(inputdir,outdir)
    analysis=joinpath(ROOT,"scripts","analysis","relaxtime")
    start_run(outdir,inputdir,@__FILE__,(purpose="public production API angular grouping, all 17 channels",
        cases=[(muB=0.0,xi=0.37),(muB=450.0,xi=0.26),(muB=900.0,xi=-0.04)],
        p_nodes=128,angle_nodes=20,phi_nodes=36,t_nodes=24,
        summation_parity_rtol=1e-10,production_write=false),
        extra_sources=[joinpath(@__DIR__,"benchmark_grouped_thermal_projection.jl"),
            [joinpath(analysis,f) for f in ("run_threshold_anomaly_validation.jl",
            "ThresholdCacheCandidates.jl","ThermalAverageCandidates.jl","GroupedThermalProjection.jl")]...])
    points=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    results=NamedTuple[]
    for (mu,xi) in ((0.0,0.37),(450.0,0.26),(900.0,-0.04))
        row=only(r for r in points if r.muB_MeV==mu && r.alpha_T==1.0 && r.xi==xi)
        state=retained_state(row)
        for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
            Main.RelaxTime.PolarizationCache.reset_cache!()
            cache_t=@timed AS.build_w0cdf_pchip_cache(process,state.quark_params,
                state.thermo_params,state.K_coeffs;N=560,p_cutoff=Main.Constants_PNJL.Λ_inv_fm,
                n_sigma_points=24,threshold_subtraction=true,sigma_cache_policy=:validated_anchored)
            cache=cache_t.value
            options=merge(rate_config(),(cs_cache=cache,))
            for method in (:none,:grouped)
                AS.average_scattering_rate(process,state.quark_params,state.thermo_params,state.K_coeffs;
                    merge(options,(p_nodes=6,angle_nodes=3,phi_nodes=5))...,angular_reduction=method)
            end
            repeat=(process===:uubar_to_uubar || (mu==0.0 && process===:us_to_us)) ? 3 : 1
            original=NaN
            for trial in 1:repeat, method in (isodd(trial) ? (:none,:grouped) : (:grouped,:none))
                GC.gc(); cpu0=grouped_process_cpu_seconds()
                measured=@timed AS.average_scattering_rate(process,state.quark_params,state.thermo_params,
                    state.K_coeffs;options...,angular_reduction=method)
                cpu=grouped_process_cpu_seconds()-cpu0
                isnan(original) && (original=measured.value)
                difference=abs(measured.value-original)/max(abs(original),1e-30)
                push!(results,(muB_MeV=mu,xi=xi,channel=String(process),trial=trial,
                    method=String(method),rate=measured.value,relative_difference=difference,
                    wall_seconds=measured.time,cpu_seconds=cpu,allocated_bytes=measured.bytes,
                    cache_wall_seconds=cache_t.time))
                write_rows(joinpath(outdir,"timings.csv"),results)
                difference<=1e-10 || error("angular grouping discrepancy exceeds numerical roundoff allowance")
            end
            println("DONE ",mu," ",xi," ",process);flush(stdout)
        end
    end
    finish_run(outdir;cases=51,measurements=length(results),numerical_parity=true,
        max_relative_difference=maximum(r.relative_difference for r in results))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: benchmark_thermal_grouped_production.jl INPUT_DIR NEW_OUTPUT_DIR")
    production_group_benchmark(abspath(ARGS[1]),abspath(ARGS[2]))
end
