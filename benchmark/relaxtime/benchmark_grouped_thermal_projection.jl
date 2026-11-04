#!/usr/bin/env julia
include(joinpath(@__DIR__, "..", "..", "scripts", "analysis", "relaxtime", "run_threshold_anomaly_validation.jl"))
if !isdefined(Main, :GroupedThermalProjection)
    include(joinpath(@__DIR__, "..", "..", "scripts", "analysis", "relaxtime", "GroupedThermalProjection.jl"))
end
using .GroupedThermalProjection

function grouped_process_cpu_seconds()
    if Sys.iswindows()
        handle=ccall((:GetCurrentProcess,"kernel32"),Ptr{Cvoid},())
        created,exited,kernel,user=(Ref{UInt64}(0) for _ in 1:4)
        ok=ccall((:GetProcessTimes,"kernel32"),Int32,
            (Ptr{Cvoid},Ref{UInt64},Ref{UInt64},Ref{UInt64},Ref{UInt64}),
            handle,created,exited,kernel,user)
        ok!=0 || error("GetProcessTimes failed")
        return Float64(kernel[]+user[])*1e-7
    end
    ticks=ccall(:clock,Clong,())
    ticks>=0 || error("Process CPU clock is unavailable")
    return ticks/1e6
end

function grouped_benchmark(inputdir,outdir)
    analysis=joinpath(ROOT,"scripts","analysis","relaxtime")
    start_run(outdir,inputdir,@__FILE__,(purpose="exact grouping of repeated angular geometry versus full union-grid projection",
        repetitions=3,p_nodes=128,angle_nodes=20,phi_nodes=36,rtol=1e-11,
        primary_reference="full original nodes with Kahan sum",historical_replay_target=1e-8,production_write=false),
        extra_sources=[joinpath(analysis,f) for f in ("ThresholdCacheCandidates.jl","ThermalAverageCandidates.jl",
            "GroupedThermalProjection.jl","run_threshold_anomaly_validation.jl")])
    points=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    source_plan=collect(CSV.File(joinpath(inputdir,"channel_plan.csv")))
    measurements=NamedTuple[]
    for (mu,xi,process) in ((0.0,0.37,:uubar_to_uubar),(450.0,0.26,:uubar_to_uubar),(0.0,0.37,:us_to_us))
        Main.RelaxTime.PolarizationCache.reset_cache!()
        row=only(r for r in points if r.muB_MeV==mu && r.alpha_T==1.0 && r.xi==xi)
        println("BUILD ",mu," ",xi," ",process);flush(stdout)
        state=retained_state(row)
        built=build_anomaly_tables(process,state)
        primary=let cache=built.stock, fixed_state=state
            s->evaluate_cache(cache,fixed_state,s)
        end
        queries=(built.candidate,built.references[1],built.references[2])
        weighted=(built.s0+built.xmax,0.0,0.0)
        knots=sort!(unique(vcat(built.s0 .+ built.local_table.nodes,built.s0 .+ built.tail_table.nodes,
            built.references[1].nodes,built.references[2].nodes)))
        bounds=ThermalAverageCandidates.process_bounds(process,state)
        small=prepare_thermal_plan(process,state;p_nodes=6,angle_nodes=4,phi_nodes=8)
        project_linear_rates(small,knots,built.s0,primary,queries,weighted,bounds)
        grouped_project_linear_rates(small,knots,built.s0,primary,queries,weighted,bounds)
        plan=prepare_thermal_plan(process,state)
        groups=angular_groups(plan)
        full_events=length(plan.cw)^2*length(plan.phiw)
        grouped_events=sum(length(g.cosines) for g in groups)
        println("ANGULAR EVENTS ",full_events," -> ",grouped_events);flush(stdout)
        source=only(unique(Float64(r.rate) for r in source_plan if r.point_id==row.point_id && Symbol(r.channel)==process))
        compensated_primary=compensated_primary_rate(plan,primary,bounds)
        println("COMPENSATED PRIMARY ",compensated_primary," source_bias=",source/compensated_primary-1);flush(stdout)
        reference=nothing
        for trial in 1:3, method in (isodd(trial) ? (:full,:grouped) : (:grouped,:full))
            GC.gc(); cpu0=grouped_process_cpu_seconds()
            t=method===:full ? @timed(project_linear_rates(plan,knots,built.s0,primary,queries,weighted,bounds)) :
                @timed(grouped_project_linear_rates(plan,knots,built.s0,primary,queries,weighted,bounds;groups=groups))
            cpu=grouped_process_cpu_seconds()-cpu0
            reference===nothing && (reference=t.value)
            legacy_difference=maximum(abs(a/b-1) for (a,b) in zip(t.value,reference))
            projection_difference=maximum(abs(a/b-1) for (a,b) in zip(t.value[2:end],reference[2:end]))
            compensated_difference=abs(t.value[1]/compensated_primary-1)
            difference=method===:grouped ? max(projection_difference,compensated_difference) : projection_difference
            replay=abs(t.value[1]/source-1)
            replay<=1e-8 || error("Historical replay failed: $replay")
            push!(measurements,(muB_MeV=mu,xi=xi,channel=String(process),trial=trial,method=String(method),
                wall_seconds=t.time,cpu_seconds=cpu,bytes=t.bytes,relative_error=difference,
                legacy_relative_error=legacy_difference,compensated_primary_error=compensated_difference,
                legacy_primary_bias=source/compensated_primary-1,
                original_relative_error=replay,old_rate=t.value[1],adaptive_rate=t.value[2],
                coarse_rate=t.value[3],fine_rate=t.value[4],full_angular_events=full_events,grouped_angular_events=grouped_events))
            write_rows(joinpath(outdir,"timings.csv"),measurements)
            difference<=1e-11 || error("Grouped projection/compensated primary parity failed: $difference")
            println(mu," ",process," ",method," wall=",t.time," cpu=",cpu," error=",difference);flush(stdout)
        end
    end
    finish_run(outdir;measurements=length(measurements),numerical_parity=true)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: benchmark_grouped_thermal_projection.jl INPUT_DIR NEW_OUTPUT_DIR")
    grouped_benchmark(abspath(ARGS[1]),abspath(ARGS[2]))
end
