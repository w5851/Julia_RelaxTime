#!/usr/bin/env julia
using Profile, CSV, JSON, Dates
const THRESHOLD_CACHE_BENCH_ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(THRESHOLD_CACHE_BENCH_ROOT,"scripts","analysis","relaxtime","run_threshold_cache_candidates.jl"))

function build_strategy(case,method;dense_level=192)
    state=case.state; process=case.process
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    # Every timed construction starts from the same production grid and fresh
    # polarization cache. JIT compilation and fixed background solving are excluded.
    stock=AS.build_w0cdf_pchip_cache(process,qp,tp,kc;N=560,
        p_cutoff=Main.Constants_PNJL.Λ_inv_fm,n_sigma_points=24,
        threshold_subtraction=true,asym_window=0.6,asym_fit_min_points=8,
        asym_extra_points=10,sigma_cache_policy=:validated_anchored)
    raw=Dict(s=>AS._exact_cached_sigma_value(stock,s,stock.sigma_vals[i])
             for (i,s) in enumerate(stock.s_vals))
    new_calls=Ref(0)
    function direct(s)
        new_calls[]+=1
        return CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
    end
    for x in (case.xmin,case.xmax)
        s=case.s0+x
        haskey(raw,s) || (raw[s]=direct(s))
    end
    if method===:dense
        for level in (96,192,384,768)
            level>dense_level && break
            xs=vcat(10.0 .^ range(log10(case.s0*1e-10),log10(case.xmax);length=level),
                    [case.xmax*(i/level)^2 for i in 1:level])
            for x in xs
                s=case.s0+x
                haskey(raw,s) || (raw[s]=direct(s))
            end
        end
        grid=sort!(collect(keys(raw)))
        cache=AS.CrossSectionCache(process,grid,[raw[s] for s in grid])
        query=s->evaluate_cache(cache,state,s)
        return (query=query,nodes=length(grid),new_evaluations=new_calls[],
            fallback_calls=Ref(0),J=local_integral(cache,state;order=12))
    end
    method in (:adaptive_raw,:adaptive_x_sigma) || error("Unknown strategy $method")
    estimates=threshold_estimates(state)
    local_case=merge(case,(stock=stock,))
    seeds=physical_seeds(local_case,estimates)
    grid=sort!(collect(keys(raw)))
    tail_cache=AS.CrossSectionCache(process,grid,[raw[s] for s in grid])
    tail(s)=evaluate_cache(tail_cache,state,s)
    representation=method===:adaptive_raw ? :raw : :x_sigma
    table=build_adaptive_sigma(x->direct(case.s0+x),seeds;
        representation=representation,rtol=1e-4,
        known_values=Dict(s-case.s0=>v for (s,v) in raw))
    table.converged || error("Adaptive benchmark table failed to converge")
    fallback=Ref(0)
    function fallback_direct(s)
        fallback[]+=1
        return CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
    end
    query=composite_sigma(table,case.s0,tail,fallback_direct)
    return (query=query,nodes=length(table.nodes)+count(s->s-case.s0>case.xmax,grid),
        new_evaluations=new_calls[],fallback_calls=fallback,J=table_weighted_integral(table))
end

function query_checksum(query,points;repeats=20)
    total=0.0
    for _ in 1:repeats, s in points
        total+=query(s)
    end
    return total
end

function process_cpu_seconds()
    if Sys.iswindows()
        handle=ccall((:GetCurrentProcess,"kernel32"),Ptr{Cvoid},())
        created,exited,kernel,user=(Ref{UInt64}(0) for _ in 1:4)
        ok=ccall((:GetProcessTimes,"kernel32"),Int32,
            (Ptr{Cvoid},Ref{UInt64},Ref{UInt64},Ref{UInt64},Ref{UInt64}),
            handle,created,exited,kernel,user)
        ok!=0 || error("GetProcessTimes failed")
        return Float64(kernel[]+user[])*1e-7
    end
    # POSIX specifies CLOCKS_PER_SEC=1e6. Keep CPU and wall time separate.
    ticks=ccall(:clock,Clong,())
    ticks>=0 || error("Process CPU clock is unavailable")
    return ticks/1e6
end

function retained_reference(refdir,outdir)
    previous=JSON.parsefile(joinpath(refdir,"manifest.json"))
    current=JSON.parsefile(joinpath(outdir,"manifest.json"))
    JSON.parsefile(joinpath(refdir,"status.json"))["complete"]===true ||
        error("Reference run is incomplete")
    current_sources=Dict(r["path"]=>r["sha256"] for r in current["source_files"])
    # This benchmark may change, but the physics/candidate implementation may
    # not change when reusing a completed dense384/dense768 reference.
    for row in previous["source_files"]
        row_path=row["path"]
        row_path==relpath(@__FILE__,ROOT) && continue
        get(current_sources,row_path,nothing)==row["sha256"] ||
            error("Reference source differs: $row_path")
    end
    sort([r["sha256"] for r in previous["input_files"]])==
        sort([r["sha256"] for r in current["input_files"]]) ||
        error("Reference inputs differ")
    return JSON.parsefile(joinpath(refdir,"reference.json"))
end

function main(ablationdir,outdir;refdir=nothing)
    start_run(outdir,joinpath(ablationdir,"inputs"),@__FILE__,(
        purpose="matched-error cache construction, lookup and full thermal-rate benchmark",
        xi=0.37,channel="uubar_to_uubar",target_relative_error=1e-4,
        rate_reference_levels=[384,768],measured_repetitions=3,
        trial_order="balanced cyclic order; all strategies warmed before measurement",
        process_cpu_clock=Sys.iswindows() ? "GetProcessTimes user+kernel" : "POSIX clock",
        retained_reference=refdir,
        background="retained fixed state; gap solve excluded",
        cache_reset="clear polarization cache before each timed construction",
        query_probe="fixed half threshold/log and half tail samples; full-rate timing is separate",
        cpu=Sys.CPU_NAME,threads=Threads.nthreads()),
        extra_sources=[joinpath(THRESHOLD_CACHE_BENCH_ROOT,"scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"),
                       joinpath(THRESHOLD_CACHE_BENCH_ROOT,"scripts","analysis","relaxtime","run_threshold_cache_candidates.jl")])
    case=load_upgrade_case(ablationdir,0.37,:uubar_to_uubar)
    reference_rates=Float64[]
    if refdir===nothing
        for level in (384,768)
            println("REFERENCE dense",level);flush(stdout)
            built=build_strategy(case,:dense;dense_level=level)
            push!(reference_rates,functional_rate(built.query,case.process,case.state))
        end
    else
        kept=retained_reference(refdir,outdir)
        kept["J"]==case.reference || error("Retained local reference differs")
        append!(reference_rates,[kept["rate384"],kept["rate768"]])
        println("REFERENCE reused ",refdir);flush(stdout)
    end
    reference=last(reference_rates)
    reference_change=reference_rates[1]/reference-1
    abs(reference_change)<2e-5 || error("Reference rate has not converged to one fifth of the target")
    write_json(joinpath(outdir,"reference.json"),(
        J=case.reference,rate384=reference_rates[1],rate768=reference,
        relative_change=reference_change,reused_from=refdir,
        reused_sha256=refdir===nothing ? nothing :
            bytes2hex(sha256(read(joinpath(refdir,"reference.json"))))))
    _,_,sup=AS._resolve_sigma_support_bounds(case.state.quark_params.m.u,
        case.state.quark_params.m.u,case.state.quark_params.m.u,case.state.quark_params.m.u,
        Main.Constants_PNJL.Λ_inv_fm,nothing,nothing)
    points=vcat(case.s0 .+ 10.0 .^ range(log10(case.xmin*2),log10(case.xmax);length=2000),
        [case.s0+case.xmax+(sup-case.s0-case.xmax)*(i-0.5)/2000 for i in 1:2000])
    rows=NamedTuple[]
    methods=(:dense,:adaptive_raw,:adaptive_x_sigma)
    process_cpu_seconds() # Compile/probe the clock before any measurement.
    for method in methods
        println("WARM ",method);flush(stdout)
        warm=build_strategy(case,method)
        functional_rate(warm.query,case.process,case.state;p_nodes=4,angle_nodes=2,phi_nodes=2)
        query_checksum(warm.query,points;repeats=1)
    end
    for trial in 1:3
        for offset in 0:2
            method=methods[mod1(trial+offset,3)]
            println("TIMED ",method," ",trial);flush(stdout)
            Main.RelaxTime.PolarizationCache.reset_cache!()
            GC.gc()
            started=string(now(UTC))
            cpu_start=process_cpu_seconds()
            construction=@timed build_strategy(case,method)
            build_cpu=process_cpu_seconds()-cpu_start
            built=construction.value
            calls_before=built.fallback_calls[]
            cpu_start=process_cpu_seconds()
            rate=@timed functional_rate(built.query,case.process,case.state)
            rate_cpu=process_cpu_seconds()-cpu_start
            rate_fallbacks=built.fallback_calls[]-calls_before
            cpu_start=process_cpu_seconds()
            lookup=@timed query_checksum(built.query,points)
            query_cpu=process_cpu_seconds()-cpu_start
            rate_error=rate.value/reference-1
            J_error=built.J/case.reference-1
            push!(rows,(method=String(method),trial=trial,measurement_order=length(rows)+1,
                started_utc=started,nodes=built.nodes,
                new_sigma_evaluations=built.new_evaluations,
                build_seconds=construction.time,build_bytes=construction.bytes,
                rate_seconds=rate.time,rate_bytes=rate.bytes,
                total_seconds=construction.time+rate.time,
                build_cpu_seconds=build_cpu,rate_cpu_seconds=rate_cpu,
                total_cpu_seconds=build_cpu+rate_cpu,query_cpu_seconds=query_cpu,
                query_seconds=lookup.time,query_count=length(points)*20,
                query_checksum=lookup.value,fallback_calls=rate_fallbacks,
                rate=rate.value,rate_relative_error=rate_error,J_relative_error=J_error,
                meets_target=abs(rate_error)<=1e-4 && abs(J_error)<=1e-4))
            write_rows(joinpath(outdir,"timings.csv"),rows)
        end
    end
    # Profiling is separate from all timed trials.
    profile_case=build_strategy(case,:adaptive_x_sigma)
    Profile.clear()
    @profile functional_rate(profile_case.query,case.process,case.state)
    open(joinpath(outdir,"profile_flat.txt"),"w") do io
        Profile.print(io;format=:flat,sortedby=:count,mincount=20)
    end
    finish_run(outdir;measurements=length(rows),all_meet_target=all(r.meets_target for r in rows),
        reference_relative_change=reference_change)
    println("COMPLETE ",outdir)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (2,3) || error("Usage: benchmark_threshold_cache.jl ABLATION_DIR NEW_OUTPUT_DIR [REFERENCE_RUN_DIR]")
    main(abspath(ARGS[1]),abspath(ARGS[2]);refdir=length(ARGS)==3 ? abspath(ARGS[3]) : nothing)
end
