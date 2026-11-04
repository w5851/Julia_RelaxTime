#!/usr/bin/env julia
# Current default CLI: one complete warmup, then repeated complete point runs.
using CSV, JSON, SHA
include(joinpath(@__DIR__,"..","..","scripts","relaxtime","run_phase_guided_transport_scan.jl"))
include(joinpath(@__DIR__,"..","..","scripts","analysis","relaxtime","ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon: start_run, finish_run, write_rows

function entrypoint_cpu_seconds()
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

function benchmark_production_entrypoint(inputdir,outdir,baseline;repeats=2)
    repeats>=2 || throw(ArgumentError("repeats must be at least two"))
    metadata_path=splitext(baseline)[1]*".json"
    metadata=JSON.parsefile(metadata_path)
    bytes2hex(sha256(read(baseline)))==metadata["baseline_sha256"] || error("baseline identity changed")
    metadata["baseline_semantics"]=="production_entrypoint" || error("a production entrypoint baseline is required")
    metadata["regression_rtol"]==1e-6 && metadata["regression_atol"]==1e-8 || error("unexpected baseline tolerance")
    targets=collect(CSV.File(baseline))
    input=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    points=[only(r for r in input if (r.muB_MeV,r.alpha_T,r.xi)==key)
        for key in ((0.0,1.0,-0.22),(900.0,1.0,-0.04))]
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    profile.source_sha256==metadata["profile_sha256"] || error("baseline sampling profile changed")
    sources=filter(p->endswith(p,".jl"),readdir(joinpath(PROJECT_ROOT,"scripts","relaxtime");join=true))
    append!(sources,[baseline,metadata_path])
    start_run(outdir,inputdir,@__FILE__,(purpose="default production entrypoint repeatability and complete-point host timing",
        points=String.(getproperty.(points,:point_id)),repeats=repeats,
        profile=Models.TransportSamplingProfiles.effective_config(profile),baseline=abspath(baseline),
        baseline_sha256=metadata["baseline_sha256"],regression_rtol=1e-6,regression_atol=1e-8,
        independent_relative_target=0.01,threads=Threads.nthreads(),cpu=Sys.CPU_NAME,
        warmup="first complete default-entrypoint run; retained and excluded from steady timings",
        timing="includes stable background, all 17 channel builds, bulk and terminal transport; host-specific",
        allocations="allocated_bytes is cumulative allocation, not peak resident memory",
        scope="two representative points; not a runner timing or a full-scan cost bound"),extra_sources=sources)
    jobs=vcat([(first(points),0,true)],[(row,repeat,false) for row in points for repeat in 1:repeats])
    timings=NamedTuple[];comparisons=NamedTuple[];failures=NamedTuple[]
    for (row,repeat,warmup) in jobs
        point=String(row.point_id);role=warmup ? "warmup" : "measurement"
        println("START ",role," ",point," repeat=",repeat);flush(stdout)
        try
            folder=joinpath(outdir,"points",point*"_"*role*"_"*string(repeat))
            # Omit --sampling-profile deliberately to exercise the actual default.
            args=["--outdir",folder,"--mode","fixed-T-sparse-muB","--channel-diagnostics",
                "--T-list",string(row.T_MeV),"--muB-list",string(row.muB_MeV),"--xi-list",string(row.xi)]
            opts=PhaseGuidedCLI.parse_args(args)
            reference=_load_runtime_phase_reference(opts)
            context=ProvenanceMetadata.new_run_context("production_entrypoint_timing",args)
            GC.gc();cpu0=entrypoint_cpu_seconds()
            measured=@timed run_phase_guided_scan(opts,context;phase_reference=reference)
            cpu=entrypoint_cpu_seconds()-cpu0
            path=joinpath(folder,"phase_guided_transport_scan.csv")
            actual=only(CSV.File(path;comment="#"))
            actual.converged && (actual.T_MeV,actual.muB_MeV,actual.xi)==(row.T_MeV,row.muB_MeV,row.xi) ||
                error("default entrypoint point is missing or unconverged")
            sampling=filter(r->r.run_id==actual.run_id,collect(CSV.File(
                TransportSamplingScanConfig.sampling_diagnostics_path((output=path,)))))
            required=Set(String.(RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
            length(sampling)==17 && Set(String.(getproperty.(sampling,:channel)))==required || error("incomplete channels")
            all(r->r.converged && r.estimated_error<=profile.sampling.atol+profile.sampling.rtol*abs(r.rate),sampling) ||
                error("channel budget exceeded")
            expected=filter(r->r.point_id==point,targets)
            length(expected)==9 || error("incomplete baseline point")
            for target in expected
                value=Float64(getproperty(actual,Symbol(target.observable)))
                passed=isfinite(value) && value>0 && isapprox(value,target.candidate;rtol=1e-6,atol=1e-8) &&
                    abs(value/target.reference-1)<=0.01
                push!(comparisons,(point_id=point,repeat=repeat,role=role,observable=String(target.observable),
                    actual=value,baseline=target.candidate,reference=target.reference,
                    baseline_relative_change=value/target.candidate-1,reference_relative_change=value/target.reference-1,passed=passed))
            end
            push!(timings,(point_id=point,repeat=repeat,role=role,total_seconds=measured.time,cpu_seconds=cpu,
                gc_seconds=measured.gctime,allocated_bytes=measured.bytes,qualified_channels=length(sampling)))
            write_rows(joinpath(outdir,"timings.csv"),timings)
            write_rows(joinpath(outdir,"comparisons.csv"),comparisons)
            println("DONE ",role," ",point," seconds=",measured.time);flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,repeat=repeat,role=role,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
        end
    end
    passed=isempty(failures) && length(timings)==length(jobs) && length(comparisons)==9length(jobs) && all(r->r.passed,comparisons)
    finish_run(outdir;points=length(points),measurements=count(r->r.role=="measurement",timings),
        warmups=count(r->r.role=="warmup",timings),comparisons=length(comparisons),failures=length(failures),
        default_entrypoint_passed=passed)
    return passed
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    3<=length(ARGS)<=4 || error("Usage: benchmark_transport_production_entrypoint.jl INPUT_DIR NEW_OUTPUT_DIR BASELINE_CSV [REPEATS]")
    benchmark_production_entrypoint(abspath(ARGS[1]),abspath(ARGS[2]),abspath(ARGS[3]);
        repeats=length(ARGS)==4 ? parse(Int,ARGS[4]) : 2) || exit(1)
end
