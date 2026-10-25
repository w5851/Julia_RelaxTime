#!/usr/bin/env julia
# Exercise the real scan entrypoint against the frozen, qualified joint results.
using CSV, JSON, SHA
include(joinpath(@__DIR__,"..","..","relaxtime","run_phase_guided_transport_scan.jl"))
include(joinpath(@__DIR__,"ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon: start_run, finish_run, write_json, write_rows

function production_validation(inputdir,outdir;scope="pilot",shard=0,shards=1,point_ids=nothing,qualified_dir=nothing,extra_sources=String[])
    0<=shard<shards || error("invalid shard")
    scope in ("pilot","remaining","all") || error("invalid production validation scope")
    qualified=qualified_dir===nothing ? joinpath(dirname(inputdir),"joint_complete_summary_v1") : abspath(qualified_dir)
    summary=JSON.parsefile(joinpath(qualified,"summary.json"))
    expected_points=point_ids===nothing ? 12 : length(point_ids)
    expected_points>0 && summary["joint_target_passed"] && summary["completed_points"]==expected_points ||
        error("joint qualification incomplete")
    bytes2hex(sha256(read(joinpath(inputdir,"review_points.csv"))))==summary["input_sha256"] ||
        error("qualification input changed")
    comparisons=collect(CSV.File(joinpath(qualified,"comparisons.csv")))
    ids=Set(String.(getproperty.(comparisons,:point_id)))
    if point_ids!==nothing
        scope=="all" || throw(ArgumentError("explicit point_ids requires scope=all"))
        length(ids)==length(point_ids) && ids==Set(String.(point_ids)) ||
            throw(ArgumentError("qualified comparisons must cover every requested point_id exactly"))
    end
    pilot="muB0.0_alpha1.0_xi-0.2200"
    rows=filter(r->String(r.point_id) in ids,collect(CSV.File(joinpath(inputdir,"review_points.csv"))))
    length(rows)==length(ids) || error("qualification points are missing or duplicated")
    scope=="pilot" && filter!(r->r.point_id==pilot,rows)
    scope=="remaining" && filter!(r->r.point_id!=pilot,rows)
    rows=[r for (i,r) in enumerate(rows) if (i-1)%shards==shard]
    isempty(rows) && error("empty validation shard")
    sources=filter(p->endswith(p,".jl"),readdir(joinpath(PROJECT_ROOT,"scripts","relaxtime");join=true))
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    start_run(outdir,inputdir,@__FILE__,(purpose="real phase-guided production scan replay at frozen qualification coordinates",
        scope=scope,shard=shard,shards=shards,points=String.(getproperty.(rows,:point_id)),requested_point_ids=point_ids,
        profile=Models.TransportSamplingProfiles.effective_config(profile),relative_target=0.01,
        replay_relative_tolerance=1e-4,reference_directory=abspath(qualified),
        reference_comparisons_sha256=bytes2hex(sha256(read(joinpath(qualified,"comparisons.csv"))))),
        extra_sources=vcat(sources,extra_sources))
    cp(joinpath(qualified,"comparisons.csv"),joinpath(outdir,"inputs","qualified_comparisons.csv"))
    output_rows=NamedTuple[];failures=NamedTuple[];timings=NamedTuple[]
    for row in rows
        point=String(row.point_id)
        println("START production ",point);flush(stdout)
        try
            folder=joinpath(outdir,"points",point)
            args=["--sampling-profile",profile.source_path,"--outdir",folder,
                "--mode","fixed-T-sparse-muB","--channel-diagnostics",
                "--T-list",string(row.T_MeV),"--muB-list",string(row.muB_MeV),"--xi-list",string(row.xi)]
            opts=PhaseGuidedCLI.parse_args(args)
            reference=_load_runtime_phase_reference(opts)
            ctx=ProvenanceMetadata.new_run_context("production_sampling_validation",args)
            elapsed=@elapsed run_phase_guided_scan(opts,ctx;phase_reference=reference)
            result_path=joinpath(folder,"phase_guided_transport_scan.csv")
            actual=collect(CSV.File(result_path;comment="#"))
            length(actual)==1 || error("production entrypoint did not produce exactly one valid point")
            actual=only(actual)
            sampling=collect(CSV.File(TransportSamplingScanConfig.sampling_diagnostics_path((output=result_path,))))
            required=Set(String.(RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
            length(sampling)==length(required) && Set(String.(getproperty.(sampling,:channel)))==required ||
                error("production channel diagnostics incomplete")
            all(s->s.converged && s.estimated_error<=profile.sampling.atol+profile.sampling.rtol*abs(s.rate),sampling) ||
                error("production channel error budget failed")
            targets=filter(c->c.point_id==point,comparisons)
            for target in targets
                value=Float64(getproperty(actual,Symbol(target.observable)))
                replay=value/target.candidate-1
                refined=value/target.reference-1
                push!(output_rows,(point_id=point,observable=String(target.observable),production=value,
                    qualified_candidate=target.candidate,qualified_reference=target.reference,
                    replay_relative_change=replay,reference_relative_change=refined,
                    passed=isfinite(value) && abs(replay)<=1e-4 && abs(refined)<=0.01))
            end
            push!(timings,(point_id=point,seconds=elapsed,includes_jit=isempty(timings),
                context="validation runtime; not a dedicated benchmark"))
            write_rows(joinpath(outdir,"comparisons.csv"),output_rows)
            write_rows(joinpath(outdir,"timings.csv"),timings)
            println("DONE production ",point," comparisons=",length(targets));flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED production ",point," ",sprint(showerror,err));flush(stdout)
        end
        write_json(joinpath(outdir,"progress.json"),(completed_points=length(timings),failures=length(failures),
            expected_points=length(rows)))
        GC.gc()
    end
    passed=isempty(failures) && length(output_rows)==9length(rows) && all(r->r.passed,output_rows)
    finish_run(outdir;points=length(rows),completed_points=length(timings),comparisons=length(output_rows),
        failures=length(failures),production_replay_passed=passed)
    passed || error("production replay did not pass; inspect comparisons and failures")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=5 || error("Usage: run_transport_production_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|remaining|all SHARD SHARDS]")
    production_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)>=3 ? ARGS[3] : "pilot",
        shard=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 0,shards=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 1)
end
