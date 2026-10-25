#!/usr/bin/env julia
using CSV, JSON, SHA
include(joinpath(@__DIR__, "run_threshold_anomaly_validation.jl"))

refinement_query(cache::C,state::S) where {C,S} = s->evaluate_cache(cache,state,s)

function verify_retained_kernel(run)
    manifest=JSON.parsefile(joinpath(run,"manifest.json"))
    for item in manifest["source_files"]
        relative=replace(item["path"],'\\'=>'/')
        # Sampler/driver changes are intentional. The underlying numerical
        # kernel and model/physics parameters must match the donor run.
        (startswith(relative,"src/") || startswith(relative,"config/models/") ||
            startswith(relative,"config/physics/")) || continue
        path=joinpath(ROOT,relative)
        bytes2hex(sha256(read(path)))==item["sha256"] || error("Retained kernel differs: $relative")
    end
end

function refinement_main(root,outdir,run_names; rtol=1e-5)
    inputdir=joinpath(root,"anomaly_selection_v1")
    records=Dict{Tuple{String,String},Any}()
    for name in run_names
        run=joinpath(root,name)
        verify_retained_kernel(run)
        for file in ("review_points.csv","channel_plan.csv")
            read(joinpath(run,"inputs",file))==read(joinpath(inputdir,file)) ||
                error("Retained background input differs: $name/$file")
        end
        file=joinpath(run,"rates.csv")
        isfile(file) || continue
        for row in CSV.File(file)
            records[(String(row.point_id),String(row.channel))]=(run=run,row=row)
        end
    end
    pending=sort!([key for (key,value) in records if !(value.row.original_replayed &&
        value.row.reference_converged && value.row.adaptive_validated)])
    isempty(pending) && error("No unresolved cases to refine")
    all(records[key].row.original_replayed for key in pending) ||
        error("Original replay mismatches require separate investigation before sampling refinement")
    start_run(outdir,inputdir,@__FILE__,(purpose="refine failed cache/reference targets while preserving previous evidence",
        donors=run_names,cases=pending,cache_rtol=rtol,p_nodes=128,angle_nodes=20,phi_nodes=36,
        original_replay_target=1e-8,reference_rate_target=2e-4,adaptive_rate_target=3e-4,
        background_solve=false,production_write=false),
        extra_sources=[joinpath(@__DIR__,f) for f in ("run_threshold_anomaly_validation.jl",
            "ThresholdCacheCandidates.jl","ThermalAverageCandidates.jl","GroupedThermalProjection.jl")])
    # Freeze all reused direct samples as inputs of this new experiment.
    manifest=JSON.parsefile(joinpath(outdir,"manifest.json"))
    for key in pending
        donor=records[key]
        table_name=join(key,"_")*".csv"
        source=joinpath(donor.run,"tables",table_name)
        saved=joinpath(outdir,"inputs","retained_tables",table_name)
        mkpath(dirname(saved));cp(source,saved)
        push!(manifest["input_files"],Dict("path"=>abspath(source),
            "sha256"=>bytes2hex(sha256(read(saved))),"snapshot"=>relpath(saved,outdir)))
    end
    write_json(joinpath(outdir,"manifest.json"),manifest)
    points=Dict(String(r.point_id)=>r for r in CSV.File(joinpath(inputdir,"review_points.csv")))
    results=NamedTuple[];failures=NamedTuple[]
    last_point="";state=nothing
    for key in pending
        point_id,channel=key;process=Symbol(channel)
        donor=records[key];old=donor.row;row=points[point_id]
        id=join(key,"_")
        println("REFINE ",id);flush(stdout)
        try
            if point_id!=last_point
                Main.RelaxTime.PolarizationCache.reset_cache!()
                state=retained_state(row);last_point=point_id
            end
            old_fine_level=hasproperty(old,:reference_fine_level) ? Int(old.reference_fine_level) : 384
            levels=old.reference_converged ? (old_fine_level,2old_fine_level) : (2old_fine_level,4old_fine_level)
            samples=Dict(Float64(r.s)=>Float64(r.sigma) for r in
                CSV.File(joinpath(outdir,"inputs","retained_tables",id*".csv")))
            built_time=@timed build_anomaly_tables(process,state;rtol=rtol,known_samples=samples,reference_levels=levels)
            built=built_time.value
            primary=refinement_query(built.stock,state)
            knots=sort!(unique(vcat(built.s0 .+ built.local_table.nodes,built.s0 .+ built.tail_table.nodes,
                built.references[1].nodes,built.references[2].nodes)))
            plan=prepare_thermal_plan(process,state)
            bounds=ThermalAverageCandidates.process_bounds(process,state)
            timed=@timed grouped_project_linear_rates(plan,knots,built.s0,primary,
                (built.candidate,built.references[1],built.references[2]),(built.s0+built.xmax,0.0,0.0),bounds)
            replay,adaptive,coarse,fine=timed.value
            source=Float64(old.source_rate)
            original_error=replay/source-1;reference_error=coarse/fine-1;adaptive_error=adaptive/fine-1
            previous_fine=hasproperty(old,:reference_fine_rate) ? Float64(old.reference_fine_rate) : Float64(old.dense384_rate)
            reference_replay_checked=levels[1]==old_fine_level
            reference_replay_error=reference_replay_checked ? coarse/previous_fine-1 : NaN
            reference_replay_checked && abs(reference_replay_error)>1e-10 &&
                error("Same fixed reference changed during refinement: $reference_replay_error")
            push!(results,(point_id=point_id,channel=channel,window_ids=String(row.window_ids),
                muB_MeV=Float64(row.muB_MeV),alpha_T=Float64(row.alpha_T),xi=Float64(row.xi),
                source_rate=source,replayed_rate=replay,adaptive_rate=adaptive,
                reference_coarse_rate=coarse,reference_fine_rate=fine,
                reference_coarse_level=levels[1],reference_fine_level=levels[2],cache_rtol=rtol,
                reference_replay_checked=reference_replay_checked,reference_replay_relative_error=reference_replay_error,
                original_relative_error=original_error,reference_relative_change=reference_error,
                adaptive_relative_error=adaptive_error,superseded_run=basename(donor.run),
                superseded_reference_relative_change=Float64(old.reference_relative_change),
                superseded_adaptive_relative_error=Float64(old.adaptive_relative_error),
                old_reference_converged=Bool(old.reference_converged),old_adaptive_validated=Bool(old.adaptive_validated),
                reused_direct_samples=length(samples),new_direct_evaluations=built.new_calls,
                local_nodes=length(built.local_table.nodes),tail_nodes=length(built.tail_table.nodes),
                build_seconds=built_time.time,thermal_seconds=timed.time,
                original_replayed=abs(original_error)<=1e-8,reference_converged=abs(reference_error)<=2e-4,
                adaptive_validated=abs(adaptive_error)<=3e-4))
            write_rows(joinpath(outdir,"rates.csv"),results)
            write_rows(joinpath(outdir,"tables",id*".csv"),
                [(s=s,sigma=value) for (s,value) in sort!(collect(built.values);by=first)])
            println("RESULT levels=",levels," replay=",original_error," ref=",reference_error,
                " adaptive=",adaptive_error," new_calls=",built.new_calls,
                " seconds=",built_time.time+timed.time);flush(stdout)
        catch err
            message=sprint(showerror,err,catch_backtrace())
            push!(failures,(point_id=point_id,channel=channel,error=message))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",id," ",message);flush(stdout)
        end
    end
    finish_run(outdir;attempted=length(pending),successful=length(results),failures=length(failures),
        original_replayed=count(r->r.original_replayed,results),reference_converged=count(r->r.reference_converged,results),
        adaptive_validated=count(r->r.adaptive_validated,results))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)>=3 || error("Usage: refine_threshold_anomaly_validation.jl RUN_ROOT NEW_OUTPUT_DIR DONOR_RUNS...")
    refinement_main(abspath(ARGS[1]),abspath(ARGS[2]),ARGS[3:end])
end
