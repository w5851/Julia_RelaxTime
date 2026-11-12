include(joinpath(@__DIR__, "run_endpoint_comparison.jl"))

module PALCPathReuseComparison

using Main.Models
using Main.PALCConstraintComparison: row_record, finite_json
using Main.PALCEndpointComparison: all_root_records, anchor_comparison_cases
using LinearAlgebra
using JSON

const HBARC = Main.Constants_PNJL.ħc_MeV_fm
const METHODS = (:baseline, :workspace, :cubic, :residual_first, :ad, :newton, :combined)

function run_method(model, path, method; previous=nothing)
    method = method === :anchor ? :baseline : method === :optimized ? :combined : method
    method in METHODS || throw(ArgumentError("unknown path reuse method"))
    workspace = method in (:workspace, :combined)
    predictor = method in (:cubic, :residual_first, :combined) ? :cubic : :linear
    residual_first = method in (:residual_first, :combined)
    jacobian = method in (:ad, :combined) ? :forwarddiff : :finite_difference
    work = Models.PALCWorkTelemetry()
    strategy = Models.PALCContinuation(; jacobian, target_predictor=predictor)
    controls = (; iterations=300, reuse_corrector=workspace, corrector_residual_first=residual_first)
    controls = (; controls..., corrector_method=method in (:newton,:combined) ? :newton : :trust_region)
    density = path isa Union{Models.FixedRhoPath,Models.FixedAsymmetricRhoPath}
    samples = parse(Int,get(ENV,"PALC_TIMING_SAMPLES","1"))
    1 <= samples <= 20 || throw(ArgumentError("PALC_TIMING_SAMPLES must be 1..20"))
    timed = @timed begin
        output = nothing
        for _ in 1:samples
            work = Models.PALCWorkTelemetry()
            result = Models.solve_path(model, path;
                continuation_strategy=density ? Models.NaturalContinuation() : strategy,
                anchor_strategy=previous === nothing ? nothing : Models._reusable_endpoint_anchor(previous, path),
                (density ? (; recovery_strategy=strategy) : (;))...,
                path_work=work, controls...)
            selected = Models.selected_path_results(model, result; path_work=work, controls...)
            output = (; result, selected)
        end
        output
    end
    return (; timed.value..., seconds=timed.time/samples, allocated_bytes=timed.bytes/samples,
        gc_seconds=timed.gctime/samples, work=Models.palc_work_snapshot(work))
end

function write_records(destination, records)
    mkpath(dirname(destination))
    data = (; schema="palc_path_reuse_v1", source_sha=get(ENV,"PALC_SOURCE_SHA",
        get(ENV,"GITHUB_SHA","working_copy")), actions_run_id=get(ENV,"GITHUB_RUN_ID",""),
        julia_version=string(VERSION), threads=Threads.nthreads(),
        timing_samples=parse(Int,get(ENV,"PALC_TIMING_SAMPLES","1")),
        timing_contract="Path solve and selected-result reconstruction, same process and BLAS thread; validation and serialization excluded. First round includes compilation; subsequent rounds rotate methods.",
        state_tolerance=1e-4, thermo_tolerance=1e-6, residual_tolerance=1e-6, records)
    open(destination,"w") do io
        JSON.print(io,finite_json(data))
        println(io)
    end
end

function record_run(model, name, path, method, repetition; previous=nothing)
    println("START ",name,"/",method," round ",repetition)
    flush(stdout)
    run = run_method(model,path,method; previous)
    rows = [row_record(model,path,value,point) for (value,point) in zip(Models._path_values(path),run.selected)]
    roots = all_root_records(model,path,run.result)
    record = (; name, method, repetition, seconds=run.seconds,
        allocated_bytes=run.allocated_bytes,gc_seconds=run.gc_seconds,
        config=Models.to_namedtuple(path), rows, all_roots=roots, work=run.work,
        diagnostics=Models.to_namedtuple(run.result.diagnostics))
    println("DONE ",name,"/",method," ",run.seconds," s, ",run.allocated_bytes,
        " bytes, iterations ",run.work.root_iterations,", roots ",sum(length(r.roots) for r in roots))
    flush(stdout)
    return record, run.result
end

function temperature_cases()
    sequences = NamedTuple[]
    for route in (:tmu,:rho,:asym), xi in (0.0,0.4), descending in (false,true)
        temperatures = [100.0,110.0,120.0,130.0,140.0,160.0,180.0,200.0]
        descending && reverse!(temperatures)
        paths = [route === :tmu ? Models.FixedMuPath(T/HBARC,collect(280.0:2.0:360.0)./HBARC;xi,p_num=24,t_num=8) :
            route === :rho ? Models.FixedRhoPath(T/HBARC,vcat(0.0,collect(0.2:0.15:3.5));xi,p_num=24,t_num=8) :
            Models.FixedAsymmetricRhoPath(T/HBARC,vcat(0.0,collect(0.2:0.15:3.5));xi,p_num=24,t_num=8) for T in temperatures]
        push!(sequences,(;name="seq_$(route)_$(descending ? "down" : "up")_xi$(Int(10xi))",temperatures,paths))
    end
    return sequences
end

function run_comparison(destination; repetitions=3, cases_filter="all")
    1 <= repetitions <= 5 || throw(ArgumentError("repetitions must be 1..5"))
    BLAS.set_num_threads(1)
    model = Models.create_model(:PNJL)
    cases = anchor_comparison_cases()
    sequences = temperature_cases()
    if cases_filter != "all"
        requested = Symbol.(split(cases_filter,','))
        all(name -> any(c -> first(c) === name,cases) || any(s -> Symbol(s.name) === name,sequences),requested) ||
            throw(ArgumentError("unknown path case"))
        cases = filter(c -> first(c) in requested,cases)
        sequences = filter(s -> Symbol(s.name) in requested,sequences)
    end
    records = NamedTuple[]
    for sequence in sequences, repetition in 1:repetitions
        for method in circshift([:baseline,:anchor,:optimized],repetition-1)
            previous = nothing
            for (T,path) in zip(sequence.temperatures,sequence.paths)
                record, result = record_run(model,"$(sequence.name)_t$(Int(T))",path,method,repetition;
                    previous=method === :baseline ? nothing : previous)
                push!(records,record)
                previous = result
            end
            write_records(destination,records)
        end
    end
    for (name,path) in cases, repetition in 1:repetitions
        methods = path isa Union{Models.FixedRhoPath,Models.FixedAsymmetricRhoPath} ?
            (:baseline,:workspace,:newton,:combined) : METHODS
        for method in circshift(collect(methods), repetition-1)
            record, _ = record_run(model,name,path,method,repetition)
            push!(records,record)
            write_records(destination,records)
        end
    end
    return records
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PALCPathReuseComparison.run_comparison(ARGS[1];
        repetitions=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 3,
        cases_filter=length(ARGS)>=3 ? ARGS[3] : "all")
end
