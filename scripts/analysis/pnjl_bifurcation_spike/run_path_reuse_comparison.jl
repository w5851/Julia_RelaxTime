include(joinpath(@__DIR__, "run_endpoint_comparison.jl"))

module PALCPathReuseComparison

using Main.Models
using Main.PALCConstraintComparison: row_record, finite_json
using Main.PALCEndpointComparison: all_root_records, anchor_comparison_cases
using LinearAlgebra
using JSON

const HBARC = Main.Constants_PNJL.ħc_MeV_fm
const METHODS = (:baseline, :workspace, :cubic, :ad, :newton, :combined, :default)

function residual_microbenchmark(model, path, result; samples=500)
    root = first(p for branch in result.branches for p in branch.points if p.converged)
    options = Models._path_options(model,path)
    residual = Models._path_residual(model,path,options)
    x, value = root.solution, root.param_value
    out = similar(x)
    fd = Models._palc_jacobian(residual,x,value,Models.PALCWorkTelemetry(),:finite_difference)
    ad = Models._palc_jacobian(residual,x,value,Models.PALCWorkTelemetry(),:forwarddiff)
    probes = ((:residual_allocating,()->residual(x,value)),
        (:residual_inplace,()->Models._path_residual!(out,residual,x,value)),
        (:jacobian_finite_difference,()->fd(x,(;value))),
        (:jacobian_forwarddiff,()->ad(x,(;value))))
    measurements = map(probes) do (name,f)
        f()
        run = @timed for _ in 1:samples
            f()
        end
        (;name,seconds=run.time/samples,allocated_bytes=run.bytes/samples)
    end
    return (;samples,coordinate=value,measurements,
        jacobian_max_abs_difference=maximum(abs.(fd(x,(;value))-ad(x,(;value)))))
end

function run_method(model, path, method; previous=nothing)
    if method === :default && !(path isa Models.FixedRhoPath)
        previous = nothing
    end
    method = method === :anchor ? :baseline : method === :optimized ? :combined : method
    method in METHODS || throw(ArgumentError("unknown path reuse method"))
    workspace = method in (:workspace, :combined)
    predictor = method in (:cubic, :combined) ? :cubic : :linear
    jacobian = method in (:ad, :combined) ? :forwarddiff : :finite_difference
    work = Models.PALCWorkTelemetry()
    strategy = method === :default ? Models.PALCContinuation() :
        Models.PALCContinuation(; jacobian, target_predictor=predictor)
    controls = method === :default ? (;iterations=300) :
        (;iterations=300,reuse_corrector=workspace,
            corrector_method=method in (:newton,:combined) ? :newton : :trust_region)
    density = path isa Union{Models.FixedRhoPath,Models.FixedAsymmetricRhoPath}
    samples = parse(Int,get(ENV,"PALC_TIMING_SAMPLES","1"))
    1 <= samples <= 20 || throw(ArgumentError("PALC_TIMING_SAMPLES must be 1..20"))
    # Discard garbage left by validation/serialization before every timed batch.
    # Collections caused by the timed work itself remain in the reported time.
    GC.gc()
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

function write_records(destination, records, plan)
    mkpath(dirname(destination))
    data = (; schema="palc_path_reuse_v1", source_sha=get(ENV,"PALC_SOURCE_SHA",
        get(ENV,"GITHUB_SHA","working_copy")), actions_run_id=get(ENV,"GITHUB_RUN_ID",""),
        julia_version=string(VERSION), threads=Threads.nthreads(),
        timing_samples=parse(Int,get(ENV,"PALC_TIMING_SAMPLES","1")),
        timing_contract="Path solve and selected-result reconstruction, same process and BLAS thread; validation, serialization and pre-batch GC excluded; in-batch GC included. First round includes compilation; subsequent rounds rotate methods.",
        state_tolerance=1e-4, thermo_tolerance=1e-6, residual_tolerance=1e-6,
        expected_runs=plan, records)
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
    microbenchmark = repetition == 1 && method === :baseline &&
        name in (:tmu_t120_fine,:sigma) ? residual_microbenchmark(model,path,run.result) : nothing
    record = (; name, method, repetition, seconds=run.seconds,
        allocated_bytes=run.allocated_bytes,gc_seconds=run.gc_seconds,
        config=Models.to_namedtuple(path), rows, all_roots=roots, work=run.work,
        diagnostics=Models.to_namedtuple(run.result.diagnostics), microbenchmark)
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

function run_comparison(destination; repetitions=3, cases_filter="all", methods_filter="all")
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
    requested_methods = methods_filter == "all" ? nothing : Tuple(Symbol.(split(methods_filter,',')))
    if requested_methods !== nothing
        :baseline in requested_methods && all(m -> m in METHODS,requested_methods) ||
            throw(ArgumentError("methods must include baseline and supported comparison methods"))
        length(unique(requested_methods)) == length(requested_methods) || throw(ArgumentError("duplicate methods"))
    end
    case_methods(path) = requested_methods !== nothing ? requested_methods :
        path isa Models.ConstraintTPath{Models.FixedEntropy} ? (:baseline,:default) :
        path isa Union{Models.FixedRhoPath,Models.FixedAsymmetricRhoPath} ?
            (:baseline,:workspace,:newton,:default) : METHODS
    sequence_methods = requested_methods === nothing ? (:baseline,:anchor,:optimized,:default) : requested_methods
    plan = NamedTuple[(;name="$(sequence.name)_t$(Int(T))",methods=sequence_methods,repetitions)
        for sequence in sequences for T in sequence.temperatures]
    append!(plan,[(;name=string(name),methods=case_methods(path),repetitions) for (name,path) in cases])
    records = NamedTuple[]
    write_records(destination,records,plan)
    for sequence in sequences, repetition in 1:repetitions
        for method in circshift(collect(sequence_methods),repetition-1)
            previous = nothing
            for (T,path) in zip(sequence.temperatures,sequence.paths)
                record, result = record_run(model,"$(sequence.name)_t$(Int(T))",path,method,repetition;
                    previous=method === :baseline ? nothing : previous)
                push!(records,record)
                previous = result
            end
            write_records(destination,records,plan)
        end
    end
    for (name,path) in cases, repetition in 1:repetitions
        for method in circshift(collect(case_methods(path)), repetition-1)
            record, _ = record_run(model,name,path,method,repetition)
            push!(records,record)
        end
        write_records(destination,records,plan)
    end
    return records
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PALCPathReuseComparison.run_comparison(ARGS[1];
        repetitions=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 3,
        cases_filter=length(ARGS)>=3 ? ARGS[3] : "all",
        methods_filter=length(ARGS)>=4 ? ARGS[4] : "all")
end
