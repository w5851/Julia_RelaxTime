include(joinpath(@__DIR__, "run_constraint_comparison.jl"))

module PALCEndpointComparison

using Main.Models
using Main.PALCConstraintComparison: row_record, finite_json
using LinearAlgebra
using JSON

const HBARC = Main.Constants_PNJL.ħc_MeV_fm
const METHODS = (:catalog, :single_low, :single_high, :dual_endpoint)

function endpoint_seed(path, value, side)
    mode = Models._path_mode(path, value)
    theta = path isa Models.FixedMuPath ? [path.T_fm, value] : [path.T_fm]
    return Models.get_seed(Models.DefaultSeed(phase_hint=side), theta, mode)
end

"""Natural-parameter sweeps, one root solve per point and no seed-catalog fallback."""
function dual_endpoint_path(model, path, work; sweep::Symbol=:both)
    sweep in (:both, :forward, :reverse) || throw(ArgumentError("unknown natural sweep $(sweep)"))
    started = time()
    options = Models._path_options(model, path; iterations=300)
    residual = Models._path_residual(model, path, options)
    values = sort!(unique(Models._path_values(path)))
    branches = Models.ContinuationBranch[]
    anchors = Models.AnchorRoot[]
    failures = 0
    messages = String[]
    for (direction, grid, side) in ((:forward, values, :hadron), (:reverse, reverse(values), :quark))
        sweep === :both || sweep === direction || continue
        seed = endpoint_seed(path, first(grid), side)
        points = Models.BranchPoint[]
        anchor_id = Symbol("endpoint_", direction)
        for (step, value) in enumerate(grid)
            candidate = Models._path_root(model, path, value, seed, residual, options, work;
                source=step == 1 ? :anchor : :correction)
            if candidate === nothing || !candidate.converged
                failures += 1
                push!(messages, "$(direction) stopped after a failed single-seed solve at $(value)")
                break
            end
            point = Models._candidate_path_point(value, candidate, Symbol("endpoint_", direction), step)
            push!(points, point)
            step == 1 && push!(anchors, Models._anchor_from_branch_point(anchor_id, point, 1))
            seed = copy(candidate.solution)
        end
        push!(branches, Models.ContinuationBranch(direction, anchor_id, :natural_parameter,
            points, length(points) == length(values) ? :complete : :partial, (;)))
    end
    selected = Models.apply_branch_policy(branches, Models.GroundStateBranchPolicy())
    diagnostics = Models.PathDiagnostics(:natural_parameter,
        sweep === :both ? :two_endpoints : :one_endpoint, :ground_state,
        time()-started, work.anchor_solves, 0, failures, 0, messages)
    return Models.PathSolveResult(path, selected.branches, selected.selections, anchors, diagnostics)
end

function run_method(model, path, method)
    work = Models.PALCWorkTelemetry()
    started = time_ns()
    result = if method in (:dual_endpoint, :natural_low, :natural_high)
        sweep = method === :dual_endpoint ? :both : method === :natural_low ? :forward : :reverse
        dual_endpoint_path(model, path, work; sweep)
    elseif method === :hot_start
        hot_start_path(model,path,work)
    elseif method === :production
        Models.solve_path(model,path;continuation_strategy=Models.NaturalContinuation(),
            path_work=work,iterations=300)
    else
        values = Models._path_values(path)
        anchor = if method === :catalog
            Models.SeedPoolAnchor()
        elseif method === :auto
            nothing
        else
            low = method === :single_low
            value = low ? minimum(values) : maximum(values)
            Models.SingleSeedAnchor(value, endpoint_seed(path, value, low ? :hadron : :quark))
        end
        Models.solve_path(model, path; continuation_strategy=Models.PALCContinuation(),
            anchor_strategy=anchor, path_work=work, iterations=300)
    end
    selected = Models.selected_path_results(model, result; iterations=300, path_work=work)
    seconds = (time_ns()-started)/1e9
    return (; result, selected, seconds, work=Models.palc_work_snapshot(work))
end

function hot_start_path(model, path::Models.FixedMuPath, work)
    started = time()
    value = maximum(path.mu_values)
    initial_T = max(250/HBARC,path.T_fm)
    grid = collect(range(initial_T,path.T_fm;length=max(2,ceil(Int,(initial_T-path.T_fm)*HBARC/20)+1)))
    seed = Models.get_seed(Models.DefaultSeed(phase_hint=:quark),[initial_T,value],Models.FixedMu())
    root = nothing
    for T in grid
        point_path = Models.FixedMuPath(T,[value];xi=path.xi,p_num=path.p_num,t_num=path.t_num)
        options = Models._path_options(model,point_path;iterations=300)
        residual = Models._path_residual(model,point_path,options)
        root = Models._path_root(model,point_path,value,seed,residual,options,work;source=:anchor)
        if root === nothing || !root.converged
            points = [Models.BranchPoint(v,Float64[],zeros(5),zeros(3),NaN,NaN,Inf,false,:failed,i)
                for (i,v) in enumerate(Models._path_values(path))]
            branches = [Models.ContinuationBranch(:hot_start,:hot_anchor,:temperature_bridge,
                points,:partial,(;))]
            selected = Models.apply_branch_policy(branches,Models.GroundStateBranchPolicy())
            diagnostics = Models.PathDiagnostics(:bifurcationkit,:hot_start,:ground_state,
                time()-started,work.anchor_solves,0,length(points),0,
                ["High-temperature anchor bridge failed at T=$(T*HBARC) MeV"])
            return Models.PathSolveResult(path,branches,selected.selections,Models.AnchorRoot[],diagnostics)
        end
        seed = root.solution
    end
    return Models._solve_palc_path(model,path,Models.PALCContinuation();
        anchor_strategy=Models.SingleSeedAnchor(value,root.solution),initial_anchor=root,
        path_work=work,iterations=300)
end

function all_root_records(model, path, result)
    options = Models._path_options(model, path; iterations=300)
    residual = Models._path_residual(model, path, options)
    rows = NamedTuple[]
    for value in sort!(unique(Models._path_values(path)))
        roots = NamedTuple[]
        for branch in result.branches, point in branch.points
            point.param_value == value && point.converged || continue
            candidate = Models._path_candidate(model, path, value, point.solution, residual, 0, options)
            Models._push_distinct_path_root!(roots, candidate)
        end
        sort!(roots; by=root -> root.rho_norm)
        push!(rows, (; coordinate=value, roots=[row_record(model,path,value,root) for root in roots]))
    end
    return rows
end

function write_records(destination, records; anchor_probes=NamedTuple[])
    mkpath(dirname(destination))
    data = (; schema="palc_endpoint_comparison_v1", source_sha=strip(read(`git rev-parse HEAD`, String)),
        actions_run_id=get(ENV,"GITHUB_RUN_ID",""), julia_version=string(VERSION),
        p_num=24,t_num=8,iterations=300,residual_tolerance=1e-6,
        timing_contract="Same process, one thread, path solve plus identical selected-result reconstruction; validation and serialization excluded; first round includes compilation.",
        endpoint_contract="Hadron seed at the lower endpoint and quark seed at the upper endpoint; one direct root solve per target, stop a direction on failure, no catalog fallback; Newton jumps are retained as actual roots.",
        records, anchor_probes)
    open(destination,"w") do io
        JSON.print(io,finite_json(data),2)
        println(io)
    end
end

"""Inspect individual start guesses outside the matched path timing."""
function probe_anchors(model, name, path)
    options = Models._path_options(model,path;iterations=300)
    residual = Models._path_residual(model,path,options)
    low, high = extrema(Models._path_values(path))
    middle = (low+high)/2
    attempts = [(;label="endpoint_$(i)",attempt...) for (i,attempt) in
        enumerate(Models._endpoint_anchor_attempts(path))]
    for (label,value) in (("middle",middle),("high",high))
        for (seed_name,base) in (("hadron",Models.SeedStrategies.HADRON_SEED_5),
                               ("medium",Models.SeedStrategies.MEDIUM_SEED_5))
            seed = Models.SeedStrategies.extend_seed(base,Models._path_mode(path,value))
            push!(attempts,(;label="$(label)_$(seed_name)",value,seed))
        end
    end
    probes = NamedTuple[]
    for attempt in attempts
        work = Models.PALCWorkTelemetry()
        started = time_ns()
        root = Models._path_root(model,path,attempt.value,attempt.seed,residual,options,work;
            source=:anchor)
        root_seconds = (time_ns()-started)/1e9
        curves = NamedTuple[]
        if root !== nothing && root.converged
            for (label,strategy) in (("default",Models.PALCContinuation()),
                ("small_step",Models.PALCContinuation(ds=0.001,dsmax=0.01,max_steps=400)))
                curve = nothing
                error_message = ""
                started = time_ns()
                try
                    curve = Models._palc_curve(residual,root,attempt.value,(low,high),strategy,work;
                        positive_parameter=path isa Models.ConstraintTPath)
                catch err
                    err isa InterruptException && rethrow()
                    error_message = sprint(showerror,err)
                end
                push!(curves,(;strategy=label,seconds=(time_ns()-started)/1e9,error=error_message,
                    point_count=curve === nothing ? 0 : length(curve),
                    target_crossings=curve === nothing ? Int[] :
                        [length(Models._palc_target_seeds(curve,value)) for value in Models._path_values(path)]))
            end
        end
        push!(probes,(;name,attempt...,root_seconds,
            root=root === nothing ? nothing : (;solution=root.solution,iterations=root.iterations,
                residual_norm=root.residual_norm,converged=root.converged,
                failed_constraints=root.failed_constraints,rho_norm=root.rho_norm),curves))
    end
    return probes
end

function anchor_comparison_cases()
    cases = Tuple{Symbol,Models.EquilibriumPath}[]
    for T in (100.0,120.0,200.0), xi in (0.0,0.4)
        push!(cases, (Symbol("tmu_t$(Int(T))_xi$(Int(10xi))"),
            Models.FixedMuPath(T/HBARC, collect(0.0:15.0:450.0)./HBARC; xi,p_num=24,t_num=8)))
    end
    push!(cases, (:tmu_t120_fine,
        Models.FixedMuPath(120/HBARC, collect(296.0:0.25:318.0)./HBARC; p_num=24,t_num=8)))
    for T in (100.0,120.0,200.0), xi in (0.0,0.4), asymmetric in (false,true)
        grid = vcat(0.0,collect(0.2:0.15:3.5))
        path = asymmetric ? Models.FixedAsymmetricRhoPath(T/HBARC,grid;xi,p_num=24,t_num=8) :
            Models.FixedRhoPath(T/HBARC,grid;xi,p_num=24,t_num=8)
        push!(cases,(Symbol("$(asymmetric ? "asym" : "rho")_t$(Int(T))_xi$(Int(10xi))"),path))
    end
    push!(cases, (:sigma, Models.ConstraintTPath(Models.FixedSigma(12.0759),
        collect(100.0:2.0:140.0)./HBARC;p_num=24,t_num=8)))
    push!(cases, (:sigma_xi4, Models.ConstraintTPath(Models.FixedSigma(12.0759),
        collect(100.0:2.0:140.0)./HBARC;xi=0.4,p_num=24,t_num=8)))
    push!(cases, (:sigma_t160, Models.ConstraintTPath(Models.FixedSigma(12.0759),
        collect(120.0:2.0:160.0)./HBARC;p_num=24,t_num=8)))
    push!(cases, (:bqs, Models.ConstraintTPath(Models.FixedMuBConservedCharges(900/HBARC),
        collect(100.0:3.0:160.0)./HBARC;p_num=24,t_num=8)))
    push!(cases, (:entropy, Models.ConstraintTPath(Models.FixedEntropy(0.5),
        collect(100.0:1.0:120.0)./HBARC;p_num=24,t_num=8)))
    return cases
end

function run_comparison(destination; repetitions=3, scope="endpoints", cases_filter="all", probe=false)
    1 <= repetitions <= 3 || throw(ArgumentError("repetitions must be 1..3"))
    BLAS.set_num_threads(1)
    model = Models.create_model(:PNJL)
    cases = scope in ("anchors","promotion") ? anchor_comparison_cases() : [
        (:tmu_t120, Models.FixedMuPath(120/HBARC, collect(260.0:2.0:380.0)./HBARC; p_num=24,t_num=8)),
        (:tmu_t120_fine, Models.FixedMuPath(120/HBARC, collect(296.0:0.25:318.0)./HBARC; p_num=24,t_num=8)),
        (:rho_t120, Models.FixedRhoPath(120/HBARC, collect(0.2:0.05:3.5); p_num=24,t_num=8)),
        (:tmu_t200, Models.FixedMuPath(200/HBARC, collect(260.0:2.0:380.0)./HBARC; p_num=24,t_num=8)),
    ]
    if cases_filter != "all"
        requested = Symbol.(split(cases_filter,','))
        all(name -> any(c -> first(c) === name,cases),requested) ||
            throw(ArgumentError("unknown endpoint comparison case"))
        cases = filter(c -> first(c) in requested,cases)
    end
    anchor_probes = NamedTuple[]
    if probe
        for (name,path) in cases
            append!(anchor_probes,probe_anchors(model,name,path))
            write_records(destination,NamedTuple[];anchor_probes)
        end
    end
    records = NamedTuple[]
    for (name,path) in cases, repetition in 1:repetitions
        methods = if scope == "promotion"
            path isa Union{Models.FixedRhoPath,Models.FixedAsymmetricRhoPath} ?
                (:catalog,:auto,:production) : (:catalog,:auto)
        elseif scope == "anchors"
            path isa Union{Models.FixedRhoPath,Models.FixedAsymmetricRhoPath} ?
                (:catalog,:auto,:natural_low,:natural_high) :
                name in (:tmu_t120_xi0,:tmu_t200_xi0) ? (:catalog,:auto,:hot_start) : (:catalog,:auto)
        else
            METHODS
        end
        for method in circshift(collect(methods), repetition-1)
            println("START ",name,"/",method," round ",repetition)
            flush(stdout)
            run = run_method(model,path,method)
            rows = [row_record(model,path,value,point) for (value,point) in zip(Models._path_values(path),run.selected)]
            all_roots = all_root_records(model,path,run.result)
            push!(records, (; name,method,repetition,seconds=run.seconds,config=Models.to_namedtuple(path),
                rows,all_roots,work=run.work,diagnostics=Models.to_namedtuple(run.result.diagnostics)))
            write_records(destination,records;anchor_probes)
            println("DONE ",name,"/",method," ",run.seconds," s, max roots ",maximum(length(r.roots) for r in all_roots))
            flush(stdout)
        end
    end
    return records
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PALCEndpointComparison.run_comparison(ARGS[1]; repetitions=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 3,
        scope=length(ARGS)>=3 ? ARGS[3] : "endpoints",
        cases_filter=length(ARGS)>=4 ? ARGS[4] : "all",probe=length(ARGS)>=5 && ARGS[5]=="true")
end
