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

"""Two natural-parameter sweeps, one root solve per point and no seed-catalog fallback."""
function dual_endpoint_path(model, path, work)
    started = time()
    options = Models._path_options(model, path; iterations=300)
    residual = Models._path_residual(model, path, options)
    values = sort!(unique(Models._path_values(path)))
    branches = Models.ContinuationBranch[]
    anchors = Models.AnchorRoot[]
    failures = 0
    messages = String[]
    for (direction, grid, side) in ((:forward, values, :hadron), (:reverse, reverse(values), :quark))
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
    diagnostics = Models.PathDiagnostics(:natural_parameter, :two_endpoints, :ground_state,
        time()-started, work.anchor_solves, 0, failures, 0, messages)
    return Models.PathSolveResult(path, selected.branches, selected.selections, anchors, diagnostics)
end

function run_method(model, path, method)
    work = Models.PALCWorkTelemetry()
    started = time_ns()
    result = if method === :dual_endpoint
        dual_endpoint_path(model, path, work)
    else
        values = Models._path_values(path)
        anchor = if method === :catalog
            Models.SeedPoolAnchor()
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

function write_records(destination, records)
    mkpath(dirname(destination))
    data = (; schema="palc_endpoint_comparison_v1", source_sha=strip(read(`git rev-parse HEAD`, String)),
        actions_run_id=get(ENV,"GITHUB_RUN_ID",""), julia_version=string(VERSION),
        p_num=24,t_num=8,iterations=300,residual_tolerance=1e-6,
        timing_contract="Same process, one thread, path solve plus identical selected-result reconstruction; validation and serialization excluded; first round includes compilation.",
        endpoint_contract="Hadron seed at the lower endpoint and quark seed at the upper endpoint; one direct root solve per target, stop a direction on failure, no catalog fallback; Newton jumps are retained as actual roots.",
        records)
    open(destination,"w") do io
        JSON.print(io,finite_json(data),2)
        println(io)
    end
end

function run_comparison(destination; repetitions=3)
    1 <= repetitions <= 3 || throw(ArgumentError("repetitions must be 1..3"))
    BLAS.set_num_threads(1)
    model = Models.create_model(:PNJL)
    cases = [
        (:tmu_t120, Models.FixedMuPath(120/HBARC, collect(260.0:2.0:380.0)./HBARC; p_num=24,t_num=8)),
        (:tmu_t120_fine, Models.FixedMuPath(120/HBARC, collect(296.0:0.25:318.0)./HBARC; p_num=24,t_num=8)),
        (:rho_t120, Models.FixedRhoPath(120/HBARC, collect(0.2:0.05:3.5); p_num=24,t_num=8)),
        (:tmu_t200, Models.FixedMuPath(200/HBARC, collect(260.0:2.0:380.0)./HBARC; p_num=24,t_num=8)),
    ]
    records = NamedTuple[]
    for (name,path) in cases, repetition in 1:repetitions
        for method in circshift(collect(METHODS), repetition-1)
            println("START ",name,"/",method," round ",repetition)
            flush(stdout)
            run = run_method(model,path,method)
            rows = [row_record(model,path,value,point) for (value,point) in zip(Models._path_values(path),run.selected)]
            all_roots = all_root_records(model,path,run.result)
            push!(records, (; name,method,repetition,seconds=run.seconds,config=Models.to_namedtuple(path),
                rows,all_roots,work=run.work,diagnostics=Models.to_namedtuple(run.result.diagnostics)))
            write_records(destination,records)
            println("DONE ",name,"/",method," ",run.seconds," s, max roots ",maximum(length(r.roots) for r in all_roots))
            flush(stdout)
        end
    end
    return records
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PALCEndpointComparison.run_comparison(ARGS[1]; repetitions=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 3)
end
