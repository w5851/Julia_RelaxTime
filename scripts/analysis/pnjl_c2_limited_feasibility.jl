#!/usr/bin/env julia

"""Production-parity evaluator for the Issue #130 C2 limited feasibility.

Numerical jobs create the fine rho pool.  This evaluator reconstructs Stage A,
Stage B and the disjoint Stage-C pool from that input, then delegates every
classification, Maxwell and geometry decision to the calculation SHA's
existing Models helpers through the versioned Stage-C v2 evaluator module.
"""

using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(PROJECT_ROOT)
using CSV
using JSON3
using SHA

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const MODELS = Main.Models
const V2_PATH = joinpath(PROJECT_ROOT, "scripts", "analysis",
    "pnjl_stagec_density_certificate_feasibility_v2.jl")
if !isdefined(Main, :StageCDensityCertificateFeasibilityV2)
    include(V2_PATH)
end
const V2 = Main.StageCDensityCertificateFeasibilityV2

module C2LimitedFeasibility

using CSV
using JSON3
using SHA

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS = Main.Models
const V2 = Main.StageCDensityCertificateFeasibilityV2
const SCHEMA_VERSION = "pnjl_c2_limited_feasibility_v1"
const JOB_SCHEMA = "pnjl_c2_limited_feasibility_job_v1"
const CALCULATION_SHA_RE = r"^[0-9a-f]{40}$"
const XI_GRID = (-0.5, -0.35, -0.25, -0.2, -0.15, -0.1, 0.0, 0.3, 0.35, 0.5)
const DENSITY_ANCHORS = V2.DENSITY_ANCHORS
const FIRST_ORDER_CONTROLS = V2.FIRST_ORDER_CONTROLS
const MONOTONE_CONTROLS = V2.MONOTONE_CONTROLS
const ALL_ANCHORS = V2.ALL_ANCHORS
const ROUTES = V2.ROUTES
const CAPS = V2.CAPS
const STAGE_A_COARSE = 0.05
const STAGE_A_FINE = 0.025
const STAGE_B_COARSE = 0.0125
const STAGE_B_FINE = 0.00625
const STAGE_C_FINE = 0.003125

@inline function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && error("missing value for $(name)")
    args[index + 1]
end

@inline function _float(value, default=NaN)
    value === nothing || value === missing ? default : try
        result = Float64(value)
        isfinite(result) ? result : default
    catch
        default
    end
end

@inline function _field(row, name::Symbol, default=nothing)
    row === nothing && return default
    try
        if row isa AbstractDict
            return get(row, String(name), get(row, name, default))
        end
        hasproperty(row, name) ? getproperty(row, name) : default
    catch
        default
    end
end

function _sha256_file(path::String)
    bytes2hex(open(sha256, path))
end

function _manifest_paths(input_dir::String)
    paths = String[]
    for (root, _dirs, files) in walkdir(input_dir)
        "manifest.json" in files && push!(paths, joinpath(root, "manifest.json"))
    end
    sort!(paths)
end

function _required_file(manifest, manifest_path, name)
    files = _field(manifest, :files, nothing)
    declared = _field(files, Symbol(name), _field(files, name, nothing))
    declared === nothing && error("manifest missing $(name): $(manifest_path)")
    path = joinpath(dirname(manifest_path), name)
    isfile(path) || error("missing $(name): $(manifest_path)")
    actual = _sha256_file(path)
    String(declared) == actual || error("hash mismatch for $(name): $(manifest_path)")
    path
end

function _on_grid(rho::Float64, step::Float64)
    isapprox(rho / step, round(rho / step); atol=3e-7, rtol=0.0)
end

function _read_job(path, expected_sha, expected_postprocess, input_root)
    manifest = JSON3.read(read(path, String))
    String(_field(manifest, :schema_version, "")) == JOB_SCHEMA ||
        error("unexpected job schema $(path)")
    String(_field(manifest, :scope, "")) == "density" || error("job scope is not density")
    calculation_sha = String(_field(manifest, :calculation_sha, ""))
    calculation_sha == expected_sha || error("calculation SHA mismatch $(path)")
    postprocess_sha = String(_field(manifest, :postprocess_sha,
        _field(manifest, :workflow_head_sha, "")))
    isempty(expected_postprocess) || postprocess_sha == expected_postprocess ||
        error("postprocess SHA mismatch $(path)")
    xi = _float(_field(manifest, :xi))
    xi in XI_GRID || error("unexpected xi $(xi) in $(path)")
    pool_path = _required_file(manifest, path, "fine_pool.csv")
    slice_path = _required_file(manifest, path, "slice_metrics.csv")
    cost_path = _required_file(manifest, path, "method_costs.csv")
    points = V2.Point[]
    for row in CSV.File(pool_path)
        rho = _float(_field(row, :rho))
        T = _float(_field(row, :T_MeV))
        mu = _float(_field(row, :muq_MeV))
        residual = _float(_field(row, :residual_norm), Inf)
        finite = lowercase(String(_field(row, :finite, false))) in ("true", "1")
        converged = lowercase(String(_field(row, :converged, false))) in ("true", "1")
        finite && converged && all(isfinite, (rho, T, mu, residual)) ||
            error("invalid fine-pool point in $(path)")
        level = _on_grid(rho, STAGE_B_FINE) ? 0 : 1
        push!(points, V2.Point(xi, T, rho, mu, residual, "independent_oracle", level))
    end
    slices = Any[collect(CSV.File(slice_path))...]
    costs = Any[collect(CSV.File(cost_path))...]
    files = Any[(path=replace(relpath(file, input_root), '\\' => '/'), sha256=_sha256_file(file))
        for file in (pool_path, slice_path, cost_path)]
    push!(files, (path=replace(relpath(path, input_root), '\\' => '/'), sha256=_sha256_file(path)))
    (; manifest, xi, points, slices, costs, files, postprocess_sha)
end

function _production_points(points)
    result = V2.Point[]
    for point in points
        _on_grid(point.rho, STAGE_A_FINE) || continue
        level = _on_grid(point.rho, STAGE_A_COARSE) ? 0 : 1
        push!(result, V2.Point(point.xi, point.T_MeV, point.rho, point.muq_MeV,
            point.residual_norm, "production_hybrid", level))
    end
    result
end

function _full_curve(points, anchor)
    selected = [point for point in points if
        isapprox(point.xi, anchor[1]; atol=1e-8, rtol=0.0) &&
        isapprox(point.T_MeV, anchor[2]; atol=1e-8, rtol=0.0)]
    V2._curve(selected)
end

function _oracle_rows(points)
    rows = Any[]
    for anchor in ALL_ANCHORS
        curve = _full_curve(points, anchor)
        result = V2._evaluate(curve)
        status = result.status == :valid ? "confirmed_first_order" :
            (result.status == :invalid && result.reason == "no_s_shape" ?
                "confirmed_monotone" : "ambiguous_near_critical")
        push!(rows, (method="independent_oracle", xi=anchor[1], T_MeV=anchor[2], result_status=status))
    end
    rows
end

function _load_source(input_dir, expected_sha, expected_postprocess)
    paths = _manifest_paths(input_dir)
    length(paths) == length(XI_GRID) || error("expected $(length(XI_GRID)) xi jobs, got $(length(paths))")
    jobs = [_read_job(path, expected_sha, expected_postprocess, input_dir) for path in paths]
    Set(job.xi for job in jobs) == Set(XI_GRID) || error("xi job matrix is incomplete")
    points = V2.Point[]
    slices = Any[]
    costs = Any[]
    manifests = Any[]
    input_files = Any[]
    for job in jobs
        append!(points, job.points)
        append!(points, _production_points(job.points))
        append!(slices, job.slices)
        append!(costs, job.costs)
        push!(manifests, job.manifest)
        append!(input_files, job.files)
    end
    # Full-fine labels are generated only after all route inputs are assembled.
    append!(slices, _oracle_rows([point for point in points if point.method == "independent_oracle"]))
    # The fine pool is the numerical cost source.  The old dense baseline is a
    # fixed diagnostic comparator, not an extra solver call.
    push!(costs, (method="memoized_dense", unique_solves=length(ALL_ANCHORS) *
        Int(round(4.0 / STAGE_B_FINE + 1))))
    data = V2.SourceData(points, slices, costs, manifests, input_files,
        length(jobs), "success")
    data
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _frontier_rows(data)
    frontier = V2._frontier(data)
    frontier
end

function _verdict(frontier)
    rows = frontier.frontiers
    cap12 = [row for row in rows if row.cap == 12 && row.feasible]
    cap_high = [row for row in rows if row.cap in (16, 24) && row.feasible]
    any(row -> row.finite_failures > 0, rows) && return "density_solver_or_curve_failure"
    any(row -> row.unsupported_confirmations > 0, rows) && return "density_oracle_inconclusive"
    any(row -> row.multiple_candidate_anchors > 0, rows) && return "density_maxwell_candidate_inconclusive"
    !isempty(cap12) && return "density_feasible_candidate"
    !isempty(cap_high) && return "density_cap_contract_inconclusive"
    any(row -> !row.cost_gate, [row for row in rows if row.cap == 12]) &&
        return "density_performance_inconclusive"
    "density_solver_or_curve_failure"
end

function _git_head()
    try
        String(readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`))
    catch
        "unknown"
    end
end

function _write_outputs(output_dir, data, frontier, verdict, expected_sha,
        expected_postprocess, source_run_id, scope)
    mkpath(output_dir)
    _write_csv(joinpath(output_dir, "route_comparison.csv"), frontier.all_rows)
    _write_csv(joinpath(output_dir, "component_geometry.csv"), frontier.component_rows)
    _write_csv(joinpath(output_dir, "selected_point_index.csv"), frontier.selected_rows)
    _write_csv(joinpath(output_dir, "candidate_point_index.csv"), frontier.candidate_rows)
    _write_csv(joinpath(output_dir, "cost_frontier.csv"), frontier.frontiers)
    policy = V2._selected_policy(frontier, verdict == "density_feasible_candidate" ?
        "feasible_candidate" : verdict)
    policy["schema_version"] = SCHEMA_VERSION
    policy["density_scope"] = scope
    policy["cost_stop_runner_minutes"] = 150
    open(joinpath(output_dir, "selected_policy.json"), "w") do io
        JSON3.pretty(io, policy); write(io, '\n')
    end
    claim_rows = [
        (claim_id="density_classification", status=verdict == "density_feasible_candidate" ? "pass" : "inconclusive", evidence="route_comparison.csv"),
        (claim_id="density_geometry", status=verdict == "density_feasible_candidate" ? "pass" : "inconclusive", evidence="component_geometry.csv"),
        (claim_id="density_cost", status=verdict == "density_feasible_candidate" ? "pass" : "inconclusive", evidence="cost_frontier.csv"),
        (claim_id="production_promotion", status="not_claimed", evidence="requires subsequent production/shadow review"),
    ]
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), claim_rows)
    audit = joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl", "c2_convergence_audit_v1")
    cep = isfile(joinpath(audit, "tables", "c1_vs_c2_cep_gates.csv"))
    crossover = isfile(joinpath(audit, "tables", "c1_vs_c2_crossover_xi_0p2875.csv"))
    write(joinpath(output_dir, "README.md"), "# C2 limited feasibility\n\n" *
        "This package is diagnostic-only. Density numerical jobs produce the fine pool; " *
        "the evaluator is solver-free after that input. CEP/crossover remain separately gated.\n")
    files = Dict{String, String}()
    for (root, _dirs, names) in walkdir(output_dir)
        for name in names
            path = joinpath(root, name)
            endswith(name, "manifest.json") && continue
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha256_file(path)
        end
    end
    manifest = Dict(
        "schema_version" => SCHEMA_VERSION, "scope" => scope, "verdict" => verdict,
        "source_run_id" => source_run_id, "source_job_count" => data.job_count,
        "source_calculation_sha" => expected_sha,
        "source_postprocess_sha" => expected_postprocess,
        "producer_head_sha" => _git_head(), "solver_called" => false,
        "c2_audit_inputs" => Dict("cep_audit_present" => cep, "crossover_audit_present" => crossover),
        "selected_policy" => policy, "input_files" => data.input_files,
        "files" => files,
    )
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest); write(io, '\n')
    end
    manifest
end

function _parse_args(args)
    values = Dict{String, String}()
    index = 1
    while index <= length(args)
        arg = args[index]
        startswith(arg, "--") || error("unexpected argument $(arg)")
        key = replace(first(split(arg, "="; limit=2)), "--" => "")
        if occursin("=", arg)
            values[key] = split(arg, "="; limit=2)[2]
        else
            index == length(args) && error("missing value for --$(key)")
            values[key] = args[index + 1]; index += 1
        end
        index += 1
    end
    haskey(values, "input-dir") || error("--input-dir is required")
    values
end

function main(args=ARGS)
    options = _parse_args(args)
    scope = get(options, "scope", "density")
    scope == "density" || error("scope $(scope) is implemented as a dispatch contract but numerical density is the first runnable scope")
    expected_sha = lowercase(get(options, "expected-calculation-sha", ""))
    occursin(CALCULATION_SHA_RE, expected_sha) || error("expected calculation SHA must be 40 lowercase hex characters")
    expected_postprocess = get(options, "expected-source-postprocess-sha", "")
    source_run_id = get(options, "source-run-id", "")
    occursin(r"^\d+$", source_run_id) || error("source-run-id must be numeric")
    data = _load_source(abspath(options["input-dir"]), expected_sha, expected_postprocess)
    frontier = _frontier_rows(data)
    verdict = _verdict(frontier)
    output_dir = abspath(get(options, "output-dir", joinpath(PROJECT_ROOT, "c2_limited_aggregate")))
    manifest = _write_outputs(output_dir, data, frontier, verdict, expected_sha,
        expected_postprocess, source_run_id, scope)
    println(JSON3.write(Dict("verdict" => verdict, "source_job_count" => data.job_count,
        "solver_called" => false, "manifest_sha256" => _sha256_file(joinpath(output_dir, "manifest.json")))) )
    verdict == "density_feasible_candidate" ? 0 : 2
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(C2LimitedFeasibility.main(ARGS))
end
