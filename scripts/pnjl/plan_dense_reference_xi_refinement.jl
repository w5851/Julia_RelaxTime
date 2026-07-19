#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using CSV
using JSON3

if !isdefined(Main, :Models)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end
using Main.Models


Base.@kwdef mutable struct DenseXiRefinementPlanConfig
    reference_root::String = ""
    tag::String = "dense"
    intervals_path::String = ""
    level::Int = 1
    max_refine_level::Int = 1
    position_tol_MeV::Float64 = 0.10
    density_tol::Float64 = 0.01
    maxwell_area_tol::Float64 = 1e-4
    response_rtol::Float64 = 0.05
    matrix_output::String = ""
    intervals_output::String = ""
    records_output::String = ""
    summary_output::String = ""
end


function _require_value(args::Vector{String}, index::Int, flag::String)
    index == length(args) && error("missing value for $flag")
    return args[index + 1]
end


function parse_refinement_args(args::Vector{String})
    cfg = DenseXiRefinementPlanConfig()
    i = 1
    while i <= length(args)
        flag = args[i]
        value = _require_value(args, i, flag)
        if flag == "--reference-root"
            cfg.reference_root = value
        elseif flag == "--tag"
            cfg.tag = value
        elseif flag == "--intervals-path"
            cfg.intervals_path = value
        elseif flag == "--level"
            cfg.level = parse(Int, value)
        elseif flag == "--max-refine-level"
            cfg.max_refine_level = parse(Int, value)
        elseif flag == "--position-tol"
            cfg.position_tol_MeV = parse(Float64, value)
        elseif flag == "--density-tol"
            cfg.density_tol = parse(Float64, value)
        elseif flag == "--maxwell-area-tol"
            cfg.maxwell_area_tol = parse(Float64, value)
        elseif flag == "--response-rtol"
            cfg.response_rtol = parse(Float64, value)
        elseif flag == "--matrix-output"
            cfg.matrix_output = value
        elseif flag == "--intervals-output"
            cfg.intervals_output = value
        elseif flag == "--records-output"
            cfg.records_output = value
        elseif flag == "--summary-output"
            cfg.summary_output = value
        else
            error("unknown argument: $flag")
        end
        i += 2
    end

    isempty(cfg.reference_root) && error("--reference-root is required")
    isempty(cfg.intervals_path) && error("--intervals-path is required")
    isempty(cfg.matrix_output) && error("--matrix-output is required")
    isempty(cfg.intervals_output) && error("--intervals-output is required")
    isempty(cfg.records_output) && error("--records-output is required")
    isempty(cfg.summary_output) && error("--summary-output is required")
    cfg.level >= 1 || error("level must be positive")
    cfg.max_refine_level >= cfg.level || error("max refine level must be >= current level")
    Models._validate_phase_geometry_tolerances(Models.PhaseGeometryTolerances(
        position_MeV=cfg.position_tol_MeV,
        density=cfg.density_tol,
        maxwell_area=cfg.maxwell_area_tol,
        response_rtol=cfg.response_rtol,
    ))
    return cfg
end


@inline _xi_key(value::Real) = round(Float64(value); digits=12)


function _read_intervals(path::String)
    payload = JSON3.read(read(path, String))
    intervals = Tuple{Float64, Float64}[]
    for item in payload
        length(item) == 2 || error("each xi interval must contain exactly two endpoints")
        left = parse(Float64, string(item[1]))
        right = parse(Float64, string(item[2]))
        isfinite(left) && isfinite(right) && left < right ||
            error("invalid xi interval: ($left, $right)")
        push!(intervals, (_xi_key(left), _xi_key(right)))
    end
    sort!(intervals)
    unique!(intervals)
    return intervals
end


function _load_reference_results(reference_root::String, tag::String)
    manifest_path = joinpath(reference_root, "phase_reference_$(tag)_manifest.json")
    isfile(manifest_path) || error("missing merged reference manifest: $manifest_path")
    manifest = JSON3.read(read(manifest_path, String))
    xis = sort(unique(_xi_key(Float64(value)) for value in manifest.config.xi_values))

    boundary_by_xi = Dict(xi => NamedTuple[] for xi in xis)
    boundary_path = joinpath(reference_root, "boundary_$(tag).csv")
    if isfile(boundary_path)
        for row in CSV.File(boundary_path)
            xi = _xi_key(row.xi)
            hasproperty(row, :area_residual) ||
                error("boundary artifact lacks area_residual required for staged xi refinement")
            push!(get!(boundary_by_xi, xi, NamedTuple[]), (
                T_MeV=Float64(row.T_MeV),
                mu_transition_MeV=Float64(row.mu_transition_MeV),
                rho_hadron=Float64(row.rho_hadron),
                rho_quark=Float64(row.rho_quark),
                area_residual=Float64(row.area_residual),
                converged=Bool(row.converged),
            ))
        end
    end

    spinodal_by_xi = Dict(xi => NamedTuple[] for xi in xis)
    spinodal_path = joinpath(reference_root, "spinodals_$(tag).csv")
    if isfile(spinodal_path)
        for row in CSV.File(spinodal_path)
            xi = _xi_key(row.xi)
            push!(get!(spinodal_by_xi, xi, NamedTuple[]), (
                T_MeV=Float64(row.T_MeV),
                mu_spinodal_hadron_MeV=Float64(row.mu_spinodal_hadron_MeV),
                mu_spinodal_quark_MeV=Float64(row.mu_spinodal_quark_MeV),
                rho_spinodal_hadron=Float64(row.rho_spinodal_hadron),
                rho_spinodal_quark=Float64(row.rho_spinodal_quark),
            ))
        end
    end

    crossover_by_xi = Dict(xi => NamedTuple[] for xi in xis)
    crossover_path = joinpath(reference_root, "crossover_$(tag).csv")
    if isfile(crossover_path)
        for row in CSV.File(crossover_path)
            xi = _xi_key(row.xi)
            push!(get!(crossover_by_xi, xi, NamedTuple[]), (
                mu_MeV=Float64(row.mu_MeV),
                T_crossover_MeV=Float64(row.T_crossover_MeV),
                rho=Float64(row.rho),
                derivative=Float64(row.derivative),
            ))
        end
    end

    cep_by_xi = Dict(xi => Models.CEPResult() for xi in xis)
    cep_path = joinpath(reference_root, "cep_$(tag).csv")
    if isfile(cep_path)
        for row in CSV.File(cep_path)
            xi = _xi_key(row.xi)
            cep_by_xi[xi] = Models.CEPResult(
                found=true,
                T_cep_MeV=Float64(row.T_CEP_MeV),
                mu_cep_MeV=Float64(row.muq_CEP_MeV),
                uncertainty_T_MeV=Float64(row.uncertainty_T_MeV),
                T_bracket_low_MeV=Float64(row.T_bracket_low_MeV),
                T_bracket_high_MeV=Float64(row.T_bracket_high_MeV),
                bracket_width_T_MeV=Float64(row.bracket_width_T_MeV),
            )
        end
    end

    results = Dict{Float64, Models.PhasePipelineResult}()
    for xi in xis
        results[xi] = Models.PhasePipelineResult(
            xi=xi,
            cep=get(cep_by_xi, xi, Models.CEPResult()),
            first_order_boundary=get(boundary_by_xi, xi, NamedTuple[]),
            spinodal=get(spinodal_by_xi, xi, NamedTuple[]),
            crossover_line=get(crossover_by_xi, xi, NamedTuple[]),
        )
    end
    return results
end


@inline _csv_value(value) = value === nothing ? "" : string(value)


function _write_records(path::String, rows::Vector{NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "axis,xi,T_MeV,level,left,right,midpoint,position_error_MeV,density_error,maxwell_area,response_rtol,converged,reason")
        for row in rows
            values = (
                row.axis,
                row.xi,
                row.T_MeV,
                row.level,
                row.left,
                row.right,
                row.midpoint,
                row.position_error_MeV,
                row.density_error,
                row.maxwell_area,
                row.response_rtol,
                row.converged,
                row.reason,
            )
            println(io, join(_csv_value.(values), ','))
        end
    end
end


function build_refinement_plan(cfg::DenseXiRefinementPlanConfig)
    intervals = _read_intervals(cfg.intervals_path)
    results = isempty(intervals) ? Dict{Float64, Models.PhasePipelineResult}() :
        _load_reference_results(cfg.reference_root, cfg.tag)
    tol = Models.PhaseGeometryTolerances(
        position_MeV=cfg.position_tol_MeV,
        density=cfg.density_tol,
        maxwell_area=cfg.maxwell_area_tol,
        response_rtol=cfg.response_rtol,
    )

    records = NamedTuple[]
    next_intervals = Tuple{Float64, Float64}[]
    for (left_xi, right_xi) in intervals
        midpoint_xi = _xi_key(0.5 * (left_xi + right_xi))
        for xi in (left_xi, midpoint_xi, right_xi)
            haskey(results, xi) || error("missing staged xi artifact for xi=$xi at level $(cfg.level)")
        end
        error = Models._phase_result_midpoint_error(
            results[left_xi],
            results[midpoint_xi],
            results[right_xi],
            tol,
        )
        push!(records, (
            axis="xi",
            xi=midpoint_xi,
            T_MeV=nothing,
            level=cfg.level,
            left=left_xi,
            right=right_xi,
            midpoint=midpoint_xi,
            position_error_MeV=isfinite(error.position_MeV) ? error.position_MeV : nothing,
            density_error=isfinite(error.density) ? error.density : nothing,
            maxwell_area=isfinite(error.maxwell_area) ? error.maxwell_area : nothing,
            response_rtol=isfinite(error.response_rtol) ? error.response_rtol : nothing,
            converged=error.converged,
            reason=error.reason,
        ))
        if !error.converged && cfg.level < cfg.max_refine_level
            push!(next_intervals, (left_xi, midpoint_xi))
            push!(next_intervals, (midpoint_xi, right_xi))
        end
    end
    sort!(next_intervals)
    unique!(next_intervals)

    next_level = cfg.level + 1
    midpoints = sort(unique(_xi_key(0.5 * (left + right)) for (left, right) in next_intervals))
    include = if isempty(midpoints)
        [Dict("stage" => "level$(next_level)", "shard_id" => "l$(next_level)_none", "xi" => "0", "enabled" => false)]
    else
        [
            Dict(
                "stage" => "level$(next_level)",
                "shard_id" => "l$(next_level)_$(lpad(index - 1, 3, '0'))",
                "xi" => string(xi),
                "enabled" => true,
            ) for (index, xi) in enumerate(midpoints)
        ]
    end

    mkpath(dirname(cfg.matrix_output))
    write(cfg.matrix_output, JSON3.write(Dict("include" => include)))
    write(cfg.intervals_output, JSON3.write([[string(left), string(right)] for (left, right) in next_intervals]))
    _write_records(cfg.records_output, records)
    summary = Dict(
        "level" => cfg.level,
        "evaluated_interval_count" => length(intervals),
        "converged_interval_count" => count(row -> row.converged, records),
        "unconverged_interval_count" => count(row -> !row.converged, records),
        "next_interval_count" => length(next_intervals),
        "next_xi_count" => length(midpoints),
    )
    open(cfg.summary_output, "w") do io
        JSON3.pretty(io, summary)
        println(io)
    end
    return summary
end


function main(args=ARGS)
    cfg = parse_refinement_args(args)
    summary = build_refinement_plan(cfg)
    println("staged xi refinement level $(cfg.level): $(JSON3.write(summary))")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
