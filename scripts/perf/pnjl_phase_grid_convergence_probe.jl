#!/usr/bin/env julia

using Dates
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
using .Models

function _parse_output(args::Vector{String})
    isempty(args) && return nothing
    length(args) == 2 && args[1] == "--output" || error(
        "usage: julia --project=. scripts/perf/pnjl_phase_grid_convergence_probe.jl [--output <json>]",
    )
    return String(args[2])
end

function _run_case(name::String; rho_geometry_convergence::Bool, adaptive_temperature::Bool)
    output_dir = mktempdir()
    result = nothing
    wall_s = @elapsed result = Models.run_production_phase_pipeline(
        :PNJL;
        T_start=110.0,
        T_end=120.0,
        dT=10.0,
        rho_grid=collect(0.1:0.2:3.0),
        xi=0.0,
        output_dir=output_dir,
        profile=:perf_probe,
        solver_backend=:models,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=20,
        compute_crossover=false,
        cep_tol=0.5,
        cep_max_bisect_iter=2,
        cep_max_refine_level_rho=1,
        rho_geometry_convergence=rho_geometry_convergence,
        rho_position_tol_MeV=10.0,
        rho_density_tol=1.0,
        rho_maxwell_area_tol=1e-3,
        adaptive_temperature=adaptive_temperature,
        temperature_max_refine_level=1,
        promote_reference=false,
    )
    return Dict(
        "name" => name,
        "wall_s" => wall_s,
        "curve_count" => result.diagnostics["curve_count"],
        "scan_total" => result.diagnostics["scan_total"],
        "scan_failure" => result.diagnostics["scan_failure"],
        "boundary_count" => result.diagnostics["boundary_count"],
        "cep_eval_count" => result.diagnostics["cep_eval_count"],
    )
end

function main(args::Vector{String})
    output = _parse_output(args)

    # Compile the stable phase entrypoint before comparing the three policies.
    Models.run_production_phase_pipeline(
        :PNJL;
        T_start=120.0,
        T_end=120.0,
        dT=10.0,
        rho_grid=[0.1, 0.2, 0.3],
        output_dir=mktempdir(),
        profile=:perf_probe_warmup,
        solver_backend=:models,
        p_num=8,
        t_num=4,
        iterations=5,
        compute_crossover=false,
        cep_max_bisect_iter=1,
        cep_max_refine_level_rho=1,
        rho_geometry_convergence=false,
        adaptive_temperature=false,
        promote_reference=false,
    )

    cases = [
        _run_case("legacy_unknown_only"; rho_geometry_convergence=false, adaptive_temperature=false),
        _run_case("rho_geometry_convergence"; rho_geometry_convergence=true, adaptive_temperature=false),
        _run_case("rho_geometry_plus_adaptive_temperature"; rho_geometry_convergence=true, adaptive_temperature=true),
    ]
    payload = Dict(
        "generated_at" => string(now()),
        "scope" => "representative T=110,120 MeV first-order slice; p_num=8, t_num=4, iterations=20",
        "note" => "performance evidence only; not a physical convergence or production artifact",
        "cases" => cases,
    )

    rendered = sprint(io -> JSON3.pretty(io, payload))
    println(rendered)
    if output !== nothing
        mkpath(dirname(abspath(output)))
        open(output, "w") do io
            write(io, rendered)
            write(io, '\n')
        end
    end
end

main(collect(String.(ARGS)))
