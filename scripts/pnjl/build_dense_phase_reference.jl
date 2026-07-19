#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using Printf
using JSON3

if !isdefined(Main, :Models)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end
using Main.Models

Base.@kwdef mutable struct DensePhaseReferenceConfig
    model_kind::Symbol = :PNJL
    mode::Symbol = :production
    xi_values::Vector{Float64} = collect(-0.5:0.05:0.5)
    T_min::Float64 = 60.0
    T_max::Float64 = 240.0
    T_step::Float64 = 5.0
    rho_min::Float64 = 0.0
    rho_max::Float64 = 4.0
    rho_step::Float64 = 0.05
    p_num::Int = 24
    t_num::Int = 8
    thermo_quadrature_policy::Symbol = :tensor_gauss
    thermo_quadrature_rtol::Float64 = 1e-8
    thermo_quadrature_atol::Float64 = 1e-10
    thermo_quadrature_maxevals::Int = 10^7
    iterations::Int = 80
    solver_backend::Symbol = :models
    seed_policy::Symbol = :hybrid_continuity
    profile::Symbol = :dense_reference
    tag::String = "dense"
    output_root::Union{Nothing, String} = nothing
    reference_root::Union{Nothing, String} = nothing
    overwrite::Bool = false
    compute_crossover::Bool = true
    crossover_method::Symbol = :peak
    crossover_variable::Symbol = :phi_u
    crossover_n_mu::Int = 16
    crossover_mu_max_MeV::Float64 = 450.0
    crossover_T_max_MeV::Float64 = NaN
    crossover_only::Bool = false
    crossover_mu_only_zero::Bool = false
    cep_tol_MeV::Float64 = 0.1
    rho_geometry_convergence::Bool = true
    rho_position_tol_MeV::Float64 = 0.05
    rho_density_tol::Float64 = 0.005
    rho_maxwell_area_tol::Float64 = 1e-4
    adaptive_temperature::Bool = true
    temperature_max_refine_level::Int = 2
    temperature_position_tol_MeV::Float64 = 0.10
    temperature_density_tol::Float64 = 0.01
    temperature_maxwell_area_tol::Float64 = 1e-4
    adaptive_xi::Bool = false
    xi_max_refine_level::Int = 2
    xi_position_tol_MeV::Float64 = 0.10
    xi_density_tol::Float64 = 0.01
    xi_maxwell_area_tol::Float64 = 1e-4
    xi_response_rtol::Float64 = 0.05
end

function usage()
    println("Usage: julia --project=. scripts/pnjl/build_dense_phase_reference.jl [options]")
    println("Options:")
    println("  --xi-list <csv>          xi list (default -0.5:0.05:0.5)")
    println("  --mode <name>            phase pipeline mode (default production)")
    println("  --xi-min <value>         xi min for ranged generation")
    println("  --xi-max <value>         xi max for ranged generation")
    println("  --xi-step <value>        xi step for ranged generation")
    println("  --T-min <MeV>            default 60")
    println("  --T-max <MeV>            default 240")
    println("  --T-step <MeV>           default 5")
    println("  --rho-min <value>        default 0.0")
    println("  --rho-max <value>        default 4.0")
    println("  --rho-step <value>       default 0.05")
    println("  --p-num <int>            gap/thermo momentum nodes (default 24)")
    println("  --t-num <int>            gap/thermo angle nodes (default 8)")
    println("  --thermo-quadrature-policy <name> tensor_gauss or rs_reduced_adaptive")
    println("  --thermo-quadrature-rtol <value>   adaptive relative tolerance (default 1e-8)")
    println("  --thermo-quadrature-atol <value>   adaptive absolute tolerance (default 1e-10)")
    println("  --thermo-quadrature-maxevals <int> adaptive evaluation cap (default 10000000)")
    println("  --iterations <int>       solver iteration cap (default 80)")
    println("  --tag <name>             output suffix, writes boundary_<tag>.csv etc. default dense")
    println("  --output-root <path>     processed run root")
    println("  --reference-root <path>  reference output directory")
    println("  --overwrite              overwrite existing aggregated outputs")
    println("  --no-crossover           skip crossover generation")
    println("  --crossover-n-mu <int>   crossover mu sampling count (default 16)")
    println("  --crossover-mu-max <MeV> crossover mu_q upper bound (default 450)")
    println("  --crossover-T-max <MeV>  explicit crossover search ceiling (default T-max)")
    println("  --cep-tol <MeV>          CEP temperature bracket width gate (default 0.1)")
    println("  --no-rho-geometry-convergence  disable coarse/fine Maxwell and spinodal gates")
    println("  --rho-position-tol <MeV> coarse/fine phase-position tolerance (default 0.05)")
    println("  --rho-density-tol <value> coarse/fine density tolerance (default 0.005)")
    println("  --rho-maxwell-area-tol <value> Maxwell diagnostic gate (default 1e-4)")
    println("  --no-adaptive-T          disable midpoint temperature refinement")
    println("  --T-refine-levels <int>  maximum adaptive temperature levels (default 2)")
    println("  --T-position-tol <MeV>   temperature interpolation position gate (default 0.10)")
    println("  --T-density-tol <value>  temperature interpolation density gate (default 0.01)")
    println("  --T-maxwell-area-tol <value> temperature Maxwell diagnostic gate (default 1e-4)")
    println("  --adaptive-xi            enable midpoint xi refinement (full reference only)")
    println("  --xi-refine-levels <int> maximum adaptive xi levels (default 2)")
    println("  --xi-position-tol <MeV>  xi interpolation position gate (default 0.10)")
    println("  --xi-density-tol <value> xi interpolation density gate (default 0.01)")
    println("  --xi-maxwell-area-tol <value> xi midpoint Maxwell diagnostic gate (default 1e-4)")
    println("  --xi-response-rtol <value> crossover response derivative gate (default 0.05)")
    println("  --crossover-only         skip phase pipeline; only generate crossover reference")
    println("  --crossover-mu0-only     with --crossover-only, compute only mu=0 crossover point")
    println("  -h, --help               show help")
end

function parse_float_list(raw::AbstractString)
    vals = Float64[]
    for token in split(raw, ',')
        s = strip(token)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && error("empty numeric list")
    return vals
end

function inclusive_step_grid(start::Real, stop::Real, step::Real; axis::AbstractString="grid")
    first_value = Float64(start)
    last_value = Float64(stop)
    spacing = Float64(step)
    all(isfinite, (first_value, last_value, spacing)) ||
        error("$axis bounds and step must be finite")
    spacing > 0 || error("$axis step must be positive")
    first_value <= last_value || error("$axis minimum must be <= maximum")

    values = collect(first_value:spacing:last_value)
    tolerance = 64 * eps(max(abs(first_value), abs(last_value), abs(spacing), 1.0))
    if last(values) < last_value - tolerance
        push!(values, last_value)
    else
        values[end] = last_value
    end
    return values
end

function parse_args(args::Vector{String})
    cfg = DensePhaseReferenceConfig()
    ranged = false
    xi_min = first(cfg.xi_values)
    xi_max = last(cfg.xi_values)
    xi_step = 0.05

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            value = args[i + 1]
            i += 1
            return value
        end

        if arg == "--xi-list"
            cfg.xi_values = unique(sort(parse_float_list(require_value())))
        elseif arg == "--mode"
            cfg.mode = Symbol(lowercase(require_value()))
        elseif arg == "--xi-min"
            xi_min = parse(Float64, require_value())
            ranged = true
        elseif arg == "--xi-max"
            xi_max = parse(Float64, require_value())
            ranged = true
        elseif arg == "--xi-step"
            xi_step = parse(Float64, require_value())
            ranged = true
        elseif arg == "--T-min"
            cfg.T_min = parse(Float64, require_value())
        elseif arg == "--T-max"
            cfg.T_max = parse(Float64, require_value())
        elseif arg == "--T-step"
            cfg.T_step = parse(Float64, require_value())
        elseif arg == "--rho-min"
            cfg.rho_min = parse(Float64, require_value())
        elseif arg == "--rho-max"
            cfg.rho_max = parse(Float64, require_value())
        elseif arg == "--rho-step"
            cfg.rho_step = parse(Float64, require_value())
        elseif arg == "--p-num"
            cfg.p_num = parse(Int, require_value())
        elseif arg == "--t-num"
            cfg.t_num = parse(Int, require_value())
        elseif arg == "--thermo-quadrature-policy"
            cfg.thermo_quadrature_policy = Symbol(lowercase(require_value()))
        elseif arg == "--thermo-quadrature-rtol"
            cfg.thermo_quadrature_rtol = parse(Float64, require_value())
        elseif arg == "--thermo-quadrature-atol"
            cfg.thermo_quadrature_atol = parse(Float64, require_value())
        elseif arg == "--thermo-quadrature-maxevals"
            cfg.thermo_quadrature_maxevals = parse(Int, require_value())
        elseif arg == "--iterations"
            cfg.iterations = parse(Int, require_value())
        elseif arg == "--tag"
            cfg.tag = String(require_value())
        elseif arg == "--output-root"
            cfg.output_root = String(require_value())
        elseif arg == "--reference-root"
            cfg.reference_root = String(require_value())
        elseif arg == "--overwrite"
            cfg.overwrite = true
        elseif arg == "--no-crossover"
            cfg.compute_crossover = false
        elseif arg == "--crossover-n-mu"
            cfg.crossover_n_mu = parse(Int, require_value())
        elseif arg == "--crossover-mu-max"
            cfg.crossover_mu_max_MeV = parse(Float64, require_value())
        elseif arg == "--crossover-T-max"
            cfg.crossover_T_max_MeV = parse(Float64, require_value())
        elseif arg == "--cep-tol"
            cfg.cep_tol_MeV = parse(Float64, require_value())
        elseif arg == "--no-rho-geometry-convergence"
            cfg.rho_geometry_convergence = false
        elseif arg == "--rho-position-tol"
            cfg.rho_position_tol_MeV = parse(Float64, require_value())
        elseif arg == "--rho-density-tol"
            cfg.rho_density_tol = parse(Float64, require_value())
        elseif arg == "--rho-maxwell-area-tol"
            cfg.rho_maxwell_area_tol = parse(Float64, require_value())
        elseif arg == "--no-adaptive-T"
            cfg.adaptive_temperature = false
        elseif arg == "--T-refine-levels"
            cfg.temperature_max_refine_level = parse(Int, require_value())
        elseif arg == "--T-position-tol"
            cfg.temperature_position_tol_MeV = parse(Float64, require_value())
        elseif arg == "--T-density-tol"
            cfg.temperature_density_tol = parse(Float64, require_value())
        elseif arg == "--T-maxwell-area-tol"
            cfg.temperature_maxwell_area_tol = parse(Float64, require_value())
        elseif arg == "--adaptive-xi"
            cfg.adaptive_xi = true
        elseif arg == "--xi-refine-levels"
            cfg.xi_max_refine_level = parse(Int, require_value())
        elseif arg == "--xi-position-tol"
            cfg.xi_position_tol_MeV = parse(Float64, require_value())
        elseif arg == "--xi-density-tol"
            cfg.xi_density_tol = parse(Float64, require_value())
        elseif arg == "--xi-maxwell-area-tol"
            cfg.xi_maxwell_area_tol = parse(Float64, require_value())
        elseif arg == "--xi-response-rtol"
            cfg.xi_response_rtol = parse(Float64, require_value())
        elseif arg == "--crossover-only"
            cfg.crossover_only = true
        elseif arg == "--crossover-mu0-only"
            cfg.crossover_mu_only_zero = true
        elseif arg in ("-h", "--help")
            usage()
            exit(0)
        else
            error("unknown option: $arg")
        end

        i += 1
    end

    if ranged
        cfg.xi_values = inclusive_step_grid(xi_min, xi_max, xi_step; axis="xi")
    end

    cfg.T_min > 0 || error("T-min must be positive for the five-variable phase solve")
    cfg.T_min <= cfg.T_max || error("T-min must be <= T-max")
    cfg.T_step > 0 || error("T-step must be positive")
    cfg.rho_min <= cfg.rho_max || error("rho-min must be <= rho-max")
    cfg.rho_step > 0 || error("rho-step must be positive")
    cfg.p_num > 0 || error("p-num must be positive")
    cfg.t_num > 0 || error("t-num must be positive")
    Models.PNJLIntegrals.validate_thermal_quadrature_policy(cfg.thermo_quadrature_policy)
    Models.PNJLIntegrals.validate_thermal_quadrature_controls(
        cfg.thermo_quadrature_rtol,
        cfg.thermo_quadrature_atol,
        cfg.thermo_quadrature_maxevals,
    )
    cfg.iterations > 0 || error("iterations must be positive")
    cfg.crossover_n_mu > 0 || error("crossover-n-mu must be positive")
    cfg.crossover_mu_max_MeV > 0 || error("crossover-mu-max must be positive")
    cfg.cep_tol_MeV > 0 || error("cep-tol must be positive")
    cfg.temperature_max_refine_level >= 0 || error("T-refine-levels must be nonnegative")
    cfg.xi_max_refine_level >= 0 || error("xi-refine-levels must be nonnegative")
    cfg.rho_position_tol_MeV > 0 || error("rho-position-tol must be positive")
    cfg.rho_density_tol > 0 || error("rho-density-tol must be positive")
    cfg.rho_maxwell_area_tol > 0 || error("rho-maxwell-area-tol must be positive")
    cfg.temperature_position_tol_MeV > 0 || error("T-position-tol must be positive")
    cfg.temperature_density_tol > 0 || error("T-density-tol must be positive")
    cfg.temperature_maxwell_area_tol > 0 || error("T-maxwell-area-tol must be positive")
    cfg.xi_position_tol_MeV > 0 || error("xi-position-tol must be positive")
    cfg.xi_density_tol > 0 || error("xi-density-tol must be positive")
    cfg.xi_maxwell_area_tol > 0 || error("xi-maxwell-area-tol must be positive")
    cfg.xi_response_rtol > 0 || error("xi-response-rtol must be positive")
    cfg.adaptive_xi && cfg.crossover_only && error("adaptive xi refinement requires full phase reference mode")
    if !isnan(cfg.crossover_T_max_MeV)
        isfinite(cfg.crossover_T_max_MeV) || error("crossover-T-max must be finite")
        cfg.crossover_T_max_MeV >= cfg.T_min || error("crossover-T-max must be >= T-min")
    end
    isempty(cfg.xi_values) && error("xi grid cannot be empty")
    return cfg
end

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function current_git_commit()
    root = project_root()
    try
        commit = readchomp(`git -C $root rev-parse HEAD`)
        return isempty(strip(commit)) ? nothing : String(strip(commit))
    catch
        return nothing
    end
end

@inline function _norm_slash(path::AbstractString)
    return replace(String(path), '\\' => '/')
end

function _repo_relpath(path::AbstractString)
    root = project_root()
    abs_path = normpath(abspath(String(path)))
    rel = try
        relpath(abs_path, root)
    catch
        nothing
    end
    if rel !== nothing
        return _norm_slash(String(rel))
    end
    return _norm_slash(abs_path)
end

function xi_token(xi::Float64)
    return replace(@sprintf("%.6f", xi), "." => "p", "-" => "m")
end

@inline resolved_crossover_T_max_MeV(cfg::DensePhaseReferenceConfig) =
    isfinite(cfg.crossover_T_max_MeV) ? cfg.crossover_T_max_MeV : cfg.T_max

function ensure_writable(path::String, overwrite::Bool)
    if isfile(path) && !overwrite
        error("output exists: $path; rerun with --overwrite")
    end
end

function write_boundary_csv(path::String, rows)
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark,area_residual,converged,curve_parameter,plot_order_key")
        for row in rows
            println(io, "$(row.xi),$(row.T_MeV),$(row.mu_transition_MeV),$(row.rho_hadron),$(row.rho_quark),$(row.area_residual),$(row.converged),$(row.T_MeV),$(row.T_MeV)")
        end
    end
end

function write_cep_csv(path::String, rows)
    open(path, "w") do io
        println(io, "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV,uncertainty_T_MeV,T_bracket_low_MeV,T_bracket_high_MeV,bracket_width_T_MeV")
        for row in rows
            println(io, "$(row.xi),$(row.T_CEP_MeV),$(row.muq_CEP_MeV),$(row.muB_CEP_MeV),$(row.uncertainty_T_MeV),$(row.T_bracket_low_MeV),$(row.T_bracket_high_MeV),$(row.bracket_width_T_MeV)")
        end
    end
end

function write_spinodal_csv(path::String, rows)
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark,curve_parameter,plot_order_key")
        for row in rows
            println(io, "$(row.xi),$(row.T_MeV),$(row.mu_spinodal_hadron_MeV),$(row.mu_spinodal_quark_MeV),$(row.rho_spinodal_hadron),$(row.rho_spinodal_quark),$(row.T_MeV),$(row.T_MeV)")
        end
    end
end

function write_crossover_csv(path::String, rows)
    open(path, "w") do io
        println(io, "xi,mu_MeV,T_crossover_MeV,rho,method,converged,derivative,variable,curve_parameter,plot_order_key")
        for row in rows
            println(io, "$(row.xi),$(row.mu_MeV),$(row.T_crossover_MeV),$(row.rho),$(row.method),$(row.converged),$(row.derivative),$(row.variable),$(row.mu_MeV),$(row.mu_MeV)")
        end
    end
end

@inline _dense_record_value(row, key::Symbol, default=nothing) =
    hasproperty(row, key) ? getproperty(row, key) : default

@inline _dense_csv_value(value) = value === nothing ? "" : string(value)
@inline _dense_record_with_xi(row::NamedTuple, xi::Real) = merge(row, (xi=Float64(xi),))

function write_grid_convergence_csv(path::String, rows)
    open(path, "w") do io
        println(io, "axis,xi,T_MeV,level,left,right,midpoint,position_error_MeV,density_error,maxwell_area,response_rtol,converged,reason")
        for row in rows
            values = (
                _dense_record_value(row, :axis, ""),
                _dense_record_value(row, :xi),
                _dense_record_value(row, :T_MeV),
                _dense_record_value(row, :level),
                _dense_record_value(row, :left),
                _dense_record_value(row, :right),
                _dense_record_value(row, :midpoint),
                _dense_record_value(row, :position_error_MeV),
                _dense_record_value(row, :density_error),
                _dense_record_value(row, :maxwell_area),
                _dense_record_value(row, :response_rtol),
                _dense_record_value(row, :converged, false),
                _dense_record_value(row, :reason, ""),
            )
            println(io, join(_dense_csv_value.(values), ','))
        end
    end
end

function _crossover_column_definitions()
    return [
        Dict("name" => "xi", "type" => "Float64", "unit" => "dimensionless", "description" => "anisotropy control parameter"),
        Dict("name" => "mu_MeV", "type" => "Float64", "unit" => "MeV", "description" => "quark chemical potential sample (mu_q) on the crossover scan"),
        Dict("name" => "T_crossover_MeV", "type" => "Float64", "unit" => "MeV", "description" => "detected crossover temperature at the sampled xi and mu_q"),
        Dict("name" => "rho", "type" => "Float64", "unit" => "rho0", "description" => "density-like coordinate returned by the crossover builder for the sampled point"),
        Dict("name" => "method", "type" => "String", "unit" => nothing, "description" => "crossover detector method"),
        Dict("name" => "converged", "type" => "Bool", "unit" => nothing, "description" => "whether the crossover detector reported a valid point"),
        Dict("name" => "derivative", "type" => "Float64", "unit" => nothing, "description" => "peak/response diagnostic returned by the detector"),
        Dict("name" => "variable", "type" => "String", "unit" => nothing, "description" => "order-parameter-like variable used by the detector"),
        Dict("name" => "curve_parameter", "type" => "Float64", "unit" => "MeV", "description" => "physical parameter used to connect this curve in plots; for crossover this is mu_q"),
        Dict("name" => "plot_order_key", "type" => "Float64", "unit" => "MeV", "description" => "stable sort key for plotting this curve; for crossover this is mu_q"),
    ]
end

function _crossover_dense_meaning(cfg::DensePhaseReferenceConfig)
    xi_diffs = length(cfg.xi_values) > 1 ? diff(cfg.xi_values) : Float64[]
    xi_uniform = !isempty(xi_diffs) && all(isapprox(d, xi_diffs[1]; atol=1e-10, rtol=1e-10) for d in xi_diffs)
    xi_strategy = xi_uniform ? "uniform_grid" : "explicit_anchor_list"
    xi_step = xi_uniform ? xi_diffs[1] : nothing

    if cfg.crossover_mu_only_zero
        return Dict(
            "kind" => "mu0_only",
            "description" => "dense in xi only; for each xi solve the crossover point at mu_q = 0 MeV",
            "xi_sampling" => Dict(
                "strategy" => xi_strategy,
                "count" => length(cfg.xi_values),
                "step" => xi_step,
                "values" => cfg.xi_values,
            ),
            "mu_q_sampling" => Dict(
                "strategy" => "fixed_single_value",
                "values_MeV" => [0.0],
            ),
        )
    end

    mu_samples = collect(range(0.0; stop=cfg.crossover_mu_max_MeV, length=cfg.crossover_n_mu))
    return Dict(
        "kind" => "xi_mu_dense_reference",
        "description" => "uniform xi grid combined with uniform mu_q sampling; each (xi, mu_q) slice runs crossover detection on a T window and stores the detected T_crossover",
        "xi_sampling" => Dict(
            "strategy" => xi_strategy,
            "count" => length(cfg.xi_values),
            "step" => xi_step,
            "values" => cfg.xi_values,
        ),
        "mu_q_sampling" => Dict(
            "strategy" => "uniform_grid",
            "count" => length(mu_samples),
            "min_MeV" => first(mu_samples),
            "max_MeV" => last(mu_samples),
            "values_MeV" => mu_samples,
        ),
        "temperature_window_MeV" => Dict(
            "min" => cfg.T_min,
            "max" => resolved_crossover_T_max_MeV(cfg),
        ),
    )
end

function write_crossover_meta(path::String, cfg::DensePhaseReferenceConfig, rows, crossover_csv_path::String)
    xi_values = sort(unique(Float64[row.xi for row in rows]))
    mu_values = sort(unique(Float64[row.mu_MeV for row in rows]))

    payload = Dict(
        "schema_version" => "v1",
        "artifact" => Dict(
            "path" => _repo_relpath(crossover_csv_path),
            "row_count" => length(rows),
        ),
        "generator" => Dict(
            "script" => _repo_relpath(joinpath(@__DIR__, "build_dense_phase_reference.jl")),
            "git_commit" => current_git_commit(),
            "generated_at" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        ),
        "xi_coverage" => Dict(
            "count" => length(xi_values),
            "min" => isempty(xi_values) ? nothing : first(xi_values),
            "max" => isempty(xi_values) ? nothing : last(xi_values),
            "values" => xi_values,
        ),
        "mu_q_coverage" => Dict(
            "count" => length(mu_values),
            "min_MeV" => isempty(mu_values) ? nothing : first(mu_values),
            "max_MeV" => isempty(mu_values) ? nothing : last(mu_values),
            "values_MeV" => mu_values,
        ),
        "column_definitions" => _crossover_column_definitions(),
        "dense_meaning" => _crossover_dense_meaning(cfg),
    )

    open(path, "w") do io
        JSON3.pretty(io, payload)
    end
end

function manifest_generator_payload(cfg::DensePhaseReferenceConfig)
    return Dict(
        "script" => _repo_relpath(joinpath(@__DIR__, "build_dense_phase_reference.jl")),
        "git_commit" => current_git_commit(),
        "generated_at" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        "crossover_only" => cfg.crossover_only,
        "crossover_mu0_only" => cfg.crossover_mu_only_zero,
    )
end

function manifest_config_payload(cfg::DensePhaseReferenceConfig)
    return Dict(
        "tag" => cfg.tag,
        "xi_values" => cfg.xi_values,
        "T_min_MeV" => cfg.T_min,
        "T_max_MeV" => cfg.T_max,
        "T_step_MeV" => cfg.T_step,
        "rho_min" => cfg.rho_min,
        "rho_max" => cfg.rho_max,
        "rho_step" => cfg.rho_step,
        "p_num" => cfg.p_num,
        "t_num" => cfg.t_num,
        "thermo_quadrature_policy" => String(cfg.thermo_quadrature_policy),
        "thermo_quadrature_rtol" => cfg.thermo_quadrature_rtol,
        "thermo_quadrature_atol" => cfg.thermo_quadrature_atol,
        "thermo_quadrature_maxevals" => cfg.thermo_quadrature_maxevals,
        "iterations" => cfg.iterations,
        "mode" => String(cfg.mode),
        "compute_crossover" => cfg.compute_crossover,
        "crossover_method" => String(cfg.crossover_method),
        "crossover_variable" => String(cfg.crossover_variable),
        "crossover_n_mu" => cfg.crossover_n_mu,
        "crossover_mu_max_MeV" => cfg.crossover_mu_max_MeV,
        "crossover_T_max_MeV" => resolved_crossover_T_max_MeV(cfg),
        "crossover_only" => cfg.crossover_only,
        "crossover_mu0_only" => cfg.crossover_mu_only_zero,
        "cep_tol_MeV" => cfg.cep_tol_MeV,
        "rho_geometry_convergence" => cfg.rho_geometry_convergence,
        "rho_position_tol_MeV" => cfg.rho_position_tol_MeV,
        "rho_density_tol" => cfg.rho_density_tol,
        "rho_maxwell_area_tol" => cfg.rho_maxwell_area_tol,
        "adaptive_temperature" => cfg.adaptive_temperature,
        "temperature_max_refine_level" => cfg.temperature_max_refine_level,
        "temperature_position_tol_MeV" => cfg.temperature_position_tol_MeV,
        "temperature_density_tol" => cfg.temperature_density_tol,
        "temperature_maxwell_area_tol" => cfg.temperature_maxwell_area_tol,
        "adaptive_xi" => cfg.adaptive_xi,
        "xi_max_refine_level" => cfg.xi_max_refine_level,
        "xi_position_tol_MeV" => cfg.xi_position_tol_MeV,
        "xi_density_tol" => cfg.xi_density_tol,
        "xi_maxwell_area_tol" => cfg.xi_maxwell_area_tol,
        "xi_response_rtol" => cfg.xi_response_rtol,
    )
end

function build_crossover_only_rows(cfg::DensePhaseReferenceConfig, xi::Float64)
    rows = NamedTuple[]
    if cfg.crossover_mu_only_zero
        T_min_fm = cfg.T_min / 197.327
        T_max_fm = resolved_crossover_T_max_MeV(cfg) / 197.327
        result = Models.detect_crossover(
            0.0,
            (T_min_fm, T_max_fm);
            method=cfg.crossover_method,
            variable=cfg.crossover_variable,
            xi=xi,
            model_kind=cfg.model_kind,
            solver_backend=cfg.solver_backend,
            p_num=cfg.p_num,
            t_num=cfg.t_num,
            thermo_quadrature_policy=cfg.thermo_quadrature_policy,
            thermo_quadrature_rtol=cfg.thermo_quadrature_rtol,
            thermo_quadrature_atol=cfg.thermo_quadrature_atol,
            thermo_quadrature_maxevals=cfg.thermo_quadrature_maxevals,
        )
        T_mev = result.T_crossover === nothing ? NaN : Float64(result.T_crossover) * 197.327
        rho = result.rho === nothing ? NaN : Float64(result.rho)
        derivative = result.derivative_value === nothing ? NaN : Float64(result.derivative_value)
        push!(rows, (
            xi=xi,
            mu_MeV=0.0,
            T_crossover_MeV=T_mev,
            rho=rho,
            method=String(cfg.crossover_method),
            converged=Bool(result.found),
            derivative=derivative,
            variable=String(cfg.crossover_variable),
        ))
        return rows
    end

    local_rows = Models.build_crossover_line(
        ;
        mu_max_MeV=cfg.crossover_mu_max_MeV,
        T_min_MeV=cfg.T_min,
        T_max_MeV=resolved_crossover_T_max_MeV(cfg),
        xi=xi,
        n_mu=cfg.crossover_n_mu,
        method=cfg.crossover_method,
        variable=cfg.crossover_variable,
        model_kind=cfg.model_kind,
        solver_backend=cfg.solver_backend,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        thermo_quadrature_policy=cfg.thermo_quadrature_policy,
        thermo_quadrature_rtol=cfg.thermo_quadrature_rtol,
        thermo_quadrature_atol=cfg.thermo_quadrature_atol,
        thermo_quadrature_maxevals=cfg.thermo_quadrature_maxevals,
    )
    for row in local_rows
        push!(rows, merge((xi=xi,), row))
    end
    return rows
end

function _adaptive_xi_refinement!(
        initial_xi::Vector{Float64},
        evaluate_xi::Function,
        cfg::DensePhaseReferenceConfig)
    resolved = sort(unique(copy(initial_xi)))
    cfg.adaptive_xi || return resolved, NamedTuple[]
    cfg.xi_max_refine_level > 0 || return resolved, NamedTuple[]

    tol = Models.PhaseGeometryTolerances(
        position_MeV=cfg.xi_position_tol_MeV,
        density=cfg.xi_density_tol,
        maxwell_area=cfg.xi_maxwell_area_tol,
        response_rtol=cfg.xi_response_rtol,
    )
    intervals = Tuple{Float64, Float64}[
        (resolved[i], resolved[i + 1]) for i in 1:(length(resolved) - 1)
    ]
    records = NamedTuple[]

    for level in 1:cfg.xi_max_refine_level
        isempty(intervals) && break
        midpoints = sort(unique(Float64[0.5 * (left + right) for (left, right) in intervals]))
        for xi in midpoints
            evaluate_xi(xi)
        end

        next_intervals = Tuple{Float64, Float64}[]
        for (left_xi, right_xi) in intervals
            midpoint_xi = 0.5 * (left_xi + right_xi)
            error = Models._phase_result_midpoint_error(
                evaluate_xi(left_xi),
                evaluate_xi(midpoint_xi),
                evaluate_xi(right_xi),
                tol,
            )
            push!(records, (
                axis="xi",
                xi=midpoint_xi,
                T_MeV=nothing,
                level=level,
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
            push!(resolved, midpoint_xi)
            if !error.converged
                push!(next_intervals, (left_xi, midpoint_xi))
                push!(next_intervals, (midpoint_xi, right_xi))
            end
        end
        sort!(resolved)
        unique!(resolved)
        intervals = next_intervals
    end
    return resolved, records
end

function build_outputs(cfg::DensePhaseReferenceConfig)
    root = project_root()
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_root = isnothing(cfg.output_root) ? joinpath(root, "data", "processed", "pnjl", "phase_diagram", "dense_reference_$(timestamp)") : cfg.output_root
    reference_root = isnothing(cfg.reference_root) ? joinpath(root, "data", "reference", "pnjl") : cfg.reference_root
    mkpath(output_root)
    mkpath(reference_root)

    boundary_rows = NamedTuple[]
    cep_rows = NamedTuple[]
    spinodal_rows = NamedTuple[]
    crossover_rows = NamedTuple[]
    grid_convergence_rows = NamedTuple[]
    requested_xi_values = sort(unique(copy(cfg.xi_values)))
    manifest = Dict{String, Any}(
        "schema_version" => "v1",
        "generator" => manifest_generator_payload(cfg),
        "config" => manifest_config_payload(cfg),
        "output_root" => _repo_relpath(output_root),
        "reference_root" => _repo_relpath(reference_root),
        "artifacts" => Dict{String, Any}(),
        "runs" => Any[],
    )

    T_grid = inclusive_step_grid(cfg.T_min, cfg.T_max, cfg.T_step; axis="temperature")
    rho_grid = inclusive_step_grid(cfg.rho_min, cfg.rho_max, cfg.rho_step; axis="rho")

    if cfg.crossover_only
        for xi in requested_xi_values
            println("phase reference xi=$(xi): crossover-only T=$(cfg.T_min):$(resolved_crossover_T_max_MeV(cfg))")
            flush(stdout)
            run_dir = joinpath(output_root, "xi_$(xi_token(xi))")
            local_rows = cfg.compute_crossover ? build_crossover_only_rows(cfg, xi) : NamedTuple[]
            append!(crossover_rows, local_rows)
            println("phase reference xi=$(xi): crossover-only rows=$(length(local_rows))")
            flush(stdout)
            push!(manifest["runs"], Dict(
                "xi" => xi,
                "run_id" => "crossover_only",
                "run_dir" => _repo_relpath(run_dir),
                "boundary_count" => 0,
                "spinodal_count" => 0,
                "crossover_count" => length(local_rows),
                "cep_found" => false,
            ))
        end
    else
        result_cache = Dict{Float64, Models.PhasePipelineResult}()
        function evaluate_xi(xi_value::Float64)
            key = round(Float64(xi_value); digits=12)
            haskey(result_cache, key) && return result_cache[key]
            println("phase reference xi=$(key): start T=$(cfg.T_min):$(cfg.T_step):$(cfg.T_max), rho=$(cfg.rho_min):$(cfg.rho_step):$(cfg.rho_max)")
            flush(stdout)
            run_dir = joinpath(output_root, "xi_$(xi_token(key))")
            result = Models.run_phase_pipeline(
                cfg.model_kind;
                mode=cfg.mode,
                T_grid=T_grid,
                rho_grid=rho_grid,
                xi=key,
                output_dir=run_dir,
                run_id="dense_$(cfg.tag)_xi_$(xi_token(key))",
                profile=cfg.profile,
                solver_backend=cfg.solver_backend,
                seed_policy=cfg.seed_policy,
                p_num=cfg.p_num,
                t_num=cfg.t_num,
                thermo_quadrature_policy=cfg.thermo_quadrature_policy,
                thermo_quadrature_rtol=cfg.thermo_quadrature_rtol,
                thermo_quadrature_atol=cfg.thermo_quadrature_atol,
                thermo_quadrature_maxevals=cfg.thermo_quadrature_maxevals,
                iterations=cfg.iterations,
                compute_crossover=cfg.compute_crossover,
                crossover_method=cfg.crossover_method,
                crossover_variable=cfg.crossover_variable,
                crossover_n_mu=cfg.crossover_n_mu,
                crossover_mu0_only=cfg.crossover_mu_only_zero,
                crossover_T_max_MeV=resolved_crossover_T_max_MeV(cfg),
                cep_tol=cfg.cep_tol_MeV,
                rho_geometry_convergence=cfg.rho_geometry_convergence,
                rho_position_tol_MeV=cfg.rho_position_tol_MeV,
                rho_density_tol=cfg.rho_density_tol,
                rho_maxwell_area_tol=cfg.rho_maxwell_area_tol,
                adaptive_temperature=cfg.adaptive_temperature,
                temperature_max_refine_level=cfg.temperature_max_refine_level,
                temperature_position_tol_MeV=cfg.temperature_position_tol_MeV,
                temperature_density_tol=cfg.temperature_density_tol,
                temperature_maxwell_area_tol=cfg.temperature_maxwell_area_tol,
                promote_reference=false,
            )
            result_cache[key] = result
            println("phase reference xi=$(key): boundary=$(length(result.first_order_boundary)) spinodal=$(length(result.spinodal)) crossover=$(length(result.crossover_line)) cep_found=$(result.cep.found)")
            flush(stdout)
            return result
        end

        for xi in requested_xi_values
            evaluate_xi(xi)
        end
        resolved_xi_values, xi_convergence_records =
            _adaptive_xi_refinement!(requested_xi_values, evaluate_xi, cfg)
        append!(grid_convergence_rows, xi_convergence_records)
        cfg.xi_values = resolved_xi_values

        for xi in resolved_xi_values
            key = round(Float64(xi); digits=12)
            result = evaluate_xi(key)
            run_dir = joinpath(output_root, "xi_$(xi_token(key))")
            for row in result.first_order_boundary
                push!(boundary_rows, merge((xi=key,), row))
            end
            if result.cep.found && isfinite(result.cep.T_cep_MeV) && isfinite(result.cep.mu_cep_MeV)
                muq_CEP_MeV = result.cep.mu_cep_MeV
                push!(cep_rows, (
                    xi=key,
                    T_CEP_MeV=result.cep.T_cep_MeV,
                    muq_CEP_MeV=muq_CEP_MeV,
                    muB_CEP_MeV=3.0 * muq_CEP_MeV,
                    uncertainty_T_MeV=result.cep.uncertainty_T_MeV,
                    T_bracket_low_MeV=result.cep.T_bracket_low_MeV,
                    T_bracket_high_MeV=result.cep.T_bracket_high_MeV,
                    bracket_width_T_MeV=result.cep.bracket_width_T_MeV,
                ))
            end
            for row in result.spinodal
                push!(spinodal_rows, merge((xi=key,), row))
            end
            for row in result.crossover_line
                push!(crossover_rows, merge((xi=key,), row))
            end
            append!(
                grid_convergence_rows,
                _dense_record_with_xi.(
                    get(result.diagnostics, "grid_convergence_records", NamedTuple[]),
                    key,
                ),
            )
            push!(manifest["runs"], Dict(
                "xi" => key,
                "run_id" => result.run_id,
                "run_dir" => _repo_relpath(run_dir),
                "boundary_count" => length(result.first_order_boundary),
                "spinodal_count" => length(result.spinodal),
                "crossover_count" => length(result.crossover_line),
                "cep_found" => result.cep.found,
                "cep_uncertainty_T_MeV" => (isfinite(result.cep.uncertainty_T_MeV) ? result.cep.uncertainty_T_MeV : nothing),
                "grid_convergence_count" => length(get(result.diagnostics, "grid_convergence_records", NamedTuple[])),
            ))
        end
    end

    manifest["config"]["requested_xi_values"] = requested_xi_values
    manifest["config"]["xi_values"] = cfg.xi_values
    manifest["grid_convergence"] = Dict(
        "record_count" => length(grid_convergence_rows),
        "unconverged_count" => count(row -> !_dense_record_value(row, :converged, false), grid_convergence_rows),
    )

    sort!(boundary_rows, by=row -> (row.xi, row.T_MeV))
    sort!(cep_rows, by=row -> row.xi)
    sort!(spinodal_rows, by=row -> (row.xi, row.T_MeV))
    sort!(crossover_rows, by=row -> (row.xi, row.mu_MeV))

    boundary_path = joinpath(reference_root, "boundary_$(cfg.tag).csv")
    cep_path = joinpath(reference_root, "cep_$(cfg.tag).csv")
    spinodal_path = joinpath(reference_root, "spinodals_$(cfg.tag).csv")
    crossover_path = joinpath(reference_root, "crossover_$(cfg.tag).csv")
    grid_convergence_path = joinpath(reference_root, "phase_grid_convergence_$(cfg.tag).csv")
    crossover_meta_path = joinpath(reference_root, "crossover_$(cfg.tag).meta.json")
    manifest_path = joinpath(reference_root, "phase_reference_$(cfg.tag)_manifest.json")

    output_paths = cfg.crossover_only ?
        (crossover_path, grid_convergence_path, crossover_meta_path, manifest_path) :
        (boundary_path, cep_path, spinodal_path, crossover_path, grid_convergence_path, crossover_meta_path, manifest_path)
    for path in output_paths
        ensure_writable(path, cfg.overwrite)
    end

    if !cfg.crossover_only
        write_boundary_csv(boundary_path, boundary_rows)
        write_cep_csv(cep_path, cep_rows)
        write_spinodal_csv(spinodal_path, spinodal_rows)
    end
    write_crossover_csv(crossover_path, crossover_rows)
    write_grid_convergence_csv(grid_convergence_path, grid_convergence_rows)
    write_crossover_meta(crossover_meta_path, cfg, crossover_rows, crossover_path)
    manifest["artifacts"]["crossover"] = Dict(
        "path" => _repo_relpath(crossover_path),
        "row_count" => length(crossover_rows),
    )
    manifest["artifacts"]["crossover_meta"] = Dict(
        "path" => _repo_relpath(crossover_meta_path),
    )
    manifest["artifacts"]["grid_convergence"] = Dict(
        "path" => _repo_relpath(grid_convergence_path),
        "row_count" => length(grid_convergence_rows),
    )
    open(manifest_path, "w") do io
        if !cfg.crossover_only
            manifest["artifacts"]["boundary"] = Dict(
                "path" => _repo_relpath(boundary_path),
                "row_count" => length(boundary_rows),
            )
            manifest["artifacts"]["cep"] = Dict(
                "path" => _repo_relpath(cep_path),
                "row_count" => length(cep_rows),
            )
            manifest["artifacts"]["spinodals"] = Dict(
                "path" => _repo_relpath(spinodal_path),
                "row_count" => length(spinodal_rows),
            )
        end
        manifest["artifacts"]["manifest"] = Dict(
            "path" => _repo_relpath(manifest_path),
        )
        JSON3.pretty(io, manifest)
    end

    println("dense reference written:")
    if !cfg.crossover_only
        println("  boundary   = $(boundary_path)")
        println("  cep        = $(cep_path)")
        println("  spinodals  = $(spinodal_path)")
    end
    println("  crossover  = $(crossover_path)")
    println("  grid-convergence = $(grid_convergence_path)")
    println("  crossover-meta = $(crossover_meta_path)")
    println("  manifest   = $(manifest_path)")
end

function main(args=ARGS)
    cfg = parse_args(args)
    build_outputs(cfg)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
