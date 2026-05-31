#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using Printf
using JSON3

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

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
    crossover_only::Bool = false
    crossover_mu_only_zero::Bool = false
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
    println("  --iterations <int>       solver iteration cap (default 80)")
    println("  --tag <name>             output suffix, writes boundary_<tag>.csv etc. default dense")
    println("  --output-root <path>     processed run root")
    println("  --reference-root <path>  reference output directory")
    println("  --overwrite              overwrite existing aggregated outputs")
    println("  --no-crossover           skip crossover generation")
    println("  --crossover-n-mu <int>   crossover mu sampling count (default 16)")
    println("  --crossover-mu-max <MeV> crossover mu_q upper bound (default 450)")
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
        xi_step > 0 || error("xi-step must be positive")
        cfg.xi_values = collect(range(xi_min; stop=xi_max, step=xi_step))
    end

    cfg.T_step > 0 || error("T-step must be positive")
    cfg.rho_step > 0 || error("rho-step must be positive")
    cfg.p_num > 0 || error("p-num must be positive")
    cfg.t_num > 0 || error("t-num must be positive")
    cfg.iterations > 0 || error("iterations must be positive")
    cfg.crossover_n_mu > 0 || error("crossover-n-mu must be positive")
    cfg.crossover_mu_max_MeV > 0 || error("crossover-mu-max must be positive")
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
    return replace(@sprintf("%.3f", xi), "." => "p", "-" => "m")
end

function ensure_writable(path::String, overwrite::Bool)
    if isfile(path) && !overwrite
        error("output exists: $path; rerun with --overwrite")
    end
end

function write_boundary_csv(path::String, rows)
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark,curve_parameter,plot_order_key")
        for row in rows
            println(io, "$(row.xi),$(row.T_MeV),$(row.mu_transition_MeV),$(row.rho_hadron),$(row.rho_quark),$(row.T_MeV),$(row.T_MeV)")
        end
    end
end

function write_cep_csv(path::String, rows)
    open(path, "w") do io
        println(io, "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV")
        for row in rows
            println(io, "$(row.xi),$(row.T_CEP_MeV),$(row.muq_CEP_MeV),$(row.muB_CEP_MeV)")
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
            "max" => min(cfg.T_max, 220.0),
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
        "iterations" => cfg.iterations,
        "mode" => String(cfg.mode),
        "compute_crossover" => cfg.compute_crossover,
        "crossover_method" => String(cfg.crossover_method),
        "crossover_variable" => String(cfg.crossover_variable),
        "crossover_n_mu" => cfg.crossover_n_mu,
        "crossover_mu_max_MeV" => cfg.crossover_mu_max_MeV,
        "crossover_only" => cfg.crossover_only,
        "crossover_mu0_only" => cfg.crossover_mu_only_zero,
    )
end

function build_crossover_only_rows(cfg::DensePhaseReferenceConfig, xi::Float64)
    rows = NamedTuple[]
    if cfg.crossover_mu_only_zero
        T_min_fm = cfg.T_min / 197.327
        T_max_fm = min(cfg.T_max, 220.0) / 197.327
        result = Models.detect_crossover(
            0.0,
            (T_min_fm, T_max_fm);
            method=cfg.crossover_method,
            variable=cfg.crossover_variable,
            xi=xi,
            model_kind=cfg.model_kind,
            solver_backend=cfg.solver_backend,
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
        T_max_MeV=min(cfg.T_max, 220.0),
        xi=xi,
        n_mu=cfg.crossover_n_mu,
        method=cfg.crossover_method,
        variable=cfg.crossover_variable,
        model_kind=cfg.model_kind,
        solver_backend=cfg.solver_backend,
    )
    for row in local_rows
        push!(rows, merge((xi=xi,), row))
    end
    return rows
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
    manifest = Dict{String, Any}(
        "schema_version" => "v1",
        "generator" => manifest_generator_payload(cfg),
        "config" => manifest_config_payload(cfg),
        "output_root" => _repo_relpath(output_root),
        "reference_root" => _repo_relpath(reference_root),
        "artifacts" => Dict{String, Any}(),
        "runs" => Any[],
    )

    T_grid = collect(cfg.T_min:cfg.T_step:cfg.T_max)
    rho_grid = collect(cfg.rho_min:cfg.rho_step:cfg.rho_max)

    for xi in cfg.xi_values
        println("phase reference xi=$(xi): start T=$(cfg.T_min):$(cfg.T_step):$(cfg.T_max), rho=$(cfg.rho_min):$(cfg.rho_step):$(cfg.rho_max)")
        flush(stdout)
        run_dir = joinpath(output_root, "xi_$(xi_token(xi))")
        if cfg.crossover_only
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
        else
            result = Models.run_phase_pipeline(
                cfg.model_kind;
                mode=cfg.mode,
                T_grid=T_grid,
                rho_grid=rho_grid,
                xi=xi,
                output_dir=run_dir,
                profile=cfg.profile,
                solver_backend=cfg.solver_backend,
                seed_policy=cfg.seed_policy,
                p_num=cfg.p_num,
                t_num=cfg.t_num,
                iterations=cfg.iterations,
                compute_crossover=cfg.compute_crossover,
                crossover_method=cfg.crossover_method,
                crossover_variable=cfg.crossover_variable,
                crossover_n_mu=cfg.crossover_n_mu,
                promote_reference=false,
            )

            for row in result.first_order_boundary
                push!(boundary_rows, merge((xi=xi,), row))
            end
            if result.cep.found && isfinite(result.cep.T_cep_MeV) && isfinite(result.cep.mu_cep_MeV)
                muq_CEP_MeV = result.cep.mu_cep_MeV
                push!(cep_rows, (xi=xi, T_CEP_MeV=result.cep.T_cep_MeV, muq_CEP_MeV=muq_CEP_MeV, muB_CEP_MeV=3.0 * muq_CEP_MeV))
            end
            for row in result.spinodal
                push!(spinodal_rows, merge((xi=xi,), row))
            end
            for row in result.crossover_line
                push!(crossover_rows, merge((xi=xi,), row))
            end

            println("phase reference xi=$(xi): boundary=$(length(result.first_order_boundary)) spinodal=$(length(result.spinodal)) crossover=$(length(result.crossover_line)) cep_found=$(result.cep.found)")
            flush(stdout)
            push!(manifest["runs"], Dict(
                "xi" => xi,
                "run_id" => result.run_id,
                "run_dir" => _repo_relpath(run_dir),
                "boundary_count" => length(result.first_order_boundary),
                "spinodal_count" => length(result.spinodal),
                "crossover_count" => length(result.crossover_line),
                "cep_found" => result.cep.found,
            ))
        end
    end

    sort!(boundary_rows, by=row -> (row.xi, row.T_MeV))
    sort!(cep_rows, by=row -> row.xi)
    sort!(spinodal_rows, by=row -> (row.xi, row.T_MeV))
    sort!(crossover_rows, by=row -> (row.xi, row.mu_MeV))

    boundary_path = joinpath(reference_root, "boundary_$(cfg.tag).csv")
    cep_path = joinpath(reference_root, "cep_$(cfg.tag).csv")
    spinodal_path = joinpath(reference_root, "spinodals_$(cfg.tag).csv")
    crossover_path = joinpath(reference_root, "crossover_$(cfg.tag).csv")
    crossover_meta_path = joinpath(reference_root, "crossover_$(cfg.tag).meta.json")
    manifest_path = joinpath(reference_root, "phase_reference_$(cfg.tag)_manifest.json")

    output_paths = cfg.crossover_only ? (crossover_path, crossover_meta_path, manifest_path) : (boundary_path, cep_path, spinodal_path, crossover_path, crossover_meta_path, manifest_path)
    for path in output_paths
        ensure_writable(path, cfg.overwrite)
    end

    if !cfg.crossover_only
        write_boundary_csv(boundary_path, boundary_rows)
        write_cep_csv(cep_path, cep_rows)
        write_spinodal_csv(spinodal_path, spinodal_rows)
    end
    write_crossover_csv(crossover_path, crossover_rows)
    write_crossover_meta(crossover_meta_path, cfg, crossover_rows, crossover_path)
    manifest["artifacts"]["crossover"] = Dict(
        "path" => _repo_relpath(crossover_path),
        "row_count" => length(crossover_rows),
    )
    manifest["artifacts"]["crossover_meta"] = Dict(
        "path" => _repo_relpath(crossover_meta_path),
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
