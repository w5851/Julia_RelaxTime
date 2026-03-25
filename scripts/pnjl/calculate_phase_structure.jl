#!/usr/bin/env julia
"""
PNJL 相结构计算脚本（薄 CLI）

职责：
- 仅做参数解析、调用 Models.run_phase_pipeline
- 可选触发晋升到 data/reference

示例：
    julia scripts/pnjl/calculate_phase_structure.jl --model_kind=PNJL --verbose
    julia scripts/pnjl/calculate_phase_structure.jl --T_min=120 --T_max=180 --T_step=20 --rho_max=0.8 --rho_step=0.2
    julia scripts/pnjl/calculate_phase_structure.jl --promote_reference
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using TOML
using JSON3

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

Base.@kwdef mutable struct PhaseCliConfig
    config_path::Union{Nothing, String} = nothing
    model_kind::Symbol = :PNJL
    mode::Symbol = :production
    xi::Float64 = 0.0
    T_min::Float64 = 30.0
    T_max::Float64 = 350.0
    T_step::Float64 = 10.0
    rho_min::Float64 = 0.0
    rho_max::Float64 = 4.0
    rho_step::Float64 = 0.05
    profile::Symbol = :default
    run_id::Union{Nothing, String} = nothing
    output_dir::Union{Nothing, String} = nothing
    solver_backend::Symbol = :models
    seed_policy::Symbol = :hybrid_continuity
    reverse_rho::Bool = true
    p_num::Int = 12
    t_num::Int = 4
    iterations::Int = 80
    compute_crossover::Bool = false
    crossover_method::Symbol = :peak
    crossover_variable::Symbol = :phi_u
    crossover_n_mu::Int = 12
    cep_strategy::Symbol = :interpolate
    cep_interpolate_use_direct_eval::Bool = false
    cep_tol::Float64 = 0.01
    cep_max_bisect_iter::Int = 20
    cep_area_tol_good::Float64 = 1e-4
    cep_area_tol_bad::Float64 = 5e-4
    cep_max_refine_level::Int = 2
    cep_adaptive_rho::Bool = true
    cep_adaptive_slope_tol::Float64 = 5.0
    cep_adaptive_min_gap::Float64 = 0.002
    cep_adaptive_max_points::Int = 32
    cep_adaptive_digits::Int = 6
    cep_direct_bracket_mode::Symbol = :directional
    cep_direct_start::Symbol = :low
    cep_direct_initial_step::Float64 = NaN
    cep_direct_expand_factor::Float64 = 2.0
    cep_direct_max_expand_steps::Int = 8
    cep_direct_fallback_scan::Bool = true
    promote_reference::Bool = false
    verbose::Bool = false
end

function _as_bool(x, name::AbstractString)
    if x isa Bool
        return x
    elseif x isa Integer
        return x != 0
    elseif x isa AbstractString
        token = lowercase(strip(x))
        token in ("1", "true", "yes", "on") && return true
        token in ("0", "false", "no", "off") && return false
    end
    throw(ArgumentError("invalid boolean value for $(name): $(repr(x))"))
end

function _phase_config_table(path::String)
    isfile(path) || throw(ArgumentError("config file does not exist: $(path)"))
    raw = TOML.parsefile(path)
    haskey(raw, "phase_pipeline") || throw(ArgumentError("missing [phase_pipeline] section in config: $(path)"))
    table = raw["phase_pipeline"]
    table isa AbstractDict || throw(ArgumentError("[phase_pipeline] must be a TOML table in config: $(path)"))
    return table
end

function _extract_model_kind_hint(args)::Symbol
    for arg in args
        startswith(arg, "--model_kind=") || continue
        return Symbol(uppercase(arg[14:end]))
    end
    return :PNJL
end

function _default_phase_config_path(model_kind::Symbol)::Union{Nothing, String}
    model_tag = lowercase(String(model_kind))
    path = joinpath(PROJECT_ROOT, "config", "models", model_tag, "phase_pipeline_default.toml")
    return isfile(path) ? path : nothing
end

function _apply_phase_config!(cfg::PhaseCliConfig, table::AbstractDict)
    haskey(table, "model_kind") && (cfg.model_kind = Symbol(uppercase(String(table["model_kind"]))))
    haskey(table, "mode") && (cfg.mode = Symbol(lowercase(String(table["mode"]))))
    haskey(table, "xi") && (cfg.xi = Float64(table["xi"]))
    haskey(table, "T_min") && (cfg.T_min = Float64(table["T_min"]))
    haskey(table, "T_max") && (cfg.T_max = Float64(table["T_max"]))
    haskey(table, "T_step") && (cfg.T_step = Float64(table["T_step"]))
    haskey(table, "rho_min") && (cfg.rho_min = Float64(table["rho_min"]))
    haskey(table, "rho_max") && (cfg.rho_max = Float64(table["rho_max"]))
    haskey(table, "rho_step") && (cfg.rho_step = Float64(table["rho_step"]))
    haskey(table, "profile") && (cfg.profile = Symbol(lowercase(String(table["profile"]))))
    haskey(table, "run_id") && (cfg.run_id = String(table["run_id"]))
    haskey(table, "output_dir") && (cfg.output_dir = String(table["output_dir"]))
    haskey(table, "solver_backend") && (cfg.solver_backend = Symbol(lowercase(String(table["solver_backend"]))))
    haskey(table, "seed_policy") && (cfg.seed_policy = Symbol(lowercase(String(table["seed_policy"]))))
    haskey(table, "reverse_rho") && (cfg.reverse_rho = _as_bool(table["reverse_rho"], "phase_pipeline.reverse_rho"))
    haskey(table, "p_num") && (cfg.p_num = Int(table["p_num"]))
    haskey(table, "t_num") && (cfg.t_num = Int(table["t_num"]))
    haskey(table, "iterations") && (cfg.iterations = Int(table["iterations"]))
    haskey(table, "compute_crossover") && (cfg.compute_crossover = _as_bool(table["compute_crossover"], "phase_pipeline.compute_crossover"))
    haskey(table, "crossover_method") && (cfg.crossover_method = Symbol(lowercase(String(table["crossover_method"]))))
    haskey(table, "crossover_variable") && (cfg.crossover_variable = Symbol(String(table["crossover_variable"])))
    haskey(table, "crossover_n_mu") && (cfg.crossover_n_mu = Int(table["crossover_n_mu"]))
    haskey(table, "cep_strategy") && (cfg.cep_strategy = Symbol(lowercase(String(table["cep_strategy"]))))
    haskey(table, "cep_interpolate_use_direct_eval") && (cfg.cep_interpolate_use_direct_eval = _as_bool(table["cep_interpolate_use_direct_eval"], "phase_pipeline.cep_interpolate_use_direct_eval"))
    haskey(table, "cep_tol") && (cfg.cep_tol = Float64(table["cep_tol"]))
    haskey(table, "cep_max_bisect_iter") && (cfg.cep_max_bisect_iter = Int(table["cep_max_bisect_iter"]))
    haskey(table, "cep_area_tol_good") && (cfg.cep_area_tol_good = Float64(table["cep_area_tol_good"]))
    haskey(table, "cep_area_tol_bad") && (cfg.cep_area_tol_bad = Float64(table["cep_area_tol_bad"]))
    haskey(table, "cep_max_refine_level") && (cfg.cep_max_refine_level = Int(table["cep_max_refine_level"]))
    haskey(table, "cep_adaptive_rho") && (cfg.cep_adaptive_rho = _as_bool(table["cep_adaptive_rho"], "phase_pipeline.cep_adaptive_rho"))
    haskey(table, "cep_adaptive_slope_tol") && (cfg.cep_adaptive_slope_tol = Float64(table["cep_adaptive_slope_tol"]))
    haskey(table, "cep_adaptive_min_gap") && (cfg.cep_adaptive_min_gap = Float64(table["cep_adaptive_min_gap"]))
    haskey(table, "cep_adaptive_max_points") && (cfg.cep_adaptive_max_points = Int(table["cep_adaptive_max_points"]))
    haskey(table, "cep_adaptive_digits") && (cfg.cep_adaptive_digits = Int(table["cep_adaptive_digits"]))
    haskey(table, "cep_direct_bracket_mode") && (cfg.cep_direct_bracket_mode = Symbol(lowercase(String(table["cep_direct_bracket_mode"]))))
    haskey(table, "cep_direct_start") && (cfg.cep_direct_start = Symbol(lowercase(String(table["cep_direct_start"]))))
    haskey(table, "cep_direct_initial_step") && (cfg.cep_direct_initial_step = Float64(table["cep_direct_initial_step"]))
    haskey(table, "cep_direct_expand_factor") && (cfg.cep_direct_expand_factor = Float64(table["cep_direct_expand_factor"]))
    haskey(table, "cep_direct_max_expand_steps") && (cfg.cep_direct_max_expand_steps = Int(table["cep_direct_max_expand_steps"]))
    haskey(table, "cep_direct_fallback_scan") && (cfg.cep_direct_fallback_scan = _as_bool(table["cep_direct_fallback_scan"], "phase_pipeline.cep_direct_fallback_scan"))
    haskey(table, "promote_reference") && (cfg.promote_reference = _as_bool(table["promote_reference"], "phase_pipeline.promote_reference"))
    haskey(table, "verbose") && (cfg.verbose = _as_bool(table["verbose"], "phase_pipeline.verbose"))
    return cfg
end

function _write_run_manifest(output_dir::String, cfg::PhaseCliConfig, args::Vector{String}, result)
    manifest_path = joinpath(output_dir, "run_manifest.json")
    git_commit = try
        readchomp(`git -C $(joinpath(@__DIR__, "..", "..")) rev-parse HEAD`)
    catch
        nothing
    end
    payload = Dict(
        "generated_at" => string(now()),
        "git_commit" => git_commit,
        "argv" => collect(String.(args)),
        "config_path" => cfg.config_path,
        "config_hash" => get(result.config_snapshot, "config_hash", nothing),
        "run_id" => result.run_id,
        "mode" => String(cfg.mode),
        "model_kind" => String(cfg.model_kind),
        "artifact_paths" => result.artifact_paths,
    )
    open(manifest_path, "w") do io
        write(io, JSON3.write(payload))
    end
    return manifest_path
end

function _usage()
    println("用法: julia scripts/pnjl/calculate_phase_structure.jl [options]")
    println("选项:")
    println("  --config=...           指定 phase pipeline TOML 配置文件")
    println("  --model_kind=PNJL      模型类型（如 PNJL/RPNJL）")
    println("  --mode=production|research  运行模式")
    println("  --xi=0.0               各向异性参数")
    println("  --T_min=30             最低温度 (MeV)")
    println("  --T_max=350            最高温度 (MeV)")
    println("  --T_step=10            温度步长 (MeV)")
    println("  --rho_min=0.0          最低密度 (ρ/ρ₀)")
    println("  --rho_max=4.0          最高密度 (ρ/ρ₀)")
    println("  --rho_step=0.05        密度步长")
    println("  --profile=default      运行 profile")
    println("  --run_id=...           指定 run_id")
    println("  --output_dir=...       指定输出目录（覆盖默认 processed 路径）")
    println("  --solver_backend=models|legacy|auto")
    println("  --seed_policy=...      扫描初值策略")
    println("  --reverse_rho=true|false")
    println("  --p_num=12             动量积分点数")
    println("  --t_num=4              角度积分点数")
    println("  --iterations=80        求解迭代上限")
    println("  --compute_crossover    计算并输出 crossover_line.csv")
    println("  --crossover_method=... crossover方法（inflection|peak）")
    println("  --crossover_variable=... crossover变量（phi_u|Phi）")
    println("  --crossover_n_mu=12    crossover扫描的mu点数")
    println("  --cep_strategy=...     CEP定位策略（interpolate|direct）")
    println("  --cep_interpolate_use_direct_eval=true|false interpolate策略下对临界二分点做direct重算")
    println("  --cep_tol=0.01         CEP二分温度容差 (MeV)")
    println("  --cep_max_bisect_iter=20 CEP二分迭代上限")
    println("  --cep_area_tol_good=1e-4 CEP判定valid阈值")
    println("  --cep_area_tol_bad=5e-4  CEP判定invalid阈值")
    println("  --cep_max_refine_level=2 CEP直接二分每点加密上限")
    println("  --cep_adaptive_rho=true  direct评估时启用rho自适应补点")
    println("  --cep_adaptive_slope_tol=5.0 自适应补点斜率阈值")
    println("  --cep_adaptive_min_gap=0.002 自适应最小区间")
    println("  --cep_adaptive_max_points=32 自适应单轮最多补点")
    println("  --cep_adaptive_digits=6 自适应去重精度")
    println("  --cep_direct_bracket_mode=... direct bracket模式（directional|scan）")
    println("  --cep_direct_start=...   directional起点（low|mid|high）")
    println("  --cep_direct_initial_step=... directional初始步长(MeV, NaN=自动)")
    println("  --cep_direct_expand_factor=2.0 directional步长扩张倍数")
    println("  --cep_direct_max_expand_steps=8 directional最大扩张步数")
    println("  --cep_direct_fallback_scan=true directional失败后是否回退扫描")
    println("  --promote_reference    运行后尝试晋升到 data/reference")
    println("  --verbose              输出详细信息")
    println("  -h, --help             显示帮助")
end

function parse_args(args)
    cfg = PhaseCliConfig()

    explicit_config = nothing
    for arg in args
        if startswith(arg, "--config=")
            explicit_config = arg[10:end]
            break
        end
    end

    if explicit_config !== nothing
        cfg.config_path = explicit_config
        _apply_phase_config!(cfg, _phase_config_table(cfg.config_path))
    else
        model_hint = _extract_model_kind_hint(args)
        default_cfg = _default_phase_config_path(model_hint)
        if default_cfg !== nothing
            cfg.config_path = default_cfg
            _apply_phase_config!(cfg, _phase_config_table(default_cfg))
        end
    end

    for arg in args
        if arg in ("-h", "--help")
            _usage()
            exit(0)
        elseif startswith(arg, "--config=")
            continue
        elseif startswith(arg, "--model_kind=")
            cfg.model_kind = Symbol(uppercase(arg[14:end]))
        elseif startswith(arg, "--mode=")
            cfg.mode = Symbol(lowercase(split(arg, "="; limit=2)[2]))
        elseif startswith(arg, "--xi=")
            cfg.xi = parse(Float64, arg[6:end])
        elseif startswith(arg, "--T_min=")
            cfg.T_min = parse(Float64, arg[9:end])
        elseif startswith(arg, "--T_max=")
            cfg.T_max = parse(Float64, arg[9:end])
        elseif startswith(arg, "--T_step=")
            cfg.T_step = parse(Float64, arg[10:end])
        elseif startswith(arg, "--rho_min=")
            cfg.rho_min = parse(Float64, arg[11:end])
        elseif startswith(arg, "--rho_max=")
            cfg.rho_max = parse(Float64, arg[11:end])
        elseif startswith(arg, "--rho_step=")
            cfg.rho_step = parse(Float64, arg[12:end])
        elseif startswith(arg, "--profile=")
            cfg.profile = Symbol(lowercase(arg[11:end]))
        elseif startswith(arg, "--run_id=")
            cfg.run_id = arg[10:end]
        elseif startswith(arg, "--output_dir=")
            cfg.output_dir = arg[14:end]
        elseif startswith(arg, "--solver_backend=")
            cfg.solver_backend = Symbol(lowercase(arg[18:end]))
        elseif startswith(arg, "--seed_policy=")
            cfg.seed_policy = Symbol(lowercase(arg[15:end]))
        elseif startswith(arg, "--reverse_rho=")
            cfg.reverse_rho = lowercase(arg[15:end]) in ("1", "true", "yes")
        elseif startswith(arg, "--p_num=")
            cfg.p_num = parse(Int, arg[9:end])
        elseif startswith(arg, "--t_num=")
            cfg.t_num = parse(Int, arg[9:end])
        elseif startswith(arg, "--iterations=")
            cfg.iterations = parse(Int, arg[14:end])
        elseif arg == "--compute_crossover"
            cfg.compute_crossover = true
        elseif startswith(arg, "--crossover_method=")
            cfg.crossover_method = Symbol(lowercase(arg[20:end]))
        elseif startswith(arg, "--crossover_variable=")
            cfg.crossover_variable = Symbol(arg[22:end])
        elseif startswith(arg, "--crossover_n_mu=")
            cfg.crossover_n_mu = parse(Int, arg[18:end])
        elseif startswith(arg, "--cep_strategy=")
            cfg.cep_strategy = Symbol(lowercase(split(arg, "="; limit=2)[2]))
        elseif startswith(arg, "--cep_interpolate_use_direct_eval=")
            cfg.cep_interpolate_use_direct_eval = lowercase(split(arg, "="; limit=2)[2]) in ("1", "true", "yes")
        elseif startswith(arg, "--cep_tol=")
            cfg.cep_tol = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_max_bisect_iter=")
            cfg.cep_max_bisect_iter = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_area_tol_good=")
            cfg.cep_area_tol_good = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_area_tol_bad=")
            cfg.cep_area_tol_bad = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_max_refine_level=")
            cfg.cep_max_refine_level = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_adaptive_rho=")
            cfg.cep_adaptive_rho = lowercase(split(arg, "="; limit=2)[2]) in ("1", "true", "yes")
        elseif startswith(arg, "--cep_adaptive_slope_tol=")
            cfg.cep_adaptive_slope_tol = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_adaptive_min_gap=")
            cfg.cep_adaptive_min_gap = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_adaptive_max_points=")
            cfg.cep_adaptive_max_points = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_adaptive_digits=")
            cfg.cep_adaptive_digits = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_direct_bracket_mode=")
            cfg.cep_direct_bracket_mode = Symbol(lowercase(split(arg, "="; limit=2)[2]))
        elseif startswith(arg, "--cep_direct_start=")
            cfg.cep_direct_start = Symbol(lowercase(split(arg, "="; limit=2)[2]))
        elseif startswith(arg, "--cep_direct_initial_step=")
            cfg.cep_direct_initial_step = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_direct_expand_factor=")
            cfg.cep_direct_expand_factor = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_direct_max_expand_steps=")
            cfg.cep_direct_max_expand_steps = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cep_direct_fallback_scan=")
            cfg.cep_direct_fallback_scan = lowercase(split(arg, "="; limit=2)[2]) in ("1", "true", "yes")
        elseif arg == "--promote_reference"
            cfg.promote_reference = true
        elseif arg == "--verbose"
            cfg.verbose = true
        else
            error("Unknown option: $arg")
        end
    end
    cfg.mode in (:production, :research) || throw(ArgumentError("invalid --mode=$(cfg.mode); accepted values: production, research"))
    return cfg
end

function main(args=ARGS)
    cfg = parse_args(args)

    T_grid = collect(cfg.T_min:cfg.T_step:cfg.T_max)
    rho_grid = collect(cfg.rho_min:cfg.rho_step:cfg.rho_max)

    println("="^60)
    println("Phase pipeline CLI")
    println("="^60)
    println("time=$(now()) model_kind=$(cfg.model_kind) mode=$(cfg.mode) profile=$(cfg.profile)")
    println("T-grid: $(first(T_grid)) -> $(last(T_grid)) (n=$(length(T_grid)))")
    println("rho-grid: $(first(rho_grid)) -> $(last(rho_grid)) (n=$(length(rho_grid)))")

    result = Models.run_phase_pipeline(
        cfg.model_kind;
        mode=cfg.mode,
        T_grid=T_grid,
        rho_grid=rho_grid,
        xi=cfg.xi,
        output_dir=cfg.output_dir,
        profile=cfg.profile,
        run_id=cfg.run_id,
        solver_backend=cfg.solver_backend,
        seed_policy=cfg.seed_policy,
        reverse_rho=cfg.reverse_rho,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        iterations=cfg.iterations,
        compute_crossover=cfg.compute_crossover,
        crossover_method=cfg.crossover_method,
        crossover_variable=cfg.crossover_variable,
        crossover_n_mu=cfg.crossover_n_mu,
        cep_strategy=cfg.cep_strategy,
        cep_interpolate_use_direct_eval=cfg.cep_interpolate_use_direct_eval,
        cep_tol=cfg.cep_tol,
        cep_max_bisect_iter=cfg.cep_max_bisect_iter,
        cep_area_tol_good=cfg.cep_area_tol_good,
        cep_area_tol_bad=cfg.cep_area_tol_bad,
        cep_max_refine_level=cfg.cep_max_refine_level,
        cep_adaptive_rho=cfg.cep_adaptive_rho,
        cep_adaptive_slope_tol=cfg.cep_adaptive_slope_tol,
        cep_adaptive_min_gap=cfg.cep_adaptive_min_gap,
        cep_adaptive_max_points=cfg.cep_adaptive_max_points,
        cep_adaptive_digits=cfg.cep_adaptive_digits,
        cep_direct_bracket_mode=cfg.cep_direct_bracket_mode,
        cep_direct_start=cfg.cep_direct_start,
        cep_direct_initial_step=cfg.cep_direct_initial_step,
        cep_direct_expand_factor=cfg.cep_direct_expand_factor,
        cep_direct_max_expand_steps=cfg.cep_direct_max_expand_steps,
        cep_direct_fallback_scan=cfg.cep_direct_fallback_scan,
        promote_reference=cfg.promote_reference,
    )

    output_root = isnothing(cfg.output_dir) ? dirname(result.artifact_paths["phase_summary"]) : cfg.output_dir
    manifest_path = _write_run_manifest(output_root, cfg, collect(String.(args)), result)

    println("\n完成:")
    println("  run_id = $(result.run_id)")
    println("  CEP found = $(result.cep.found)")
    println("  boundary_count = $(length(result.first_order_boundary))")
    println("  artifacts = $(result.artifact_paths)")
    println("  manifest = $(manifest_path)")
    if cfg.promote_reference
        println("  promotion = passed=$(result.promotion_status.passed), baseline_id=$(result.promotion_status.baseline_id)")
        if !isempty(result.promotion_status.failed_checks)
            println("  promotion_failed_checks = $(result.promotion_status.failed_checks)")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
