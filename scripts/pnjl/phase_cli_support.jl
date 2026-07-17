module PhaseCliSupport

using Dates
using TOML
using JSON3

export PhaseCliConfig, parse_args, main

const DEFAULT_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

Base.@kwdef mutable struct PhaseCliConfig
    config_path::Union{Nothing, String} = nothing
    preset::Union{Nothing, Symbol} = nothing
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
    thermo_quadrature_policy::Symbol = :tensor_gauss
    thermo_quadrature_rtol::Float64 = 1e-8
    thermo_quadrature_atol::Float64 = 1e-10
    thermo_quadrature_maxevals::Int = 10^7
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

@inline function _norm_slash(path::AbstractString)
    return replace(String(path), '\\' => '/')
end

function _repo_relpath(path::AbstractString, project_root::AbstractString)
    abs_path = normpath(abspath(String(path)))
    rel = try
        relpath(abs_path, project_root)
    catch
        nothing
    end
    if rel !== nothing
        return _norm_slash(String(rel))
    end
    return _norm_slash(abs_path)
end

function _normalize_artifact_paths(paths, project_root::AbstractString)
    normalized = Dict{String, Any}()
    for (key, value) in pairs(paths)
        key_str = String(key)
        if value isa AbstractString
            normalized[key_str] = _repo_relpath(String(value), project_root)
        else
            normalized[key_str] = value
        end
    end
    return normalized
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

function _default_phase_config_path(model_kind::Symbol, project_root::AbstractString)::Union{Nothing, String}
    model_tag = lowercase(String(model_kind))
    path = joinpath(project_root, "config", "models", model_tag, "phase_pipeline_default.toml")
    return isfile(path) ? path : nothing
end

function _apply_preset!(cfg::PhaseCliConfig, preset::Symbol)
    if preset == :smoke
        cfg.mode = :research
        cfg.profile = :smoke
        cfg.solver_backend = :models
        cfg.T_min = 150.0
        cfg.T_max = 150.0
        cfg.T_step = 10.0
        cfg.rho_min = 0.1
        cfg.rho_max = 0.3
        cfg.rho_step = 0.1
        cfg.p_num = 12
        cfg.t_num = 4
        cfg.iterations = 10
        return cfg
    end
    throw(ArgumentError("invalid --preset=$(preset); accepted values: smoke"))
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
    haskey(table, "thermo_quadrature_policy") && (cfg.thermo_quadrature_policy = Symbol(lowercase(String(table["thermo_quadrature_policy"]))))
    haskey(table, "thermo_quadrature_rtol") && (cfg.thermo_quadrature_rtol = Float64(table["thermo_quadrature_rtol"]))
    haskey(table, "thermo_quadrature_atol") && (cfg.thermo_quadrature_atol = Float64(table["thermo_quadrature_atol"]))
    haskey(table, "thermo_quadrature_maxevals") && (cfg.thermo_quadrature_maxevals = Int(table["thermo_quadrature_maxevals"]))
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

function _write_run_manifest(output_dir::String, cfg::PhaseCliConfig, args::Vector{String}, result, project_root::AbstractString)
    manifest_path = joinpath(output_dir, "run_manifest.json")
    existing_manifest = if isfile(manifest_path)
        try
            JSON3.read(read(manifest_path, String))
        catch err
            @warn "Failed to parse existing run_manifest.json; writing fresh CLI manifest projection" manifest_path exception=(err, catch_backtrace())
            nothing
        end
    else
        nothing
    end

    existing_pipeline_payload = if existing_manifest !== nothing && haskey(existing_manifest, :pipeline)
        pipe = existing_manifest[:pipeline]
        same_model_kind = haskey(pipe, :model_kind) && String(pipe[:model_kind]) == String(cfg.model_kind)
        if same_model_kind
            Dict(
                "name" => get(pipe, :name, "phase_pipeline_runner"),
                "version" => get(pipe, :version, "v1"),
                "model_kind" => String(cfg.model_kind),
                "run_id" => String(result.run_id),
                "git_commit" => get(pipe, :git_commit, nothing),
                "manifest_schema_version" => get(pipe, :manifest_schema_version, "phase_manifest_v1"),
                "timestamp" => get(pipe, :timestamp, nothing),
                "success" => get(pipe, :success, true),
                "failed_stage" => get(pipe, :failed_stage, nothing),
                "error_kind" => get(pipe, :error_kind, nothing),
                "error_msg" => get(pipe, :error_msg, nothing),
                "config_hash" => get(result.config_snapshot, "config_hash", nothing),
                "artifact_hash" => get(pipe, :artifact_hash, nothing),
            )
        else
            nothing
        end
    else
        nothing
    end

    git_commit = try
        readchomp(`git -C $(project_root) rev-parse HEAD`)
    catch
        nothing
    end
    effective_config = Dict(
        "model_kind" => String(cfg.model_kind),
        "mode" => String(cfg.mode),
        "profile" => String(cfg.profile),
        "xi" => cfg.xi,
        "T_min" => cfg.T_min,
        "T_max" => cfg.T_max,
        "T_step" => cfg.T_step,
        "rho_min" => cfg.rho_min,
        "rho_max" => cfg.rho_max,
        "rho_step" => cfg.rho_step,
        "solver_backend" => String(cfg.solver_backend),
        "seed_policy" => String(cfg.seed_policy),
        "reverse_rho" => cfg.reverse_rho,
        "p_num" => cfg.p_num,
        "t_num" => cfg.t_num,
        "thermo_quadrature_policy" => String(cfg.thermo_quadrature_policy),
        "thermo_quadrature_rtol" => cfg.thermo_quadrature_rtol,
        "thermo_quadrature_atol" => cfg.thermo_quadrature_atol,
        "thermo_quadrature_maxevals" => cfg.thermo_quadrature_maxevals,
        "iterations" => cfg.iterations,
        "compute_crossover" => cfg.compute_crossover,
        "promote_reference" => cfg.promote_reference,
    )

    payload = Dict(
        "generated_at" => string(now()),
        "git_commit" => git_commit,
        "argv" => collect(String.(args)),
        "config_path" => isnothing(cfg.config_path) ? nothing : _repo_relpath(cfg.config_path, project_root),
        "preset" => isnothing(cfg.preset) ? nothing : String(cfg.preset),
        "config_hash" => get(result.config_snapshot, "config_hash", nothing),
        "run_id" => result.run_id,
        "mode" => String(cfg.mode),
        "model_kind" => String(cfg.model_kind),
        "artifact_paths" => _normalize_artifact_paths(result.artifact_paths, project_root),
        "effective_config" => effective_config,
    )

    if existing_pipeline_payload !== nothing
        payload["pipeline"] = existing_pipeline_payload
        haskey(existing_manifest, :completed_stages) && (payload["completed_stages"] = existing_manifest[:completed_stages])
        haskey(existing_manifest, :stage_records) && (payload["stage_records"] = existing_manifest[:stage_records])
    end

    open(manifest_path, "w") do io
        write(io, JSON3.write(payload))
    end
    return manifest_path
end

function _usage()
    println("用法: julia scripts/pnjl/calculate_phase_structure.jl [options]")
    println("选项:")
    println("  --config=...           指定 phase pipeline TOML 配置文件")
    println("  --preset=smoke         使用内置快速复现参数模板（可被后续CLI参数覆盖）")
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
    println("  --solver_backend=models|auto")
    println("  --seed_policy=...      扫描初值策略")
    println("  --reverse_rho=true|false")
    println("  --p_num=12             动量积分点数")
    println("  --t_num=4              角度积分点数")
    println("  --thermo_quadrature_policy=tensor_gauss|rs_reduced_adaptive")
    println("  --thermo_quadrature_rtol=1e-8")
    println("  --thermo_quadrature_atol=1e-10")
    println("  --thermo_quadrature_maxevals=10000000")
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

function parse_args(args, project_root::AbstractString)
    cfg = PhaseCliConfig()

    explicit_config = nothing
    preset = nothing
    for arg in args
        if startswith(arg, "--config=")
            explicit_config = arg[10:end]
        elseif startswith(arg, "--preset=")
            preset = Symbol(lowercase(split(arg, "="; limit=2)[2]))
        end
    end

    if explicit_config !== nothing
        cfg.config_path = explicit_config
        _apply_phase_config!(cfg, _phase_config_table(cfg.config_path))
    else
        model_hint = _extract_model_kind_hint(args)
        default_cfg = _default_phase_config_path(model_hint, project_root)
        if default_cfg !== nothing
            cfg.config_path = default_cfg
            _apply_phase_config!(cfg, _phase_config_table(default_cfg))
        end
    end

    if preset !== nothing
        cfg.preset = preset
        _apply_preset!(cfg, preset)
    end

    for arg in args
        if arg in ("-h", "--help")
            _usage()
            return nothing
        elseif startswith(arg, "--config=")
            continue
        elseif startswith(arg, "--preset=")
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
        elseif startswith(arg, "--thermo_quadrature_policy=")
            cfg.thermo_quadrature_policy = Symbol(lowercase(split(arg, "="; limit=2)[2]))
        elseif startswith(arg, "--thermo_quadrature_rtol=")
            cfg.thermo_quadrature_rtol = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--thermo_quadrature_atol=")
            cfg.thermo_quadrature_atol = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--thermo_quadrature_maxevals=")
            cfg.thermo_quadrature_maxevals = parse(Int, split(arg, "="; limit=2)[2])
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
    cfg.solver_backend in (:models, :auto) || throw(ArgumentError("invalid --solver_backend=$(cfg.solver_backend); accepted values: models, auto"))
    cfg.thermo_quadrature_policy in (:tensor_gauss, :rs_reduced_adaptive) || throw(ArgumentError(
        "unsupported thermo_quadrature_policy=$(cfg.thermo_quadrature_policy)",
    ))
    if cfg.thermo_quadrature_policy === :rs_reduced_adaptive && cfg.model_kind !== :PNJL
        throw(ArgumentError("rs_reduced_adaptive thermal quadrature is supported only for model_kind=PNJL"))
    end
    isfinite(cfg.thermo_quadrature_rtol) && cfg.thermo_quadrature_rtol > 0 || throw(ArgumentError("thermo_quadrature_rtol must be finite and positive"))
    isfinite(cfg.thermo_quadrature_atol) && cfg.thermo_quadrature_atol >= 0 || throw(ArgumentError("thermo_quadrature_atol must be finite and nonnegative"))
    cfg.thermo_quadrature_maxevals > 0 || throw(ArgumentError("thermo_quadrature_maxevals must be positive"))
    return cfg
end

parse_args(args) = parse_args(args, DEFAULT_PROJECT_ROOT)

function main(models_module, project_root::AbstractString, args::Vector{String}=collect(String.(ARGS)))
    cfg = parse_args(args, project_root)
    cfg === nothing && return 0

    T_grid = collect(cfg.T_min:cfg.T_step:cfg.T_max)
    rho_grid = collect(cfg.rho_min:cfg.rho_step:cfg.rho_max)

    println("="^60)
    println("Phase pipeline CLI")
    println("="^60)
    println("time=$(now()) model_kind=$(cfg.model_kind) mode=$(cfg.mode) profile=$(cfg.profile)")
    println("T-grid: $(first(T_grid)) -> $(last(T_grid)) (n=$(length(T_grid)))")
    println("rho-grid: $(first(rho_grid)) -> $(last(rho_grid)) (n=$(length(rho_grid)))")

    result = models_module.run_phase_pipeline(
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
        thermo_quadrature_policy=cfg.thermo_quadrature_policy,
        thermo_quadrature_rtol=cfg.thermo_quadrature_rtol,
        thermo_quadrature_atol=cfg.thermo_quadrature_atol,
        thermo_quadrature_maxevals=cfg.thermo_quadrature_maxevals,
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
    manifest_path = _write_run_manifest(output_root, cfg, collect(String.(args)), result, project_root)

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
    return 0
end

end # module PhaseCliSupport
