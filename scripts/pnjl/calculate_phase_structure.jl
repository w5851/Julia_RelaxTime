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

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

Base.@kwdef mutable struct PhaseCliConfig
    model_kind::Symbol = :PNJL
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
    crossover_method::Symbol = :inflection
    crossover_variable::Symbol = :phi_u
    crossover_n_mu::Int = 12
    cep_strategy::Symbol = :interpolate
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

function _usage()
    println("用法: julia scripts/pnjl/calculate_phase_structure.jl [options]")
    println("选项:")
    println("  --model_kind=PNJL      模型类型（如 PNJL/RPNJL）")
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
    for arg in args
        if arg in ("-h", "--help")
            _usage()
            exit(0)
        elseif startswith(arg, "--model_kind=")
            cfg.model_kind = Symbol(uppercase(arg[14:end]))
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
    return cfg
end

function main(args=ARGS)
    cfg = parse_args(args)

    T_grid = collect(cfg.T_min:cfg.T_step:cfg.T_max)
    rho_grid = collect(cfg.rho_min:cfg.rho_step:cfg.rho_max)

    println("="^60)
    println("Phase pipeline CLI")
    println("="^60)
    println("time=$(now()) model_kind=$(cfg.model_kind) profile=$(cfg.profile)")
    println("T-grid: $(first(T_grid)) -> $(last(T_grid)) (n=$(length(T_grid)))")
    println("rho-grid: $(first(rho_grid)) -> $(last(rho_grid)) (n=$(length(rho_grid)))")

    result = Models.run_phase_pipeline(
        cfg.model_kind;
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

    println("\n完成:")
    println("  run_id = $(result.run_id)")
    println("  CEP found = $(result.cep.found)")
    println("  boundary_count = $(length(result.first_order_boundary))")
    println("  artifacts = $(result.artifact_paths)")
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
