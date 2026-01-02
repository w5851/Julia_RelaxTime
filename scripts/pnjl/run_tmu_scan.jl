#!/usr/bin/env julia
#=
PNJL T-μ 参数空间扫描脚本

使用 PhaseAwareContinuitySeed 进行 T-μ 参数空间扫描，
自动加载相变线数据实现相变感知的初值选择。

输出热力学量：P, s, ε, ρ, M_u, M_d, M_s 等

用法：
    julia scripts/pnjl/run_tmu_scan.jl [options]

选项：
    --xi=0.0          各向异性参数
    --T_min=50        最低温度 (MeV)
    --T_max=200       最高温度 (MeV)
    --T_step=10       温度步长 (MeV)
    --mu_min=0        最低化学势 (MeV)
    --mu_max=400      最高化学势 (MeV)
    --mu_step=10      化学势步长 (MeV)
    --output=...      输出文件路径（默认 data/outputs/results/pnjl/tmu_scan.csv）
    --resume          断点续扫（默认启用）
    --overwrite       覆盖已有文件
    --no_phase_aware  禁用相变感知策略
    --p_num=24        动量积分节点数
    --t_num=8         角度积分节点数
    --verbose         详细输出
    --help            显示帮助

示例：
    # 基本扫描
    julia scripts/pnjl/run_tmu_scan.jl

    # 自定义参数范围
    julia scripts/pnjl/run_tmu_scan.jl --xi=0.2 --T_min=50 --T_max=150 --mu_max=350

    # 高分辨率扫描
    julia scripts/pnjl/run_tmu_scan.jl --T_step=5 --mu_step=5 --verbose
=#

const HELP_TEXT = """
PNJL T-μ 参数空间扫描脚本

用法：julia scripts/pnjl/run_tmu_scan.jl [options]

选项：
    --xi=0.0          各向异性参数
    --T_min=50        最低温度 (MeV)
    --T_max=200       最高温度 (MeV)
    --T_step=10       温度步长 (MeV)
    --mu_min=0        最低化学势 (MeV)
    --mu_max=400      最高化学势 (MeV)
    --mu_step=10      化学势步长 (MeV)
    --output=...      输出文件路径
    --resume          断点续扫（默认启用）
    --overwrite       覆盖已有文件
    --no_phase_aware  禁用相变感知策略
    --p_num=24        动量积分节点数
    --t_num=8         角度积分节点数
    --verbose         详细输出
    --help            显示帮助
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Dates

# 加载模块
include(joinpath(@__DIR__, "..", "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pnjl", "PNJL.jl"))

using .PNJL
using .PNJL.TmuScan

# ============================================================================
# 配置
# ============================================================================

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "data", "outputs", "results", "pnjl")

struct TmuScanConfig
    xi::Float64
    T_min::Float64
    T_max::Float64
    T_step::Float64
    mu_min::Float64
    mu_max::Float64
    mu_step::Float64
    output_path::String
    resume::Bool
    overwrite::Bool
    use_phase_aware::Bool
    p_num::Int
    t_num::Int
    verbose::Bool
end

function parse_args(args)
    xi = 0.0
    T_min = 50.0
    T_max = 200.0
    T_step = 10.0
    mu_min = 0.0
    mu_max = 400.0
    mu_step = 10.0
    output_path = ""
    resume = true
    overwrite = false
    use_phase_aware = true
    p_num = 24
    t_num = 8
    verbose = false
    
    for arg in args
        if startswith(arg, "--xi=")
            xi = parse(Float64, arg[6:end])
        elseif startswith(arg, "--T_min=")
            T_min = parse(Float64, arg[9:end])
        elseif startswith(arg, "--T_max=")
            T_max = parse(Float64, arg[9:end])
        elseif startswith(arg, "--T_step=")
            T_step = parse(Float64, arg[10:end])
        elseif startswith(arg, "--mu_min=")
            mu_min = parse(Float64, arg[10:end])
        elseif startswith(arg, "--mu_max=")
            mu_max = parse(Float64, arg[10:end])
        elseif startswith(arg, "--mu_step=")
            mu_step = parse(Float64, arg[11:end])
        elseif startswith(arg, "--output=")
            output_path = arg[10:end]
        elseif arg == "--resume"
            resume = true
        elseif arg == "--overwrite"
            overwrite = true
            resume = false
        elseif arg == "--no_phase_aware"
            use_phase_aware = false
        elseif startswith(arg, "--p_num=")
            p_num = parse(Int, arg[9:end])
        elseif startswith(arg, "--t_num=")
            t_num = parse(Int, arg[9:end])
        elseif arg == "--verbose"
            verbose = true
        elseif arg == "--help" || arg == "-h"
            println(HELP_TEXT)
            exit(0)
        else
            @warn "未知参数: $arg"
        end
    end
    
    # 默认输出路径
    if isempty(output_path)
        output_path = joinpath(DEFAULT_OUTPUT_DIR, "tmu_scan_xi$(xi).csv")
    end
    
    return TmuScanConfig(
        xi, T_min, T_max, T_step, mu_min, mu_max, mu_step,
        output_path, resume, overwrite, use_phase_aware, p_num, t_num, verbose
    )
end

# ============================================================================
# 主函数
# ============================================================================

function main(args=ARGS)
    config = parse_args(args)
    
    println("=" ^ 60)
    println("PNJL T-μ 参数空间扫描")
    println("=" ^ 60)
    println("时间: $(now())")
    println()
    println("参数配置:")
    println("  ξ = $(config.xi)")
    println("  T 范围: $(config.T_min) - $(config.T_max) MeV (步长 $(config.T_step))")
    println("  μ 范围: $(config.mu_min) - $(config.mu_max) MeV (步长 $(config.mu_step))")
    println("  积分节点: p_num=$(config.p_num), t_num=$(config.t_num)")
    println("  相变感知: $(config.use_phase_aware ? "启用" : "禁用")")
    println("  断点续扫: $(config.resume ? "启用" : "禁用")")
    println("  输出文件: $(config.output_path)")
    println()
    
    # 构建参数网格
    T_values = collect(config.T_min:config.T_step:config.T_max)
    mu_values = collect(config.mu_min:config.mu_step:config.mu_max)
    
    total_points = length(T_values) * length(mu_values)
    println("扫描网格:")
    println("  温度点数: $(length(T_values))")
    println("  化学势点数: $(length(mu_values))")
    println("  总点数: $(total_points)")
    println()
    
    # 检查相变线数据
    if config.use_phase_aware
        boundary_path = joinpath(@__DIR__, "..", "..", "data", "reference", "pnjl", "boundary.csv")
        if isfile(boundary_path)
            println("相变线数据: $(boundary_path) ✓")
        else
            @warn "相变线数据不存在: $(boundary_path)，将使用普通连续性跟踪"
        end
    end
    println()
    
    # 进度回调
    progress_cb = nothing
    completed = Ref(0)
    failed = Ref(0)
    
    if config.verbose
        progress_cb = (point, result) -> begin
            completed[] += 1
            if result !== nothing && result.converged
                print(".")
            else
                failed[] += 1
                print("x")
            end
            if completed[] % 50 == 0
                println(" [$(completed[])/$(total_points)]")
            end
        end
    end
    
    # 执行扫描
    println("开始扫描...")
    println("-" ^ 40)
    t_start = time()
    
    result = run_tmu_scan(
        T_values = T_values,
        mu_values = mu_values,
        xi_values = [config.xi],
        output_path = config.output_path,
        overwrite = config.overwrite,
        resume = config.resume,
        use_phase_aware = config.use_phase_aware,
        p_num = config.p_num,
        t_num = config.t_num,
        progress_cb = progress_cb
    )
    
    t_elapsed = time() - t_start
    
    if config.verbose
        println()  # 换行
    end
    
    # 输出统计
    println("-" ^ 40)
    println("扫描完成!")
    println()
    println("统计:")
    println("  总点数: $(result.total)")
    println("  成功: $(result.success) ($(round(100*result.success/max(1,result.total-result.skipped), digits=1))%)")
    println("  失败: $(result.failure)")
    println("  跳过: $(result.skipped)")
    println("  耗时: $(round(t_elapsed, digits=1)) 秒")
    if result.total - result.skipped > 0
        avg_time = t_elapsed / (result.total - result.skipped)
        println("  平均: $(round(avg_time * 1000, digits=1)) ms/点")
    end
    println()
    println("输出文件: $(result.output)")
    println("=" ^ 60)
    
    return result
end

# ============================================================================
# 入口
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
