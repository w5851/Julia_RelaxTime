#!/usr/bin/env julia
#=
PNJL 热力学导数计算脚本

计算质量导数和体粘滞系数相关量：
- ∂M/∂T, ∂M/∂μ（夸克有效质量导数）
- v_n² (等熵声速平方)
- ∂μ_B/∂T|_σ (等比熵线斜率)
- ζ/s 相关系数

用法：
    julia scripts/pnjl/calculate_derivatives.jl [options]

选项：
    --xi=0.0          各向异性参数
    --T_min=100       最低温度 (MeV)
    --T_max=200       最高温度 (MeV)
    --T_step=10       温度步长 (MeV)
    --mu_min=0        最低化学势 (MeV)
    --mu_max=300      最高化学势 (MeV)
    --mu_step=50      化学势步长 (MeV)
    --output_dir=...  输出目录
    --verbose         详细输出
    --help            显示帮助

输出文件：
    - derivatives.csv: 质量导数和热力学导数
    - bulk_viscosity.csv: 体粘滞系数相关量
=#

const HELP_TEXT = """
PNJL 热力学导数计算脚本

用法：julia scripts/pnjl/calculate_derivatives.jl [options]

选项：
    --xi=0.0          各向异性参数
    --T_min=100       最低温度 (MeV)
    --T_max=200       最高温度 (MeV)
    --T_step=10       温度步长 (MeV)
    --mu_min=0        最低化学势 (MeV)
    --mu_max=300      最高化学势 (MeV)
    --mu_step=50      化学势步长 (MeV)
    --output_dir=...  输出目录
    --p_num=24        动量积分节点数
    --t_num=8         角度积分节点数
    --verbose         详细输出
    --help            显示帮助

输出文件：
    - derivatives.csv: 质量导数和热力学导数
    - bulk_viscosity.csv: 体粘滞系数相关量
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Dates

# 加载模块
include(joinpath(@__DIR__, "..", "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pnjl", "PNJL.jl"))

using .PNJL
using .PNJL.ThermoDerivatives

# 单位转换
const hbarc = 197.327  # MeV·fm

# ============================================================================
# 配置
# ============================================================================

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "data", "outputs", "results", "pnjl")

struct DerivativesConfig
    xi::Float64
    T_min::Float64
    T_max::Float64
    T_step::Float64
    mu_min::Float64
    mu_max::Float64
    mu_step::Float64
    output_dir::String
    p_num::Int
    t_num::Int
    verbose::Bool
end

function parse_args(args)
    xi = 0.0
    T_min = 100.0
    T_max = 200.0
    T_step = 10.0
    mu_min = 0.0
    mu_max = 300.0
    mu_step = 50.0
    output_dir = DEFAULT_OUTPUT_DIR
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
        elseif startswith(arg, "--output_dir=")
            output_dir = arg[14:end]
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
            @warn "未知参数: $(arg)"
        end
    end
    
    return DerivativesConfig(
        xi, T_min, T_max, T_step, mu_min, mu_max, mu_step,
        output_dir, p_num, t_num, verbose
    )
end

# ============================================================================
# 主函数
# ============================================================================

function main(args=ARGS)
    config = parse_args(args)
    
    println("=" ^ 60)
    println("PNJL 热力学导数计算")
    println("=" ^ 60)
    println("时间: $(now())")
    println()
    println("参数配置:")
    println("  ξ = $(config.xi)")
    println("  T 范围: $(config.T_min) - $(config.T_max) MeV (步长 $(config.T_step))")
    println("  μ 范围: $(config.mu_min) - $(config.mu_max) MeV (步长 $(config.mu_step))")
    println("  积分节点: p_num=$(config.p_num), t_num=$(config.t_num)")
    println("  输出目录: $(config.output_dir)")
    println()
    
    mkpath(config.output_dir)
    
    # 构建参数网格
    T_values = collect(config.T_min:config.T_step:config.T_max)
    mu_values = collect(config.mu_min:config.mu_step:config.mu_max)
    
    total_points = length(T_values) * length(mu_values)
    println("扫描网格:")
    println("  温度点数: $(length(T_values))")
    println("  化学势点数: $(length(mu_values))")
    println("  总点数: $(total_points)")
    println()
    
    # Step 1: 计算质量导数
    derivatives_path = joinpath(config.output_dir, "derivatives_xi$(config.xi).csv")
    step1_mass_derivatives(config, T_values, mu_values, derivatives_path)
    
    # Step 2: 计算体粘滞系数
    bulk_path = joinpath(config.output_dir, "bulk_viscosity_xi$(config.xi).csv")
    step2_bulk_viscosity(config, T_values, mu_values, bulk_path)
    
    println("\n" * "=" ^ 60)
    println("计算完成!")
    println("=" ^ 60)
end


# ============================================================================
# Step 1: 质量导数计算
# ============================================================================

const DERIVATIVES_HEADER = join([
    "T_MeV", "mu_MeV", "xi",
    "M_u_MeV", "M_d_MeV", "M_s_MeV",
    "dM_u_dT", "dM_d_dT", "dM_s_dT",
    "dM_u_dmu", "dM_d_dmu", "dM_s_dmu",
    "dP_dT", "dP_dmu",
    "dEpsilon_dT", "dEpsilon_dmu",
    "dn_dT", "dn_dmu",
    "pressure_fm4", "energy_fm4", "entropy_fm3", "rho_norm",
    "converged"
], ",")

function step1_mass_derivatives(config::DerivativesConfig, T_values, mu_values, output_path)
    println("[Step 1] 质量导数计算")
    println("-" ^ 40)
    
    t_start = time()
    success_count = 0
    failure_count = 0
    
    open(output_path, "w") do io
        println(io, DERIVATIVES_HEADER)
        
        for T_MeV in T_values
            for μ_MeV in mu_values
                # 转换单位
                T_fm = T_MeV / hbarc
                μ_fm = μ_MeV / hbarc
                
                try
                    # 计算热力学导数
                    result = thermo_derivatives(T_fm, μ_fm;
                        xi=config.xi,
                        p_num=config.p_num,
                        t_num=config.t_num
                    )
                    
                    # 转换质量单位
                    M_MeV = result.masses .* hbarc
                    
                    # 写入结果
                    values = (
                        @sprintf("%.2f", T_MeV),
                        @sprintf("%.2f", μ_MeV),
                        @sprintf("%.2f", config.xi),
                        @sprintf("%.6f", M_MeV[1]),
                        @sprintf("%.6f", M_MeV[2]),
                        @sprintf("%.6f", M_MeV[3]),
                        @sprintf("%.6e", result.dM_dT[1]),
                        @sprintf("%.6e", result.dM_dT[2]),
                        @sprintf("%.6e", result.dM_dT[3]),
                        @sprintf("%.6e", result.dM_dmu[1]),
                        @sprintf("%.6e", result.dM_dmu[2]),
                        @sprintf("%.6e", result.dM_dmu[3]),
                        @sprintf("%.6e", result.dP_dT),
                        @sprintf("%.6e", result.dP_dmu),
                        @sprintf("%.6e", result.dEpsilon_dT),
                        @sprintf("%.6e", result.dEpsilon_dmu),
                        @sprintf("%.6e", result.dn_dT),
                        @sprintf("%.6e", result.dn_dmu),
                        @sprintf("%.6e", result.pressure),
                        @sprintf("%.6e", result.energy),
                        @sprintf("%.6e", result.entropy),
                        @sprintf("%.6f", result.rho_norm),
                        string(result.converged)
                    )
                    println(io, join(values, ","))
                    flush(io)
                    
                    success_count += 1
                    if config.verbose
                        print(".")
                    end
                    
                catch e
                    failure_count += 1
                    if config.verbose
                        print("x")
                        @warn "T=$(T_MeV), μ=$(μ_MeV) 失败: $(e)"
                    end
                end
            end
            
            if config.verbose
                println(" T=$(T_MeV) MeV")
            end
        end
    end
    
    t_elapsed = time() - t_start
    
    println("完成:")
    println("  成功: $(success_count) / $(success_count + failure_count)")
    println("  失败: $(failure_count)")
    println("  耗时: $(round(t_elapsed, digits=1)) 秒")
    println("  输出: $(output_path)")
end

# ============================================================================
# Step 2: 体粘滞系数计算
# ============================================================================

const BULK_HEADER = join([
    "T_MeV", "mu_MeV", "xi",
    "v_n_sq", "dmuB_dT_sigma",
    "M_u_MeV", "M_d_MeV", "M_s_MeV",
    "dM_u_dT", "dM_d_dT", "dM_s_dT",
    "dM_u_dmuB", "dM_d_dmuB", "dM_s_dmuB",
    "entropy_fm3", "n_B_fm3",
    "success"
], ",")

function step2_bulk_viscosity(config::DerivativesConfig, T_values, mu_values, output_path)
    println("\n[Step 2] 体粘滞系数计算")
    println("-" ^ 40)
    
    t_start = time()
    success_count = 0
    failure_count = 0
    
    open(output_path, "w") do io
        println(io, BULK_HEADER)
        
        for T_MeV in T_values
            for μ_MeV in mu_values
                # 转换单位
                T_fm = T_MeV / hbarc
                μ_fm = μ_MeV / hbarc
                
                try
                    # 计算体粘滞系数
                    result = bulk_viscosity_coefficients(T_fm, μ_fm;
                        xi=config.xi,
                        p_num=config.p_num,
                        t_num=config.t_num
                    )
                    
                    # 转换质量单位
                    M_MeV = result.masses .* hbarc
                    
                    # 检查结果有效性
                    is_valid = isfinite(result.v_n_sq) && isfinite(result.dμB_dT_sigma)
                    
                    # 写入结果
                    values = (
                        @sprintf("%.2f", T_MeV),
                        @sprintf("%.2f", μ_MeV),
                        @sprintf("%.2f", config.xi),
                        @sprintf("%.6e", result.v_n_sq),
                        @sprintf("%.6e", result.dμB_dT_sigma),
                        @sprintf("%.6f", M_MeV[1]),
                        @sprintf("%.6f", M_MeV[2]),
                        @sprintf("%.6f", M_MeV[3]),
                        @sprintf("%.6e", result.dM_dT[1]),
                        @sprintf("%.6e", result.dM_dT[2]),
                        @sprintf("%.6e", result.dM_dT[3]),
                        @sprintf("%.6e", result.dM_dμB[1]),
                        @sprintf("%.6e", result.dM_dμB[2]),
                        @sprintf("%.6e", result.dM_dμB[3]),
                        @sprintf("%.6e", result.s),
                        @sprintf("%.6e", result.n_B),
                        string(is_valid)
                    )
                    println(io, join(values, ","))
                    flush(io)
                    
                    success_count += 1
                    if config.verbose
                        print(".")
                    end
                    
                catch e
                    failure_count += 1
                    if config.verbose
                        print("x")
                        @warn "T=$(T_MeV), μ=$(μ_MeV) 失败: $(e)"
                    end
                end
            end
            
            if config.verbose
                println(" T=$(T_MeV) MeV")
            end
        end
    end
    
    t_elapsed = time() - t_start
    
    println("完成:")
    println("  成功: $(success_count) / $(success_count + failure_count)")
    println("  失败: $(failure_count)")
    println("  耗时: $(round(t_elapsed, digits=1)) 秒")
    println("  输出: $(output_path)")
end

# ============================================================================
# 入口
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
