#!/usr/bin/env julia
"""
PNJL 相结构计算脚本

计算并输出：
1. CEP (临界终点)
2. 一阶相变线 (T-ρ 扫描 + Maxwell 构造)
3. Spinodal 数据（亚稳态边界）
4. Crossover 线（手征和退禁闭）

输出文件保存到 data/reference/pnjl/

用法：
    julia scripts/pnjl/calculate_phase_structure.jl [options]

选项：
    --xi=0.0          各向异性参数
    --T_min=50        最低温度 (MeV)
    --T_max=200       最高温度 (MeV)
    --T_step=10       温度步长 (MeV)
    --rho_max=4.0     最大密度 (ρ/ρ₀)
    --rho_step=0.05   密度步长
    --output_dir=...  输出目录
    --skip_trho       跳过 T-ρ 扫描（使用已有数据）
    --skip_crossover  跳过 crossover 计算
    --verbose         详细输出
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Dates

# 加载模块
include(joinpath(@__DIR__, "..", "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pnjl", "PNJL.jl"))

using .PNJL
using .PNJL.TrhoScan
using .PNJL.PhaseTransition: SShapeResult, detect_s_shape, MaxwellResult, maxwell_construction
using .PNJL.PhaseTransition: detect_crossover, CrossoverResult, scan_crossover_line

# 单位转换
const hbarc = 197.327  # MeV·fm

# ============================================================================
# 配置
# ============================================================================

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "data", "reference", "pnjl")

struct PhaseStructureConfig
    xi::Float64
    T_min::Float64
    T_max::Float64
    T_step::Float64
    rho_max::Float64
    rho_step::Float64
    output_dir::String
    skip_trho::Bool
    skip_crossover::Bool
    verbose::Bool
end

function parse_args(args)
    xi = 0.0
    T_min = 50.0
    T_max = 200.0
    T_step = 10.0
    rho_max = 4.0
    rho_step = 0.05
    output_dir = DEFAULT_OUTPUT_DIR
    skip_trho = false
    skip_crossover = false
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
        elseif startswith(arg, "--rho_max=")
            rho_max = parse(Float64, arg[11:end])
        elseif startswith(arg, "--rho_step=")
            rho_step = parse(Float64, arg[12:end])
        elseif startswith(arg, "--output_dir=")
            output_dir = arg[14:end]
        elseif arg == "--skip_trho"
            skip_trho = true
        elseif arg == "--skip_crossover"
            skip_crossover = true
        elseif arg == "--verbose"
            verbose = true
        elseif arg == "--help" || arg == "-h"
            println("用法: julia calculate_phase_structure.jl [options]")
            println("选项:")
            println("  --xi=0.0          各向异性参数")
            println("  --T_min=50        最低温度 (MeV)")
            println("  --T_max=200       最高温度 (MeV)")
            println("  --T_step=10       温度步长 (MeV)")
            println("  --rho_max=4.0     最大密度 (ρ/ρ₀)")
            println("  --rho_step=0.05   密度步长")
            println("  --output_dir=...  输出目录")
            println("  --skip_trho       跳过 T-ρ 扫描")
            println("  --skip_crossover  跳过 crossover 计算")
            println("  --verbose         详细输出")
            exit(0)
        end
    end
    
    return PhaseStructureConfig(xi, T_min, T_max, T_step, rho_max, rho_step, 
                                 output_dir, skip_trho, skip_crossover, verbose)
end

# ============================================================================
# 主函数
# ============================================================================

function main(args=ARGS)
    config = parse_args(args)
    
    println("=" ^ 60)
    println("PNJL 相结构计算")
    println("=" ^ 60)
    println("时间: $(now())")
    println("参数:")
    println("  xi = $(config.xi)")
    println("  T 范围: $(config.T_min) - $(config.T_max) MeV (步长 $(config.T_step))")
    println("  ρ 范围: 0 - $(config.rho_max) ρ₀ (步长 $(config.rho_step))")
    println("  输出目录: $(config.output_dir)")
    println()
    
    mkpath(config.output_dir)
    
    # Step 1: T-ρ 扫描
    trho_path = joinpath(config.output_dir, "trho_scan_xi$(config.xi).csv")
    curves = step1_trho_scan(config, trho_path)
    
    # Step 2: CEP 搜索
    cep_result = step2_find_cep(config, curves)
    
    # Step 3: Maxwell 构造
    boundary_results = step3_maxwell_construction(config, curves)
    
    # Step 4: 保存结果（包括 spinodal）
    step4_save_results(config, cep_result, boundary_results, curves)
    
    # Step 5: Crossover 计算
    if !config.skip_crossover
        step5_crossover_calculation(config, cep_result)
    else
        println("\n[Step 5] 跳过 Crossover 计算")
    end
    
    println("\n" * "=" ^ 60)
    println("计算完成!")
    println("=" ^ 60)
end

# ============================================================================
# Step 1: T-ρ 扫描
# ============================================================================

function step1_trho_scan(config::PhaseStructureConfig, output_path::String)
    println("\n[Step 1] T-ρ 扫描")
    println("-" ^ 40)
    
    if config.skip_trho && isfile(output_path)
        println("跳过扫描，加载已有数据: $output_path")
        return load_curves_from_csv(output_path, config.xi)
    end
    
    T_values = collect(config.T_min:config.T_step:config.T_max)
    rho_values = collect(0.0:config.rho_step:config.rho_max)
    
    println("温度点数: $(length(T_values))")
    println("密度点数: $(length(rho_values))")
    println("总点数: $(length(T_values) * length(rho_values))")
    
    # 进度回调
    progress_cb = nothing
    if config.verbose
        progress_cb = (point, result) -> begin
            if result !== nothing && result.converged
                print(".")
            else
                print("x")
            end
        end
    end
    
    println("\n开始扫描...")
    t_start = time()
    
    result = run_trho_scan(
        T_values = T_values,
        rho_values = rho_values,
        xi_values = [config.xi],
        output_path = output_path,
        overwrite = true,
        reverse_rho = true,  # 反向扫描避免 ρ=0 奇异点问题
        progress_cb = progress_cb
    )
    
    t_elapsed = time() - t_start
    
    println("\n扫描完成:")
    println("  成功: $(result.success) / $(result.total)")
    println("  失败: $(result.failure)")
    println("  耗时: $(round(t_elapsed, digits=1)) 秒")
    
    # 加载曲线数据
    return load_curves_from_csv(output_path, config.xi)
end

function load_curves_from_csv(path::String, xi::Float64)
    curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    isfile(path) || return curves
    
    # 按温度分组
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    
    for line in eachline(path)
        startswith(line, "T_MeV") && continue
        isempty(strip(line)) && continue
        
        cols = split(line, ',')
        length(cols) < 21 && continue
        
        # TrhoScan 输出格式:
        # T_MeV(1), rho(2), xi(3), mu_u(4), mu_d(5), mu_s(6), mu_avg(7), ...
        # converged(21)
        T = tryparse(Float64, cols[1])
        rho = tryparse(Float64, cols[2])
        xi_val = tryparse(Float64, cols[3])
        mu = tryparse(Float64, cols[7])  # mu_avg_MeV
        converged = lowercase(strip(cols[21])) in ("true", "1")
        
        T === nothing && continue
        rho === nothing && continue
        mu === nothing && continue
        !converged && continue
        xi_val === nothing && continue
        abs(xi_val - xi) > 1e-6 && continue
        
        if !haskey(grouped, T)
            grouped[T] = Tuple{Float64, Float64}[]
        end
        push!(grouped[T], (mu, rho))
    end
    
    # 转换为 curves 格式
    for (T, samples) in grouped
        length(samples) < 3 && continue
        mu_vals = Float64[s[1] for s in samples]
        rho_vals = Float64[s[2] for s in samples]
        curves[T] = (mu_vals, rho_vals)
    end
    
    println("  加载了 $(length(curves)) 个温度点的曲线数据")
    
    return curves
end

# ============================================================================
# Step 2: CEP 搜索
# ============================================================================

"""
CEP 搜索：使用二分法细化 CEP 位置

算法：
1. 找到初始区间 [T_low, T_high]，其中 T_low 有 S 形，T_high 无 S 形
2. 二分法细化：计算 T_mid 的曲线，判断是否有 S 形
3. 重复直到区间宽度 < tol
"""
function step2_find_cep(config::PhaseStructureConfig, curves;
                        tol::Float64=0.01,  # 温度精度 (MeV)
                        max_bisect_iter::Int=20)
    println("\n[Step 2] CEP 搜索")
    println("-" ^ 40)
    
    if isempty(curves)
        println("警告: 无曲线数据，跳过 CEP 搜索")
        return (has_cep=false, T_cep=NaN, mu_cep=NaN)
    end
    
    println("可用温度点: $(length(curves))")
    
    # 按温度排序
    temperatures = sort(collect(keys(curves)))
    
    # 找到初始区间：最后一个有 S 形的温度和第一个没有 S 形的温度
    last_with_s = nothing
    first_without_s = nothing
    s_shape_cache = Dict{Float64, SShapeResult}()
    
    for T in temperatures
        mu_vals, rho_vals = curves[T]
        result = detect_s_shape(mu_vals, rho_vals)
        s_shape_cache[T] = result
        
        if result.has_s_shape
            last_with_s = (T, result)
        elseif last_with_s !== nothing && first_without_s === nothing
            first_without_s = T
        end
    end
    
    if last_with_s === nothing
        println("未找到 S 形曲线 (可能全为 crossover)")
        return (has_cep=false, T_cep=NaN, mu_cep=NaN)
    end
    
    T_low, res_low = last_with_s
    T_high = first_without_s !== nothing ? first_without_s : T_low + config.T_step
    
    println("初始区间: [$T_low, $T_high] MeV")
    
    # 二分法细化
    bisect_count = 0
    while (T_high - T_low) > tol && bisect_count < max_bisect_iter
        T_mid = (T_low + T_high) / 2
        
        # 检查是否已有该温度的曲线
        if haskey(curves, T_mid)
            mu_vals, rho_vals = curves[T_mid]
        else
            # 需要计算新的曲线（通过线性插值估计，或者实际计算）
            # 这里使用线性插值作为近似
            mu_vals, rho_vals = _interpolate_curve(curves, T_mid)
            if mu_vals === nothing
                println("  无法插值 T=$T_mid MeV 的曲线，停止二分")
                break
            end
            curves[T_mid] = (mu_vals, rho_vals)
        end
        
        result = detect_s_shape(mu_vals, rho_vals)
        s_shape_cache[T_mid] = result
        
        if result.has_s_shape
            T_low = T_mid
            res_low = result
            if config.verbose
                println("  T=$T_mid MeV: 有 S 形 → 更新下界")
            end
        else
            T_high = T_mid
            if config.verbose
                println("  T=$T_mid MeV: 无 S 形 → 更新上界")
            end
        end
        
        bisect_count += 1
    end
    
    # CEP 估计为区间中点
    T_cep = (T_low + T_high) / 2
    
    # μ_CEP 估计为 spinodal 中点
    mu_hadron = something(res_low.mu_spinodal_hadron, NaN)
    mu_quark = something(res_low.mu_spinodal_quark, NaN)
    mu_cep = (mu_hadron + mu_quark) / 2
    
    println("找到 CEP ($(bisect_count) 次二分):")
    println("  T_CEP ≈ $(round(T_cep, digits=2)) MeV (区间: $(round(T_low, digits=2)) - $(round(T_high, digits=2)))")
    println("  μ_CEP ≈ $(round(mu_cep, digits=2)) MeV")
    
    return (has_cep=true, T_cep=T_cep, mu_cep=mu_cep)
end

"""线性插值获取中间温度的曲线"""
function _interpolate_curve(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T_target::Float64)
    temps = sort(collect(keys(curves)))
    
    # 找到相邻的温度点
    T_below = nothing
    T_above = nothing
    for T in temps
        if T < T_target
            T_below = T
        elseif T > T_target && T_above === nothing
            T_above = T
        end
    end
    
    (T_below === nothing || T_above === nothing) && return nothing, nothing
    
    mu_below, rho_below = curves[T_below]
    mu_above, rho_above = curves[T_above]
    
    # 需要相同的 ρ 采样点
    length(rho_below) == length(rho_above) || return nothing, nothing
    
    # 线性插值 μ
    α = (T_target - T_below) / (T_above - T_below)
    mu_interp = mu_below .+ α .* (mu_above .- mu_below)
    
    return mu_interp, rho_below
end

# ============================================================================
# Step 3: Maxwell 构造
# ============================================================================

function step3_maxwell_construction(config::PhaseStructureConfig, curves)
    println("\n[Step 3] Maxwell 构造")
    println("-" ^ 40)
    
    if isempty(curves)
        println("警告: 无曲线数据，跳过 Maxwell 构造")
        return Dict{Float64, MaxwellResult}()
    end
    
    results = Dict{Float64, MaxwellResult}()
    
    for T in sort(collect(keys(curves)))
        mu_vals, rho_vals = curves[T]
        # 使用模块化的 maxwell_construction 算法
        result = maxwell_construction(mu_vals, rho_vals;
            min_samples=8,           # 降低最小样本要求
            detect_min_points=5,     # S 形检测最小点数
            detect_eps=1e-6,
            candidate_steps=64,
            max_iter=60,
            tol_area=1e-4
        )
        results[T] = result
    end
    
    success_count = count(r -> r.converged, values(results))
    println("Maxwell 构造结果: $(success_count) / $(length(results)) 成功")
    
    if config.verbose
        println("\n详细结果:")
        for T in sort(collect(keys(results)))
            r = results[T]
            if r.converged
                area_str = r.area_residual !== nothing ? @sprintf("%.2e", r.area_residual) : "N/A"
                println("  T=$(T) MeV: μ_c=$(round(r.mu_transition, digits=2)) MeV, " *
                       "ρ_hadron=$(round(r.rho_hadron, digits=3)), ρ_quark=$(round(r.rho_quark, digits=3)), " *
                       "残差=$(area_str), 迭代=$(r.iterations)")
            else
                reason = get(r.details, :reason, "unknown")
                println("  T=$(T) MeV: 失败 ($(reason))")
            end
        end
    end
    
    return results
end

# ============================================================================
# Step 4: 保存结果
# ============================================================================

function step4_save_results(config::PhaseStructureConfig, cep_result, boundary_results, curves)
    println("\n[Step 4] 保存结果")
    println("-" ^ 40)
    
    # 保存 CEP
    cep_path = joinpath(config.output_dir, "cep.csv")
    save_cep(cep_path, config.xi, cep_result)
    println("CEP 数据: $cep_path")
    
    # 保存相变线
    boundary_path = joinpath(config.output_dir, "boundary.csv")
    save_boundary(boundary_path, config.xi, boundary_results)
    println("相变线数据: $boundary_path")
    
    # 保存 spinodal（亚稳态边界）
    spinodal_path = joinpath(config.output_dir, "spinodals.csv")
    save_spinodals(spinodal_path, config.xi, curves)
    println("Spinodal 数据: $spinodal_path")
end

function save_cep(path::String, xi::Float64, result)
    # 读取现有数据（如果存在）
    existing = Dict{Float64, Tuple{Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 3 || continue
            xi_val = tryparse(Float64, cols[1])
            T_cep = tryparse(Float64, cols[2])
            mu_cep = tryparse(Float64, cols[3])
            xi_val === nothing && continue
            existing[xi_val] = (T_cep, mu_cep)
        end
    end
    
    # 更新当前 xi 的数据
    if result.has_cep
        existing[xi] = (result.T_cep, result.mu_cep)
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,T_CEP_MeV,mu_CEP_MeV")
        for xi_val in sort(collect(keys(existing)))
            T_cep, mu_cep = existing[xi_val]
            println(io, "$xi_val,$T_cep,$mu_cep")
        end
    end
end

function save_boundary(path::String, xi::Float64, results::Dict{Float64, MaxwellResult})
    # 读取现有数据（如果存在）
    existing = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 5 || continue
            xi_val = tryparse(Float64, cols[1])
            T = tryparse(Float64, cols[2])
            mu = tryparse(Float64, cols[3])
            rho_hadron = tryparse(Float64, cols[4])
            rho_quark = tryparse(Float64, cols[5])
            xi_val === nothing && continue
            existing[(xi_val, T)] = (mu, rho_hadron, rho_quark)
        end
    end
    
    # 更新当前 xi 的数据
    for (T, r) in results
        if r.converged && r.mu_transition !== nothing && r.rho_hadron !== nothing && r.rho_quark !== nothing
            existing[(xi, T)] = (r.mu_transition, r.rho_hadron, r.rho_quark)
        end
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark")
        for key in sort(collect(keys(existing)))
            xi_val, T = key
            mu, rho_hadron, rho_quark = existing[key]
            println(io, "$xi_val,$T,$mu,$rho_hadron,$rho_quark")
        end
    end
end

"""保存 spinodal（亚稳态边界）数据"""
function save_spinodals(path::String, xi::Float64, curves)
    # 读取现有数据（如果存在）
    existing = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 6 || continue
            xi_val = tryparse(Float64, cols[1])
            T = tryparse(Float64, cols[2])
            mu_low = tryparse(Float64, cols[3])
            mu_high = tryparse(Float64, cols[4])
            rho_low = tryparse(Float64, cols[5])
            rho_high = tryparse(Float64, cols[6])
            xi_val === nothing && continue
            existing[(xi_val, T)] = (mu_low, mu_high, rho_low, rho_high)
        end
    end
    
    # 从曲线数据提取 spinodal
    for (T, (mu_vals, rho_vals)) in curves
        result = detect_s_shape(mu_vals, rho_vals)
        if result.has_s_shape && 
           result.mu_spinodal_hadron !== nothing && 
           result.mu_spinodal_quark !== nothing &&
           result.rho_spinodal_hadron !== nothing &&
           result.rho_spinodal_quark !== nothing
            existing[(xi, T)] = (
                result.mu_spinodal_hadron,
                result.mu_spinodal_quark,
                result.rho_spinodal_hadron,
                result.rho_spinodal_quark
            )
        end
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark")
        for key in sort(collect(keys(existing)))
            xi_val, T = key
            mu_low, mu_high, rho_hadron, rho_quark = existing[key]
            println(io, "$xi_val,$T,$mu_low,$mu_high,$rho_hadron,$rho_quark")
        end
    end
end

# ============================================================================
# Step 5: Crossover 计算
# ============================================================================

"""
Crossover 计算：扫描 μ = 0 到 μ_CEP 的 crossover 线

计算手征 crossover (φ_u) 和退禁闭 crossover (Φ)
"""
function step5_crossover_calculation(config::PhaseStructureConfig, cep_result)
    println("\n[Step 5] Crossover 计算")
    println("-" ^ 40)
    
    # 确定 μ 扫描范围
    if cep_result.has_cep
        μ_max_MeV = cep_result.mu_cep
        println("μ 范围: 0 - $(round(μ_max_MeV, digits=1)) MeV (CEP)")
    else
        # 如果没有 CEP，使用默认范围
        μ_max_MeV = 300.0
        println("μ 范围: 0 - $(μ_max_MeV) MeV (默认)")
    end
    
    # 转换为 fm⁻¹
    μ_max_fm = μ_max_MeV / hbarc
    
    # 温度搜索范围 (fm⁻¹)
    T_min_fm = config.T_min / hbarc
    T_max_fm = config.T_max / hbarc
    
    # μ 采样点数
    n_mu = max(10, Int(ceil(μ_max_MeV / 20)))  # 约每 20 MeV 一个点
    
    println("温度搜索范围: $(config.T_min) - $(config.T_max) MeV")
    println("μ 采样点数: $n_mu")
    println()
    
    # 扫描手征 crossover
    println("计算手征 crossover (φ_u)...")
    t_start = time()
    chiral_results = scan_crossover_line(
        (0.0, μ_max_fm, n_mu), (T_min_fm, T_max_fm);
        method=:inflection, variable=:phi_u, xi=config.xi
    )
    t_chiral = time() - t_start
    chiral_success = count(r -> r.converged, chiral_results)
    println("  完成: $(chiral_success)/$n_mu 成功, 耗时 $(round(t_chiral, digits=1)) 秒")
    
    # 扫描退禁闭 crossover
    println("计算退禁闭 crossover (Φ)...")
    t_start = time()
    deconf_results = scan_crossover_line(
        (0.0, μ_max_fm, n_mu), (T_min_fm, T_max_fm);
        method=:inflection, variable=:Phi, xi=config.xi
    )
    t_deconf = time() - t_start
    deconf_success = count(r -> r.converged, deconf_results)
    println("  完成: $(deconf_success)/$n_mu 成功, 耗时 $(round(t_deconf, digits=1)) 秒")
    
    # 保存结果
    crossover_path = joinpath(config.output_dir, "crossover.csv")
    save_crossover(crossover_path, config.xi, chiral_results, deconf_results)
    println("\nCrossover 数据: $crossover_path")
end

"""保存 crossover 数据"""
function save_crossover(path::String, xi::Float64, 
                        chiral_results::Vector, deconf_results::Vector)
    # 读取现有数据（如果存在）
    existing = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 4 || continue
            xi_val = tryparse(Float64, cols[1])
            mu = tryparse(Float64, cols[2])
            T_chiral = tryparse(Float64, cols[3])
            T_deconf = tryparse(Float64, cols[4])
            xi_val === nothing && continue
            existing[(xi_val, mu)] = (something(T_chiral, NaN), something(T_deconf, NaN))
        end
    end
    
    # 更新当前 xi 的数据
    n = min(length(chiral_results), length(deconf_results))
    for i in 1:n
        μ_MeV = chiral_results[i].mu_fm * hbarc
        T_chiral_MeV = chiral_results[i].converged ? 
                       chiral_results[i].T_crossover_fm * hbarc : NaN
        T_deconf_MeV = deconf_results[i].converged ? 
                       deconf_results[i].T_crossover_fm * hbarc : NaN
        existing[(xi, μ_MeV)] = (T_chiral_MeV, T_deconf_MeV)
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,mu_MeV,T_crossover_chiral_MeV,T_crossover_deconf_MeV")
        for key in sort(collect(keys(existing)))
            xi_val, mu = key
            T_chiral, T_deconf = existing[key]
            T_chiral_str = isnan(T_chiral) ? "" : string(T_chiral)
            T_deconf_str = isnan(T_deconf) ? "" : string(T_deconf)
            println(io, "$xi_val,$mu,$T_chiral_str,$T_deconf_str")
        end
    end
end

# ============================================================================
# 入口
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
