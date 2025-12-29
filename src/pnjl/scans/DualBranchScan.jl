"""
    DualBranchScan

双分支扫描模块，用于处理一阶相变区域。

## 核心思想
在一阶相变区域，能隙方程存在多值解（强子相和夸克相）。
通过双向扫描可以：
1. 分别追踪两个分支的解
2. 比较 Ω 选择物理稳定解
3. 定位相变点（Ω 交叉点）

## 扫描策略
- 强子分支：从低 μ 向高 μ 扫描，使用强子相初值
- 夸克分支：从高 μ 向低 μ 扫描，使用夸克相初值
- 每个 μ 点比较两分支的 Ω，选择较小者

## 使用示例
```julia
using PNJL.DualBranchScan

# 双分支扫描
result = run_dual_branch_scan(T_mev=100.0, mu_range=0.0:10.0:400.0)

# 获取相变点
transition = find_phase_transition(result)
println("相变点: μ_c = \$(transition.mu_c) MeV")
```
"""
module DualBranchScan

using Printf
using StaticArrays

# 导入新架构模块
using ..Constants_PNJL: ħc_MeV_fm
using ..ConstraintModes: FixedMu, ConstraintMode
using ..SeedStrategies: SeedStrategy, DefaultSeed, ContinuitySeed
using ..SeedStrategies: get_seed, update!, reset!
using ..SeedStrategies: HADRON_SEED_5, QUARK_SEED_5, HIGH_TEMP_SEED_5
using ..ImplicitSolver: solve, SolverResult

export run_dual_branch_scan, find_phase_transition, merge_branches
export DualBranchResult, BranchPoint, PhaseTransitionInfo

# ============================================================================
# 数据结构
# ============================================================================

"""单个分支点的结果"""
struct BranchPoint
    mu_mev::Float64
    converged::Bool
    omega::Float64
    pressure::Float64
    rho_norm::Float64
    entropy::Float64
    energy::Float64
    x_state::SVector{5, Float64}
    masses::SVector{3, Float64}
    iterations::Int
    residual_norm::Float64
end

"""双分支扫描结果"""
struct DualBranchResult
    T_mev::Float64
    xi::Float64
    mu_values::Vector{Float64}
    hadron_branch::Vector{Union{Nothing, BranchPoint}}
    quark_branch::Vector{Union{Nothing, BranchPoint}}
    physical_branch::Vector{Union{Nothing, BranchPoint}}  # Ω 最小的分支
end

"""相变信息"""
struct PhaseTransitionInfo
    found::Bool
    mu_c::Float64                 # 相变点化学势 (MeV) - Ω 交叉点
    omega_at_transition::Float64
    mu_hadron_spinodal::Float64   # 强子分支亚稳态边界 (MeV) - 强子相消失点
    mu_quark_spinodal::Float64    # 夸克分支亚稳态边界 (MeV) - 夸克相出现点
    coexistence_region::Tuple{Float64, Float64}  # 共存区 (μ_min, μ_max)
end

# ============================================================================
# 主扫描函数
# ============================================================================

"""
    run_dual_branch_scan(; T_mev, mu_range, xi=0.0, kwargs...) -> DualBranchResult

执行双分支扫描。

# 参数
- `T_mev`: 温度 (MeV)
- `mu_range`: 化学势范围 (MeV)，如 0.0:10.0:400.0
- `xi`: 各向异性参数
- `p_num`, `t_num`: 积分节点数
- `verbose`: 是否输出进度信息

# 返回
DualBranchResult 包含两个分支的完整结果和物理分支选择。
"""
function run_dual_branch_scan(;
    T_mev::Real,
    mu_range,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    verbose::Bool=false,
    nlsolve_kwargs...
)
    mu_values = collect(Float64, mu_range)
    n_points = length(mu_values)
    
    T_fm = T_mev / ħc_MeV_fm
    
    # 初始化结果数组
    hadron_branch = Vector{Union{Nothing, BranchPoint}}(nothing, n_points)
    quark_branch = Vector{Union{Nothing, BranchPoint}}(nothing, n_points)
    
    # ========== 强子分支：从低 μ 向高 μ ==========
    verbose && println("扫描强子分支 (μ: $(mu_values[1]) → $(mu_values[end]) MeV)...")
    
    hadron_tracker = ContinuitySeed(fallback=DefaultSeed(phase_hint=:hadron))
    prev_hadron_result = nothing
    
    for (i, mu_mev) in enumerate(mu_values)
        μ_fm = mu_mev / ħc_MeV_fm
        result = _solve_point(T_fm, μ_fm, xi, hadron_tracker; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
        
        if result !== nothing && result.converged
            # 检查是否发生了解跳跃（序参量突变）
            if prev_hadron_result !== nothing && _is_solution_jump(prev_hadron_result, result)
                verbose && println("  强子分支在 μ = $mu_mev MeV 发生跳跃，停止追踪")
                break
            end
            
            hadron_branch[i] = _to_branch_point(mu_mev, result)
            update!(hadron_tracker, result.solution)
            prev_hadron_result = result
        else
            # 分支断裂，停止追踪
            verbose && println("  强子分支在 μ = $mu_mev MeV 断裂")
            break
        end
    end
    
    # ========== 夸克分支：从高 μ 向低 μ ==========
    verbose && println("扫描夸克分支 (μ: $(mu_values[end]) → $(mu_values[1]) MeV)...")
    
    quark_tracker = ContinuitySeed(fallback=DefaultSeed(phase_hint=:quark))
    prev_quark_result = nothing
    
    for i in n_points:-1:1
        mu_mev = mu_values[i]
        μ_fm = mu_mev / ħc_MeV_fm
        result = _solve_point(T_fm, μ_fm, xi, quark_tracker; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
        
        if result !== nothing && result.converged
            # 检查是否发生了解跳跃
            if prev_quark_result !== nothing && _is_solution_jump(prev_quark_result, result)
                verbose && println("  夸克分支在 μ = $mu_mev MeV 发生跳跃，停止追踪")
                break
            end
            
            quark_branch[i] = _to_branch_point(mu_mev, result)
            update!(quark_tracker, result.solution)
            prev_quark_result = result
        else
            # 分支断裂，停止追踪
            verbose && println("  夸克分支在 μ = $mu_mev MeV 断裂")
            break
        end
    end
    
    # ========== 选择物理分支（Ω 最小） ==========
    physical_branch = _select_physical_branch(hadron_branch, quark_branch)
    
    return DualBranchResult(
        Float64(T_mev),
        Float64(xi),
        mu_values,
        hadron_branch,
        quark_branch,
        physical_branch
    )
end

# ============================================================================
# 相变点查找
# ============================================================================

"""
    find_phase_transition(result::DualBranchResult) -> PhaseTransitionInfo

从双分支扫描结果中查找相变点和亚稳态边界。

返回三个关键化学势：
- `mu_c`: 一阶相变点（Ω 交叉点）
- `mu_hadron_spinodal`: 强子分支亚稳态边界（强子相消失点，spinodal line）
- `mu_quark_spinodal`: 夸克分支亚稳态边界（夸克相出现点，spinodal line）

共存区定义为 [mu_quark_spinodal, mu_hadron_spinodal]。
"""
function find_phase_transition(result::DualBranchResult)
    mu_values = result.mu_values
    hadron = result.hadron_branch
    quark = result.quark_branch
    
    # 找强子分支的亚稳态边界（最后一个有效点）
    mu_hadron_spinodal = _find_branch_endpoint(hadron, mu_values, :forward)
    
    # 找夸克分支的亚稳态边界（最后一个有效点，从高μ向低μ扫描）
    mu_quark_spinodal = _find_branch_endpoint(quark, mu_values, :backward)
    
    # 找到两分支都有解的区域（共存区）
    coexist_indices = Int[]
    for i in eachindex(mu_values)
        if hadron[i] !== nothing && quark[i] !== nothing
            # 额外检查：确保两分支不是同一个解
            if !_is_same_solution(hadron[i], quark[i])
                push!(coexist_indices, i)
            end
        end
    end
    
    if isempty(coexist_indices)
        return PhaseTransitionInfo(
            false, NaN, NaN,
            mu_hadron_spinodal,
            mu_quark_spinodal,
            (NaN, NaN)
        )
    end
    
    # 在共存区寻找 Ω 交叉点
    mu_c, omega_c = _find_omega_crossing(
        mu_values[coexist_indices],
        [hadron[i].omega for i in coexist_indices],
        [quark[i].omega for i in coexist_indices]
    )
    
    coexist_region = (mu_values[coexist_indices[1]], mu_values[coexist_indices[end]])
    
    return PhaseTransitionInfo(
        !isnan(mu_c),
        mu_c,
        omega_c,
        mu_hadron_spinodal,
        mu_quark_spinodal,
        coexist_region
    )
end

"""
    merge_branches(result::DualBranchResult; output_path=nothing) -> DataFrame-like

合并双分支结果，输出物理解。
"""
function merge_branches(result::DualBranchResult; output_path::Union{Nothing, String}=nothing)
    rows = NamedTuple[]
    
    for (i, mu_mev) in enumerate(result.mu_values)
        pt = result.physical_branch[i]
        h_pt = result.hadron_branch[i]
        q_pt = result.quark_branch[i]
        
        if pt === nothing
            continue
        end
        
        # 判断当前点属于哪个分支
        branch = :unknown
        if h_pt !== nothing && q_pt !== nothing
            # 检查两分支是否是同一个解
            if _is_same_solution(h_pt, q_pt)
                branch = :same  # 两分支收敛到同一解
            elseif abs(pt.omega - h_pt.omega) < 1e-10
                branch = :hadron
            else
                branch = :quark
            end
        elseif h_pt !== nothing
            branch = :hadron
        elseif q_pt !== nothing
            branch = :quark
        end
        
        # 计算 ΔΩ（仅当两分支不同时有意义）
        delta_omega = NaN
        if h_pt !== nothing && q_pt !== nothing && !_is_same_solution(h_pt, q_pt)
            delta_omega = h_pt.omega - q_pt.omega
        end
        
        row = (
            T_MeV = result.T_mev,
            mu_MeV = mu_mev,
            xi = result.xi,
            branch = branch,
            omega = pt.omega,
            pressure = pt.pressure,
            rho = pt.rho_norm,
            entropy = pt.entropy,
            energy = pt.energy,
            phi_u = pt.x_state[1],
            phi_d = pt.x_state[2],
            phi_s = pt.x_state[3],
            Phi1 = pt.x_state[4],
            Phi2 = pt.x_state[5],
            M_u_MeV = pt.masses[1] * ħc_MeV_fm,
            M_d_MeV = pt.masses[2] * ħc_MeV_fm,
            M_s_MeV = pt.masses[3] * ħc_MeV_fm,
            delta_omega = delta_omega,
        )
        push!(rows, row)
    end
    
    # 可选：写入文件
    if output_path !== nothing
        _write_merged_csv(output_path, rows)
    end
    
    return rows
end

# ============================================================================
# 内部辅助函数
# ============================================================================

"""固定种子策略"""
struct _FixedSeedStrategy <: SeedStrategy
    seed::Vector{Float64}
end

import ..SeedStrategies: get_seed
function get_seed(s::_FixedSeedStrategy, ::AbstractVector, ::ConstraintMode)
    return copy(s.seed)
end

"""单点求解"""
function _solve_point(T_fm, μ_fm, xi, tracker::ContinuitySeed; p_num, t_num, nlsolve_kwargs...)
    seed = get_seed(tracker, [T_fm, μ_fm], FixedMu())
    strategy = _FixedSeedStrategy(seed)
    
    try
        result = solve(FixedMu(), T_fm, μ_fm;
            xi=xi,
            seed_strategy=strategy,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...
        )
        return result
    catch
        return nothing
    end
end

"""转换为 BranchPoint"""
function _to_branch_point(mu_mev::Float64, result::SolverResult)
    return BranchPoint(
        mu_mev,
        result.converged,
        result.omega,
        result.pressure,
        result.rho_norm,
        result.entropy,
        result.energy,
        result.x_state,
        result.masses,
        result.iterations,
        result.residual_norm
    )
end

"""选择物理分支（Ω 最小）"""
function _select_physical_branch(hadron, quark)
    n = length(hadron)
    physical = Vector{Union{Nothing, BranchPoint}}(nothing, n)
    
    for i in eachindex(hadron)
        h = hadron[i]
        q = quark[i]
        
        if h === nothing && q === nothing
            physical[i] = nothing
        elseif h === nothing
            physical[i] = q
        elseif q === nothing
            physical[i] = h
        else
            # 两分支都有解
            # 检查是否实际上是同一个解（通过比较序参量）
            if _is_same_solution(h, q)
                # 同一个解，优先选强子分支（保持连续性）
                physical[i] = h
            else
                # 不同解，选 Ω 较小的
                physical[i] = h.omega <= q.omega ? h : q
            end
        end
    end
    
    return physical
end

"""判断两个解是否实际上是同一个解"""
function _is_same_solution(h::BranchPoint, q::BranchPoint; tol=0.01)
    # 比较序参量（φ_u 和 Φ）
    # 如果差异小于阈值，认为是同一个解
    delta_phi = abs(h.x_state[1] - q.x_state[1])
    delta_Phi = abs(h.x_state[4] - q.x_state[4])
    
    # 也比较有效质量
    delta_M = abs(h.masses[1] - q.masses[1])
    
    return delta_phi < tol && delta_Phi < tol && delta_M < tol * 197.327  # 质量用 MeV
end

"""检测解是否发生跳跃（用于连续性追踪）"""
function _is_solution_jump(prev::SolverResult, curr::SolverResult; 
                           phi_threshold=0.5, mass_threshold=50.0)
    # 比较 φ_u 的变化
    delta_phi = abs(prev.x_state[1] - curr.x_state[1])
    
    # 比较有效质量的变化 (MeV)
    delta_M = abs(prev.masses[1] - curr.masses[1]) * 197.327
    
    # 如果变化超过阈值，认为发生了跳跃
    return delta_phi > phi_threshold || delta_M > mass_threshold
end

"""找分支端点"""
function _find_branch_endpoint(branch, mu_values, direction::Symbol)
    if direction == :forward
        for i in reverse(eachindex(branch))
            if branch[i] !== nothing
                return mu_values[i]
            end
        end
    else
        for i in eachindex(branch)
            if branch[i] !== nothing
                return mu_values[i]
            end
        end
    end
    return NaN
end

"""找 Ω 交叉点（线性插值）"""
function _find_omega_crossing(mu_vals, omega_hadron, omega_quark)
    n = length(mu_vals)
    
    # 计算 Δω = ω_hadron - ω_quark
    delta = omega_hadron .- omega_quark
    
    # 过滤掉 Δω 太小的点（两分支实际上是同一个解）
    # 使用相对阈值：|Δω| > 1e-6 * |ω|
    significant_indices = Int[]
    for i in 1:n
        avg_omega = 0.5 * (abs(omega_hadron[i]) + abs(omega_quark[i]))
        threshold = max(1e-6 * avg_omega, 1e-10)
        if abs(delta[i]) > threshold
            push!(significant_indices, i)
        end
    end
    
    if length(significant_indices) < 2
        return NaN, NaN
    end
    
    # 在显著差异的点中寻找符号变化
    for j in 1:(length(significant_indices)-1)
        i1 = significant_indices[j]
        i2 = significant_indices[j+1]
        
        if delta[i1] * delta[i2] < 0
            # 线性插值找交叉点
            t = delta[i1] / (delta[i1] - delta[i2])
            mu_c = mu_vals[i1] + t * (mu_vals[i2] - mu_vals[i1])
            omega_c = omega_hadron[i1] + t * (omega_hadron[i2] - omega_hadron[i1])
            return mu_c, omega_c
        end
    end
    
    return NaN, NaN
end

"""写入合并后的 CSV"""
function _write_merged_csv(path, rows)
    mkpath(dirname(path))
    
    header = join([
        "T_MeV", "mu_MeV", "xi", "branch",
        "omega", "pressure", "rho", "entropy", "energy",
        "phi_u", "phi_d", "phi_s", "Phi1", "Phi2",
        "M_u_MeV", "M_d_MeV", "M_s_MeV", "delta_omega"
    ], ",")
    
    open(path, "w") do io
        println(io, header)
        for row in rows
            values = [
                @sprintf("%.6f", row.T_MeV),
                @sprintf("%.6f", row.mu_MeV),
                @sprintf("%.6f", row.xi),
                string(row.branch),
                @sprintf("%.10f", row.omega),
                @sprintf("%.10f", row.pressure),
                @sprintf("%.10f", row.rho),
                @sprintf("%.10f", row.entropy),
                @sprintf("%.10f", row.energy),
                @sprintf("%.10f", row.phi_u),
                @sprintf("%.10f", row.phi_d),
                @sprintf("%.10f", row.phi_s),
                @sprintf("%.10f", row.Phi1),
                @sprintf("%.10f", row.Phi2),
                @sprintf("%.6f", row.M_u_MeV),
                @sprintf("%.6f", row.M_d_MeV),
                @sprintf("%.6f", row.M_s_MeV),
                isnan(row.delta_omega) ? "NaN" : @sprintf("%.10f", row.delta_omega),
            ]
            println(io, join(values, ","))
        end
    end
end

# ============================================================================
# 便捷函数：多温度扫描
# ============================================================================

"""
    scan_phase_diagram(; T_range, mu_range, kwargs...) -> Vector{PhaseTransitionInfo}

扫描多个温度下的相变点，构建相变线。
"""
function scan_phase_diagram(;
    T_range,
    mu_range,
    xi::Real=0.0,
    output_dir::Union{Nothing, String}=nothing,
    verbose::Bool=true,
    kwargs...
)
    transitions = PhaseTransitionInfo[]
    
    for T_mev in T_range
        verbose && println("\n处理 T = $T_mev MeV...")
        
        result = run_dual_branch_scan(T_mev=T_mev, mu_range=mu_range, xi=xi, verbose=verbose; kwargs...)
        transition = find_phase_transition(result)
        push!(transitions, transition)
        
        if transition.found
            verbose && println("  相变点: μ_c = $(round(transition.mu_c, digits=2)) MeV")
        else
            verbose && println("  未找到相变点（可能是 crossover）")
        end
        
        # 可选：保存每个温度的详细结果
        if output_dir !== nothing
            output_path = joinpath(output_dir, "dual_branch_T$(Int(T_mev)).csv")
            merge_branches(result; output_path=output_path)
        end
    end
    
    return transitions
end

export scan_phase_diagram

end # module DualBranchScan
