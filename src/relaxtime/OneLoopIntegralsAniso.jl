"""
单圈积分在各向异性一阶修正下的修正项B0_correction和A_correction的计算函数
"""
module OneLoopIntegralsCorrection

export B0_correction, A_correction, A_aniso,
    IntegrationDiagnostics

include("../QuarkDistribution_Aniso.jl")
include("OneLoopIntegrals.jl")
include("../integration/GaussLegendre.jl")
include("../integration/IntervalQuadratureStrategies.jl")
using .GaussLegendre: transform_standard16, transform_standard32, gauleg, gausslegendre
using .PNJLQuarkDistributions_Aniso: correction_cos_theta_coefficient, distribution_aniso
using .OneLoopIntegrals: internal_momentum, EPS_K,
    energy_cutoff, singularity_k_positive, const_integral_term_A

""" k>0时
对x积分所得结果的Sokhotski–Plemelj formula实部部分(柯西主值积分)
被积函数展开到一阶,与x相关的项只多了x^2
对x的积分的被积函数可以抽象为x^2/(Ax+B)
积分变量x=cosθ从-1积到1
A对应的形参名:coeff_x
B对应的形参名:denominator_const
"""
function real_integral_tool(coeff_x::Float64, denominator_const::Float64)
    term1 = -2*denominator_const/coeff_x^2
    term2 = denominator_const^2/coeff_x^3
    log_arg = (coeff_x + denominator_const) / (coeff_x - denominator_const)
    term3 = log(abs(log_arg))
    result = term1 + term2 * term3
    return result
end

"""对x积分的虚部部分"""
function imag_integral_tool(coeff_x::Float64, denominator_const::Float64)
    return π *denominator_const ^ 2 / coeff_x^3
end

"""阶跃函数,确保奇点出现在积分范围内时才贡献虚部"""
@inline function heaviside_step(coeff_x::Float64, denominator_const::Float64, x_min::Float64, x_max::Float64)
    # 计算奇点位置
    x = -denominator_const / coeff_x 
    return (x >= x_min && x <= x_max) ? 1.0 : 0.0
end

"""coeff_x和denominator_const的计算函数"""
@inline function compute_coefficients(λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64)
    if k>EPS_K # k>0
        denominator_const = λ^2+2*λ*E + m^2 - m_prime^2-k^2
        p = internal_momentum(E,m)
        coeff_x = 2*p*k
        return coeff_x, denominator_const
    else # k=0
        denominator_const = λ^2+m^2-m_prime^2
        coeff_E = 2*λ
        return coeff_E, denominator_const
    end
end

"""k=0 时被积函数的分子（不含分母）

返回值 f(E) = (4/3) * p(E) * correction_cos_theta_coefficient(...)
此函数仅计算分子部分，便于在外部进行奇点减法或在不同分母下复用。
"""
function real_integrand_k_zero_numer(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    p = internal_momentum(E, m)
    return (4.0/3.0) * p * correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
end

"""k=0 时被积函数的分母

返回值 denom(E) = coeff_E * E + denominator_const
此处复用 `compute_coefficients` 的 k=0 分支以确保与其它代码一致。
"""
@inline function real_integrand_k_zero_denom(λ::Float64, m::Float64, m_prime::Float64, E::Float64)
    coeff_E, denominator_const = compute_coefficients(λ, 0.0, m, m_prime, E)
    return coeff_E * E + denominator_const
end

"""k=0时的积分实部被积函数（保留原名，用分子/分母组合）"""
function real_integrand_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    num = real_integrand_k_zero_numer(sign_, λ, m, m_prime, E, ξ, T, μ, Φ, Φbar)
    den = real_integrand_k_zero_denom(λ, m, m_prime, E)
    return isfinite(den) && isfinite(num) ? num / den : 0.0
end

"""k>0时的积分实部被积函数"""
function real_integrand_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    coeff_x, denominator_const = compute_coefficients(λ, k, m, m_prime, E)
    p = internal_momentum(E,m)
    # 含各向异性修正项
    real_ = 2.0*p*real_integral_tool(coeff_x, denominator_const)*
        correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
    return real_
end

"""k=0时的积分虚部函数"""
function imag_integrand_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    coeff_E, denominator_const = compute_coefficients(λ, 0.0, m, m_prime, 0.0)
    E_pole = -denominator_const / coeff_E
    Θ = heaviside_step(coeff_E, denominator_const, m, energy_cutoff(m))
    if Θ == 0.0
        return 0.0
    else
        p_pole = internal_momentum(E_pole,m)
        # 含各向异性修正项
        imag_ = π*(4.0/3.0)*p_pole/coeff_E*
            correction_cos_theta_coefficient(sign_, p_pole, m, μ, T, Φ, Φbar, ξ)
        return imag_
    end
end

"""k>0时的积分虚部被积函数"""
function imag_integrand_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)

    coeff_x, denominator_const = compute_coefficients(λ, k, m, m_prime, E)
    
    p = internal_momentum(E,m)
    # 含各向异性修正项
    imag_ = 2.0*p*imag_integral_tool(coeff_x, denominator_const)*
            correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
    return imag_
end

"""k=0时的 B0分量 含各向异性修正项的积分计算"""
function tilde_B0_correction_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64)
    Emin = m
    Emax = energy_cutoff(m)
    # 采用奇点减法：f(E) = (4/3)*p(E)*correction(...) ，integrand = f(E) / (coeff_E*E + denominator_const)
    coeff_E = 2.0 * λ
    denominator_const = λ^2 + m^2 - m_prime^2

    # 寻找奇点 E0，当 coeff_E != 0 时
    E0 = nothing
    if coeff_E != 0.0
        E0_val = -denominator_const / coeff_E
        if isfinite(E0_val) && (E0_val >= Emin) && (E0_val <= Emax)
            E0 = E0_val
        end
    end

    # 使用已拆分的分子函数 `real_integrand_k_zero_numer` 来获取 f(E)
    # 注意：不要在此重复实现分子，直接调用模块级 helper。 

    # 使用模块内预计算的标准 32 节点并仿射映射到 [Emin, Emax]
    # 提升默认节点数到 32 以提高在奇点附近的近似精度
    nodes, weights = transform_standard32(Emin, Emax)

    # 如果存在奇点，仍做奇点减法（计算 A），但不在节点处特殊替代极限值。
    # 使用 GL 节点处直接评估剩余函数 (f(E)-A)/den，节点恰好命中奇点的概率极小。
    if E0 !== nothing
        # 直接调用分子函数计算 A = f(E0)
        A = real_integrand_k_zero_numer(sign_, λ, m, m_prime, E0, ξ, T, μ, Φ, Φbar)
        vals = similar(nodes)
        @inbounds for i in eachindex(nodes)
            E = nodes[i]
            den = coeff_E * E + denominator_const
            v = (real_integrand_k_zero_numer(sign_, λ, m, m_prime, E, ξ, T, μ, Φ, Φbar) - A)
            vals[i] = isfinite(den) && isfinite(v) ? v / den : 0.0
        end
        real_part = sum(weights .* vals)
        # 将常数项 A 的主值积分解析部分加回：
        # PV ∫_Emin^Emax A / (coeff_E * E + denominator_const) dE = A/coeff_E * ln|coeff_E*E + denominator_const| |_Emin^Emax
        if coeff_E != 0.0
            real_part += A/coeff_E * (log(abs(coeff_E * Emax + denominator_const)) - log(abs(coeff_E * Emin + denominator_const)))
        end
        imag_part = imag_integrand_k_zero(sign_, λ, m, m_prime,
            ξ, T, μ, Φ, Φbar)
        return real_part, imag_part
    else
        # 无奇点或 coeff_E == 0：直接用 GL 对原始实部被积函数求积
        vals = similar(nodes)
        @inbounds for i in eachindex(nodes)
            E = nodes[i]
            vals[i] = real_integrand_k_zero(sign_, λ, m, m_prime, E,
                ξ, T, μ, Φ, Φbar)
        end
        real_part = sum(weights .* vals)
        imag_part = imag_integrand_k_zero(sign_, λ, m, m_prime,
            ξ, T, μ, Φ, Φbar)
        return real_part , imag_part
    end
end

# ============================================================================
# 辅助函数：根查找与区间构建
# ============================================================================
"""解析求解 A±B=0 的根

方程推导：
  A = 2*k*sqrt(E² - m²)
  B = λ² + 2λE + m² - m'² - k²
  
  A ± B = 0 两边平方后化简为二次方程：
  (k² - λ²)E² - λCE - (k²m² + C²/4) = 0
  其中 C = λ² + m² - m'² - k²

返回在 [Emin, Emax] 范围内的有效根。
"""
function find_roots_AB(λ::Float64, k::Float64, m::Float64, m_prime::Float64, Emin::Float64, Emax::Float64)
    C = λ^2 + m^2 - m_prime^2 - k^2
    
    a = k^2 - λ^2
    b = -λ * C
    c = -(k^2 * m^2 + C^2 / 4)
    
    # 处理 a ≈ 0 的情况 (k ≈ |λ|)
    if abs(a) < 1e-14
        if abs(b) < 1e-14
            return Float64[]  # 无解或恒等式
        end
        E = -c / b
        # 检查是否在范围内且满足原方程
        if E >= Emin && E <= Emax && E >= m
            # 验证原方程
            p = sqrt(max(0.0, E^2 - m^2))
            A = 2*k*p
            B = λ^2 + 2*λ*E + m^2 - m_prime^2 - k^2
            if abs(A + B) < 1e-10 || abs(A - B) < 1e-10
                return [E]
            end
        end
        return Float64[]
    end
    
    Δ = b^2 - 4*a*c
    
    if Δ < 0
        return Float64[]  # 无实根
    end
    
    sqrt_Δ = sqrt(Δ)
    E1 = (-b - sqrt_Δ) / (2*a)
    E2 = (-b + sqrt_Δ) / (2*a)
    
    # 筛选有效根：
    # 1. E >= m (物理约束，确保 p = sqrt(E²-m²) 为实数)
    # 2. E 在 [Emin, Emax] 范围内
    # 3. 验证原方程（平方可能引入伪根）
    valid_roots = Float64[]
    
    for E in [E1, E2]
        if E >= m && E >= Emin && E <= Emax
            # 验证原方程 A ± B = 0
            p = sqrt(E^2 - m^2)
            A = 2*k*p
            B = λ^2 + 2*λ*E + m^2 - m_prime^2 - k^2
            
            # A + B = 0 或 A - B = 0
            if abs(A + B) < 1e-10 || abs(A - B) < 1e-10
                push!(valid_roots, E)
            end
        end
    end
    
    return sort(unique(valid_roots))
end

"""根据根位置构建分割区间"""
function build_intervals_from_roots(roots::Vector{Float64}, Emin::Float64, Emax::Float64)
    δ = max(1e-10, 1e-8*(Emax-Emin))
    endpoints = [Emin]
    for r in roots
        push!(endpoints, r-δ)
        push!(endpoints, r+δ)
    end
    push!(endpoints, Emax)
    
    pairs = Tuple{Float64,Float64}[]
    for i in 1:2:length(endpoints)-1
        a, b = endpoints[i], endpoints[i+1]
        if b > a + 1e-16
            push!(pairs, (a, b))
        end
    end
    return pairs
end

# ============================================================================
# k>0 时的 B0 分量积分计算（增强版）
# ============================================================================
"""k>0时的 B0分量 含各向异性修正项的积分计算

参数:
- sign_: 符号类型 (:quark 或 :antiquark)
- λ, k, m, m_prime: 物理参数
- μ, T, Φ, Φbar, ξ: 热力学参数
- diagnostics: 是否返回诊断信息 (默认 false)

返回:
- (real_part, imag_part) 或 (real_part, imag_part, diagnostics)
"""
function tilde_B0_correction_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64;
    diagnostics::Bool=false)
    
    t_start = time()
    Emin = m
    Emax = energy_cutoff(m)

    # 被积函数闭包
    integrand_real(E) = real_integrand_k_positive(sign_, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar)
    integrand_imag(E) = imag_integrand_k_positive(sign_, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar)

    # 获取虚部积分区间
    intervals_imag, _ = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
    
    # 固定 hybrid：根据 A±B=0 的根分段，并按奇点位置选择映射
    real_part = 0.0
    roots = Float64[]
    intervals_real = Tuple{Float64,Float64}[]

    # 查找 A±B=0 的根并构造分段
    roots = find_roots_AB(λ, k, m, m_prime, Emin, Emax)
    intervals_real = build_intervals_from_roots(roots, Emin, Emax)

    # 混合策略：根据奇点位置选择映射与节点数
    n_intervals = length(intervals_real)
    for (idx, (a, b)) in enumerate(intervals_real)
        sing_pos = if n_intervals == 1
            isempty(roots) ? SING_NONE : SING_BOTH
        elseif idx == 1
            SING_RIGHT
        elseif idx == n_intervals
            SING_LEFT
        else
            SING_BOTH
        end

        n_nodes = if isempty(roots)
            16
        elseif sing_pos == SING_LEFT
            16
        else
            32
        end

        real_part += integrate_hybrid_interval(integrand_real, a, b, sing_pos; n=n_nodes)
    end
    
    # 计算虚部（使用标准 16 节点 GL，无分配）
    imag_part = 0.0
    if !isempty(intervals_imag)
        for (E1, E2) in intervals_imag
            half = (E2 - E1) / 2
            center = (E1 + E2) / 2
            @inbounds @simd for i in eachindex(_STD_16_NODES)
                x = center + half * _STD_16_NODES[i]
                w = _STD_16_WEIGHTS[i] * half
                v = integrand_imag(x)
                if isfinite(v)
                    imag_part += w * v
                end
            end
        end
        imag_part *= sign(λ)
    end
    
    t_elapsed = (time() - t_start) * 1000  # ms
    
    if diagnostics
        diag = IntegrationDiagnostics(
            STRATEGY_HYBRID, length(roots), roots,
            length(intervals_real), intervals_real,
            real_part, imag_part, t_elapsed
        )
        return real_part, imag_part, diag
    else
        return real_part, imag_part
    end
end

"""含各向异性修正项的 B0分量 积分计算"""
function tilde_B0_correction(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64)
    if k > EPS_K
        return tilde_B0_correction_k_positive(sign_, λ, k, m, m_prime, μ, T,
            Φ, Φbar, ξ)
    else
        return tilde_B0_correction_k_zero(sign_, λ, m, m_prime, μ, T,
            Φ, Φbar, ξ)
    end
end

"""
    B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
动量各向异性下B0的一阶修正
"""
function B0_correction(λ::Float64, k::Float64, m1::Float64, m2::Float64, μ1::Float64, μ2::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64)

    term1 = tilde_B0_correction(:quark, -λ, k, m1, m2, μ1, T, Φ, Φbar, ξ)
    term2 = tilde_B0_correction(:antiquark, λ, k, m1, m2, μ1, T, Φ, Φbar, ξ)
    term3 = tilde_B0_correction(:quark, λ, k, m2, m1, μ2, T, Φ, Φbar, ξ)
    term4 = tilde_B0_correction(:antiquark, -λ, k, m2, m1, μ2, T, Φ, Φbar, ξ)

    # 注意：B0 的 pm=-1(antiquark) 分支语义是 NJL 中的f(-E-μ)=1-f(E+μ)
    # 这里f(E+μ)是NJL中反粒子的分布函数，为避免 PNJL 分布在负能量下的数值溢出，
    # 实现上采将严格恒等式变换后的NJL形式转换到PNJL的分布函数f^-(E,μ)上计算：
    # 一阶各向异性修正只取 δf（对 ξ 的线性项），常数 1 无修正：
    #   δ f(−E-μ) = δ(1 − f^−(E, μ)) = − δ f^−(E, μ)
    # 因此在最终组合 B0 = B0^+(...) − B0^−(...) + ... 中，对应 term2/term4 出现“负负得正”。
    real_part = term1[1] + term2[1] + term3[1] + term4[1]
    imag_part = term1[2] + term2[2] + term3[2] + term4[2]
    return real_part, imag_part
end

# -------------------------------------
"""
    A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)
计算单线积分函数A对ξ的一阶修正项(需要加上各向同性项才能进行后续计算，即零阶项),需要传入预生成的动量的积分节点与权重
可以使用 GaussLegendre 模块中的 DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS (积分上限 10.0 fm⁻¹)
"""
function A_correction(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    ξ::Float64, nodes_p::Vector{Float64}, weights_p::Vector{Float64})
    integral = 0.0
    @inbounds for i in eachindex(nodes_p)
        node_p = nodes_p[i]
        weight_p = weights_p[i]
        E = sqrt(node_p^2 + m^2)
        quark_correction = correction_cos_theta_coefficient(:quark, node_p, m, μ, T, Φ, Φbar, ξ)
        antiquark_correction = correction_cos_theta_coefficient(:antiquark, node_p, m, μ, T, Φ, Φbar, ξ)
        integral += weight_p * node_p^2 / E * (quark_correction + antiquark_correction)
    end
    return 4.0 * integral / 3.0
end

"""
    A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
计算单线积分函数A在动量各向异性下的完整形式,需要传入预生成的动量和角度的积分节点与权重
可以使用 GaussLegendre 模块中的常量：
- DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS (p ∈ [0, 10] fm⁻¹)
- DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS (cosθ ∈ [-1, 1])
"""
function A_aniso(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    ξ::Float64, nodes_p::Vector{Float64}, weights_p::Vector{Float64},
    nodes_cosθ::Vector{Float64}, weights_cosθ::Vector{Float64})
    # 只计算分布函数的积分部分，常数项单独处理
    integral = 0.0
    @inbounds @simd for i in eachindex(nodes_p)
        node_p = nodes_p[i]
        weight_p = weights_p[i]
        # 对角度进行积分
        angle_integral_quark = 0.0
        angle_integral_antiquark = 0.0
        E = sqrt(node_p^2 + m^2)  # 使用各向同性能量作为分母
        @inbounds @simd for j in eachindex(nodes_cosθ)
            cosθ = nodes_cosθ[j]
            weight_cosθ = weights_cosθ[j]
            dist_quark = distribution_aniso(:quark, node_p, m, μ, T, Φ, Φbar, ξ, cosθ)
            dist_antiquark = distribution_aniso(:antiquark, node_p, m, μ, T, Φ, Φbar, ξ, cosθ)
            angle_integral_quark += weight_cosθ * node_p^2 / E * dist_quark
            angle_integral_antiquark += weight_cosθ * node_p^2 / E * dist_antiquark
        end
        integral += weight_p * (angle_integral_quark + angle_integral_antiquark)
    end
    # 分布函数项：angle integration给出因子2（cosθ权重和），φ积分贡献2π→归一化为2，总共4
    # 常数项：无角度依赖，需要完整4π因子
    return 2.0 * integral - 4.0 * const_integral_term_A(m)
end
end # module OneLoopIntegralsCorrection