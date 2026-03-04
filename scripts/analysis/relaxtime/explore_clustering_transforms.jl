# 探索不同的聚簇变换函数
# 目标：找到更适合单侧奇点的变换，进一步提升精度
using CSV, DataFrames, Printf
using SpecialFunctions: erf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
include(joinpath(@__DIR__, "../../../src/integration/GaussLegendre.jl"))
using .GaussLegendre: gauleg
using QuadGK: quadgk

# 测试参数
λ = -1.0; k = 0.01; m = 0.3; m_prime = 0.3; ξ = -0.2
T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)

integrand_full(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar)
val_ref, _ = quadgk(integrand_full, Emin, Emax; rtol=1e-12)
println(@sprintf("Reference: %.14g", val_ref))

# 找根
function find_roots()
    f_minus(E) = begin
        coeff_x, denom = OneLoopIntegralsCorrection.compute_coefficients(λ, k, m, m_prime, E)
        return denom - coeff_x
    end
    f_plus(E) = begin
        coeff_x, denom = OneLoopIntegralsCorrection.compute_coefficients(λ, k, m, m_prime, E)
        return denom + coeff_x
    end
    
    function bracket(f, a0, b0)
        xs = range(a0, stop=b0, length=2000)
        vals = [f(x) for x in xs]
        roots = Float64[]
        for i in 1:length(xs)-1
            if isfinite(vals[i]) && isfinite(vals[i+1]) && vals[i]*vals[i+1] < 0
                lo, hi = xs[i], xs[i+1]
                for _ in 1:60
                    mid = (lo+hi)/2
                    vm = f(mid)
                    if !isfinite(vm); lo = mid; continue; end
                    if sign(vm) == sign(vals[i]); lo = mid
                    else; hi = mid; end
                    if abs(hi-lo) < 1e-14 break end
                end
                push!(roots, (lo+hi)/2)
            end
        end
        return roots
    end
    return sort(unique(vcat(bracket(f_minus, Emin, Emax), bracket(f_plus, Emin, Emax))))
end

roots = find_roots()
println("Roots: ", roots)

# 构建区间（奇点在区间端点）
δ = 1e-10
intervals = [
    (Emin, roots[1] - δ),           # 奇点在右端
    (roots[1] + δ, roots[2] - δ),   # 奇点在两端
    (roots[2] + δ, Emax)            # 奇点在左端
]
println("Intervals: ", intervals)

# ============================================================================
# 变换函数定义
# ============================================================================

"""1. 标准 tanh 映射（对称聚簇）"""
function transform_tanh_symmetric(a, b, n; beta=8.0)
    us, ws = gauleg(-1.0, 1.0, n)
    tanh_beta = tanh(beta)
    phi(u) = tanh(beta*u) / tanh_beta
    phi_prime(u) = beta * (1 - tanh(beta*u)^2) / tanh_beta
    half = (b - a) / 2
    center = (a + b) / 2
    xs = center .+ half .* phi.(us)
    wx = ws .* half .* phi_prime.(us)
    return xs, wx
end

"""2. 单侧 tanh 映射（聚簇于右端 b）"""
function transform_tanh_right(a, b, n; beta=4.0)
    us, ws = gauleg(-1.0, 1.0, n)
    # 映射 u ∈ [-1,1] → t ∈ [0,1]，然后用 tanh 聚簇于 t=1
    # t = (1 + tanh(β*(u+1)/2 - β/2)) / (1 + tanh(β/2))
    # 简化：使用 sinh 变换
    tanh_beta = tanh(beta)
    # 将 u ∈ [-1,1] 映射到 [0,1]，聚簇于 1
    phi(u) = (1 + tanh(beta * (u + 1) / 2)) / (1 + tanh_beta) - 1 + 1
    # 更简单的方法：幂次变换
    # t = ((u+1)/2)^α, α > 1 聚簇于 0，α < 1 聚簇于 1
    alpha = 0.3  # 聚簇于右端
    t(u) = 1 - (1 - (u+1)/2)^(1/alpha)
    t_prime(u) = (1/alpha) * (1 - (u+1)/2)^(1/alpha - 1) / 2
    xs = a .+ (b - a) .* t.(us)
    wx = ws .* (b - a) .* t_prime.(us)
    return xs, wx
end

"""3. 单侧 tanh 映射（聚簇于左端 a）"""
function transform_tanh_left(a, b, n; beta=4.0)
    us, ws = gauleg(-1.0, 1.0, n)
    alpha = 0.3  # 聚簇于左端
    t(u) = ((u+1)/2)^(1/alpha)
    t_prime(u) = (1/alpha) * ((u+1)/2)^(1/alpha - 1) / 2
    xs = a .+ (b - a) .* t.(us)
    wx = ws .* (b - a) .* t_prime.(us)
    return xs, wx
end

"""4. 双指数 (DE/tanh-sinh) 变换"""
function transform_de(a, b, n; h=0.1)
    # 使用等距点而非 GL 点
    ks = -n÷2:n÷2
    ts = ks .* h
    # DE 变换: x = (a+b)/2 + (b-a)/2 * tanh(π/2 * sinh(t))
    half = (b - a) / 2
    center = (a + b) / 2
    xs = center .+ half .* tanh.(π/2 .* sinh.(ts))
    # 权重: w = h * (b-a)/2 * π/2 * cosh(t) * sech²(π/2 * sinh(t))
    wx = h .* half .* (π/2) .* cosh.(ts) .* (sech.(π/2 .* sinh.(ts))).^2
    return xs, wx
end

"""5. Gauss-Jacobi 风格变换（端点加权）"""
function transform_jacobi_style(a, b, n; alpha_left=0.0, alpha_right=0.5)
    # 使用幂次变换模拟 Jacobi 权重
    us, ws = gauleg(-1.0, 1.0, n)
    # 变换: t = (1+u)/2, 然后 x = a + (b-a) * t^(1+α_left) * (1-t)^(-α_right) 归一化
    # 简化版本：使用 Beta 分布 CDF 风格
    # 这里用简单的双幂次变换
    function transform(u)
        t = (u + 1) / 2  # t ∈ [0, 1]
        # 使用 sinh 变换聚簇两端
        s = sinh(3.0 * (t - 0.5)) / sinh(1.5) / 2 + 0.5
        return a + (b - a) * s
    end
    function transform_deriv(u)
        t = (u + 1) / 2
        return (b - a) * 3.0 * cosh(3.0 * (t - 0.5)) / sinh(1.5) / 2 / 2
    end
    xs = transform.(us)
    wx = ws .* transform_deriv.(us)
    return xs, wx
end

"""6. 指数聚簇（更激进的端点聚簇）"""
function transform_exp_cluster(a, b, n; gamma=3.0)
    us, ws = gauleg(-1.0, 1.0, n)
    # 使用 erf 变换
    # t = (1 + erf(γ*u)) / 2
    t(u) = (1 + erf(gamma * u)) / 2
    t_prime(u) = gamma * 2 / sqrt(π) * exp(-gamma^2 * u^2) / 2
    xs = a .+ (b - a) .* t.(us)
    wx = ws .* (b - a) .* t_prime.(us)
    return xs, wx
end

"""7. 自适应聚簇（根据奇点位置调整）"""
function transform_adaptive(a, b, n, singularity_pos::Symbol; beta=6.0)
    us, ws = gauleg(-1.0, 1.0, n)
    
    if singularity_pos == :left
        # 聚簇于左端 a - 使用幂次变换
        alpha = 0.3
        t = ((us .+ 1) ./ 2).^(1/alpha)
        t_prime = (1/alpha) .* ((us .+ 1) ./ 2).^(1/alpha - 1) ./ 2
        xs = a .+ (b - a) .* t
        wx = ws .* (b - a) .* t_prime
    elseif singularity_pos == :right
        # 聚簇于右端 b
        alpha = 0.3
        t = 1 .- (1 .- (us .+ 1) ./ 2).^(1/alpha)
        t_prime = (1/alpha) .* (1 .- (us .+ 1) ./ 2).^(1/alpha - 1) ./ 2
        xs = a .+ (b - a) .* t
        wx = ws .* (b - a) .* t_prime
    else  # :both
        # 对称聚簇（原始 tanh）
        return transform_tanh_symmetric(a, b, n; beta=beta)
    end
    return xs, wx
end

# ============================================================================
# 测试各种变换
# ============================================================================

println("\n" * "=" ^ 70)
println("Testing different clustering transforms")
println("=" ^ 70)

rows = DataFrame(transform=String[], n=Int[], param=Float64[], value=Float64[], relerr=Float64[])

# 对每个区间测试
for (idx, (a, b)) in enumerate(intervals)
    # 确定奇点位置
    sing_pos = if idx == 1
        :right  # 第一个区间，奇点在右端
    elseif idx == 3
        :left   # 最后一个区间，奇点在左端
    else
        :both   # 中间区间，两端都有奇点
    end
    
    # 区间参考值
    sub_ref, _ = quadgk(integrand_full, a, b; rtol=1e-12)
    
    println(@sprintf("\nInterval %d: [%.6f, %.6f] (singularity: %s)", idx, a, b, sing_pos))
    println(@sprintf("  Reference: %.12g", sub_ref))
    
    for n in [32, 64]
        # 1. 标准 tanh
        for beta in [4.0, 8.0]
            xs, wx = transform_tanh_symmetric(a, b, n; beta=beta)
            val = sum(wx .* integrand_full.(xs))
            err = abs((val - sub_ref) / sub_ref)
            push!(rows, ("tanh_symmetric", n, beta, val, err))
            if n == 64
                println(@sprintf("  tanh_sym(β=%.0f, n=%d): relerr=%.2e", beta, n, err))
            end
        end
        
        # 2. 自适应聚簇
        for beta in [4.0, 6.0, 8.0]
            xs, wx = transform_adaptive(a, b, n, sing_pos; beta=beta)
            val = sum(wx .* integrand_full.(xs))
            err = abs((val - sub_ref) / sub_ref)
            push!(rows, ("adaptive_$(sing_pos)", n, beta, val, err))
            if n == 64
                println(@sprintf("  adaptive(β=%.0f, n=%d): relerr=%.2e", beta, n, err))
            end
        end
        
        # 3. DE 变换
        for h in [0.08, 0.1, 0.12]
            xs, wx = transform_de(a, b, n; h=h)
            val = sum(wx .* integrand_full.(xs))
            err = abs((val - sub_ref) / sub_ref)
            push!(rows, ("DE", n, h, val, err))
            if n == 64 && h == 0.1
                println(@sprintf("  DE(h=%.2f, n=%d): relerr=%.2e", h, n, err))
            end
        end
        
        # 4. erf 聚簇
        for gamma in [2.0, 3.0, 4.0]
            xs, wx = transform_exp_cluster(a, b, n; gamma=gamma)
            val = sum(wx .* integrand_full.(xs))
            err = abs((val - sub_ref) / sub_ref)
            push!(rows, ("erf_cluster", n, gamma, val, err))
            if n == 64 && gamma == 3.0
                println(@sprintf("  erf(γ=%.0f, n=%d): relerr=%.2e", gamma, n, err))
            end
        end
    end
end

# 全局积分测试
println("\n" * "=" ^ 70)
println("Full integral comparison (all intervals combined)")
println("=" ^ 70)

function integrate_all(transform_func, n, param; adaptive_positions=[:right, :both, :left])
    total = 0.0
    for (idx, (a, b)) in enumerate(intervals)
        if transform_func == transform_adaptive
            xs, wx = transform_func(a, b, n, adaptive_positions[idx]; beta=param)
        elseif transform_func == transform_de
            xs, wx = transform_func(a, b, n; h=param)
        elseif transform_func == transform_exp_cluster
            xs, wx = transform_func(a, b, n; gamma=param)
        else
            xs, wx = transform_func(a, b, n; beta=param)
        end
        total += sum(wx .* integrand_full.(xs))
    end
    return total
end

println(@sprintf("\nReference: %.14g\n", val_ref))

for n in [32, 64, 128]
    println(@sprintf("n = %d:", n))
    
    # tanh symmetric
    for beta in [4.0, 8.0, 12.0]
        val = integrate_all(transform_tanh_symmetric, n, beta)
        err = abs((val - val_ref) / val_ref)
        println(@sprintf("  tanh_symmetric(β=%.0f): relerr=%.2e", beta, err))
    end
    
    # adaptive
    for beta in [4.0, 6.0, 8.0]
        val = integrate_all(transform_adaptive, n, beta)
        err = abs((val - val_ref) / val_ref)
        println(@sprintf("  adaptive(β=%.0f):       relerr=%.2e", beta, err))
    end
    
    # DE
    for h in [0.08, 0.1]
        val = integrate_all(transform_de, n, h)
        err = abs((val - val_ref) / val_ref)
        println(@sprintf("  DE(h=%.2f):            relerr=%.2e", h, err))
    end
    
    # erf
    for gamma in [2.5, 3.0]
        val = integrate_all(transform_exp_cluster, n, gamma)
        err = abs((val - val_ref) / val_ref)
        println(@sprintf("  erf(γ=%.1f):           relerr=%.2e", gamma, err))
    end
    
    println()
end

# 保存结果
outcsv = joinpath(@__DIR__, "../../data/processed/results/relaxtime/clustering_transforms_comparison.csv")
mkpath(dirname(outcsv))
CSV.write(outcsv, rows)
println("Wrote results to ", outcsv)
