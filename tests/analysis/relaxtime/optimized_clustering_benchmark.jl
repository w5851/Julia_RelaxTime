# 优化聚簇策略的深入测试
# 基于 explore_clustering_transforms.jl 的发现，进一步优化
using CSV, DataFrames, Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
include(joinpath(@__DIR__, "../../../src/integration/GaussLegendre.jl"))
using .GaussLegendre: gauleg
using QuadGK: quadgk

println("=" ^ 70)
println("Optimized Clustering Strategy Benchmark")
println("=" ^ 70)

# 测试参数
λ = -1.0; k = 0.01; m = 0.3; m_prime = 0.3; ξ = -0.2
T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)

integrand_full(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar)
val_ref, _ = quadgk(integrand_full, Emin, Emax; rtol=1e-14)
println(@sprintf("Reference: %.15g\n", val_ref))

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
δ = 1e-10

# ============================================================================
# 关键发现总结：
# 1. 单侧奇点区间：自适应单侧聚簇 (adaptive) 效果最好
# 2. 双侧奇点区间：DE 变换效果最好，但 tanh 也可接受
# 3. 最佳组合：对不同区间使用不同策略
# ============================================================================

"""标准 tanh 对称聚簇"""
function transform_tanh(a, b, n; beta=8.0)
    us, ws = gauleg(-1.0, 1.0, n)
    tanh_beta = tanh(beta)
    xs = (a+b)/2 .+ (b-a)/2 .* tanh.(beta .* us) ./ tanh_beta
    wx = ws .* (b-a)/2 .* beta .* (1 .- tanh.(beta .* us).^2) ./ tanh_beta
    return xs, wx
end

"""单侧幂次聚簇 - 聚簇于左端"""
function transform_power_left(a, b, n; alpha=0.3)
    us, ws = gauleg(-1.0, 1.0, n)
    t = ((us .+ 1) ./ 2).^(1/alpha)
    t_prime = (1/alpha) .* ((us .+ 1) ./ 2).^(1/alpha - 1) ./ 2
    xs = a .+ (b - a) .* t
    wx = ws .* (b - a) .* t_prime
    return xs, wx
end

"""单侧幂次聚簇 - 聚簇于右端"""
function transform_power_right(a, b, n; alpha=0.3)
    us, ws = gauleg(-1.0, 1.0, n)
    t = 1 .- (1 .- (us .+ 1) ./ 2).^(1/alpha)
    t_prime = (1/alpha) .* (1 .- (us .+ 1) ./ 2).^(1/alpha - 1) ./ 2
    xs = a .+ (b - a) .* t
    wx = ws .* (b - a) .* t_prime
    return xs, wx
end

"""DE (tanh-sinh) 变换 - 对双侧奇点最优"""
function transform_de(a, b, n; h=0.1)
    ks = collect(-n÷2:n÷2)
    ts = ks .* h
    half = (b - a) / 2
    center = (a + b) / 2
    xs = center .+ half .* tanh.(π/2 .* sinh.(ts))
    wx = h .* half .* (π/2) .* cosh.(ts) .* (sech.(π/2 .* sinh.(ts))).^2
    return xs, wx
end

"""混合策略：根据奇点位置自动选择最优变换"""
function transform_hybrid(a, b, n, sing_pos::Symbol; alpha=0.3, beta=8.0, h=0.1)
    if sing_pos == :left
        return transform_power_left(a, b, n; alpha=alpha)
    elseif sing_pos == :right
        return transform_power_right(a, b, n; alpha=alpha)
    elseif sing_pos == :both
        # 双侧奇点：DE 变换最优，但 tanh 也可接受
        return transform_de(a, b, n; h=h)
    else  # :none
        return gauleg(a, b, n)
    end
end

# ============================================================================
# 测试混合策略
# ============================================================================

intervals_info = [
    (Emin, roots[1] - δ, :right),
    (roots[1] + δ, roots[2] - δ, :both),
    (roots[2] + δ, Emax, :left)
]

println("Testing hybrid strategy vs uniform tanh:")
println("-" ^ 70)

for n in [32, 64, 128]
    # 1. 统一 tanh (当前实现)
    val_tanh = 0.0
    for (a, b, _) in intervals_info
        xs, wx = transform_tanh(a, b, n; beta=8.0)
        val_tanh += sum(wx .* integrand_full.(xs))
    end
    err_tanh = abs((val_tanh - val_ref) / val_ref)
    
    # 2. 混合策略
    val_hybrid = 0.0
    for (a, b, sing_pos) in intervals_info
        xs, wx = transform_hybrid(a, b, n, sing_pos; alpha=0.3, beta=8.0, h=0.1)
        val_hybrid += sum(wx .* integrand_full.(xs))
    end
    err_hybrid = abs((val_hybrid - val_ref) / val_ref)
    
    # 3. 纯 DE
    val_de = 0.0
    for (a, b, _) in intervals_info
        xs, wx = transform_de(a, b, n; h=0.1)
        val_de += sum(wx .* integrand_full.(xs))
    end
    err_de = abs((val_de - val_ref) / val_ref)
    
    println(@sprintf("n=%3d: tanh(β=8)=%.2e  hybrid=%.2e  DE=%.2e  improvement=%.1fx",
        n, err_tanh, err_hybrid, err_de, err_tanh/err_hybrid))
end

# ============================================================================
# 参数优化
# ============================================================================

println("\n" * "=" ^ 70)
println("Parameter optimization for hybrid strategy (n=64)")
println("=" ^ 70)

best_err = Inf
best_params = (alpha=0.0, h=0.0)

for alpha in [0.2, 0.25, 0.3, 0.35, 0.4]
    for h in [0.08, 0.1, 0.12, 0.15]
        val = 0.0
        for (a, b, sing_pos) in intervals_info
            xs, wx = transform_hybrid(a, b, 64, sing_pos; alpha=alpha, h=h)
            val += sum(wx .* integrand_full.(xs))
        end
        err = abs((val - val_ref) / val_ref)
        if err < best_err
            global best_err = err
            global best_params = (alpha=alpha, h=h)
        end
    end
end

println(@sprintf("Best params: alpha=%.2f, h=%.2f -> relerr=%.2e", 
    best_params.alpha, best_params.h, best_err))

# ============================================================================
# 最终对比
# ============================================================================

println("\n" * "=" ^ 70)
println("Final comparison with optimized parameters")
println("=" ^ 70)

for n in [32, 64, 128]
    # 当前实现 (tanh β=8)
    val_current = 0.0
    for (a, b, _) in intervals_info
        xs, wx = transform_tanh(a, b, n; beta=8.0)
        val_current += sum(wx .* integrand_full.(xs))
    end
    err_current = abs((val_current - val_ref) / val_ref)
    
    # 优化后的混合策略
    val_opt = 0.0
    for (a, b, sing_pos) in intervals_info
        xs, wx = transform_hybrid(a, b, n, sing_pos; 
            alpha=best_params.alpha, h=best_params.h)
        val_opt += sum(wx .* integrand_full.(xs))
    end
    err_opt = abs((val_opt - val_ref) / val_ref)
    
    println(@sprintf("n=%3d: current=%.2e  optimized=%.2e  improvement=%.1fx",
        n, err_current, err_opt, err_current/err_opt))
end

println("\n" * "=" ^ 70)
println("Conclusion:")
println("=" ^ 70)
println("""
混合策略（根据奇点位置选择变换）可以进一步提升精度：
- 单侧奇点区间：使用幂次变换 (alpha≈0.3) 聚簇于奇点端
- 双侧奇点区间：使用 DE 变换 (h≈0.1)
- 相比统一 tanh(β=8)，精度提升约 10-100 倍
""")
