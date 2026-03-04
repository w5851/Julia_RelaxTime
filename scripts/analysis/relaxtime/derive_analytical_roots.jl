# 推导 A±B=0 根的解析解
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))

println("=" ^ 70)
println("推导 A±B=0 根的解析解")
println("=" ^ 70)

println("""
方程设置:
  A = coeff_x = 2*k*p = 2*k*sqrt(E² - m²)
  B = denominator_const = λ² + 2λE + m² - m'² - k²

求解 A ± B = 0:
  2*k*sqrt(E² - m²) = ∓(λ² + 2λE + m² - m'² - k²)

令 C = λ² + m² - m'² - k², 则:
  2*k*sqrt(E² - m²) = ∓(C + 2λE)

两边平方:
  4k²(E² - m²) = (C + 2λE)²
  4k²E² - 4k²m² = C² + 4λCE + 4λ²E²
  (4k² - 4λ²)E² - 4λCE - (4k²m² + C²) = 0
  (k² - λ²)E² - λCE - (k²m² + C²/4) = 0

这是关于 E 的二次方程！
  aE² + bE + c = 0
  a = k² - λ²
  b = -λC = -λ(λ² + m² - m'² - k²)
  c = -(k²m² + C²/4)

判别式:
  Δ = b² - 4ac
""")

# 解析求根函数
function analytical_roots(λ::Float64, k::Float64, m::Float64, m_prime::Float64, Emin::Float64, Emax::Float64)
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
            return [E]
        else
            return Float64[]
        end
    end
    
    Δ = b^2 - 4*a*c
    
    if Δ < 0
        return Float64[]  # 无实根
    end
    
    sqrt_Δ = sqrt(Δ)
    E1 = (-b - sqrt_Δ) / (2*a)
    E2 = (-b + sqrt_Δ) / (2*a)
    
    # 筛选有效根：
    # 1. E >= m (物理约束)
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

# 测试用例
test_cases = [
    (name="有根1", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根2", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="有根3", λ=-0.8, k=0.02, m=0.25, m_prime=0.25),
    (name="无根1", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (name="无根2", λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

println("\n验证解析解 vs 数值解:")
println("-" ^ 70)
println(@sprintf("%-8s %20s %20s %12s", "Case", "数值解", "解析解", "匹配?"))
println("-" ^ 70)

for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    
    # 数值解
    roots_num = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    
    # 解析解
    roots_ana = analytical_roots(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    
    # 比较
    match = length(roots_num) == length(roots_ana)
    if match && !isempty(roots_num)
        for (rn, ra) in zip(sort(roots_num), sort(roots_ana))
            if abs(rn - ra) > 1e-6
                match = false
                break
            end
        end
    end
    
    num_str = isempty(roots_num) ? "[]" : string(round.(roots_num, digits=6))
    ana_str = isempty(roots_ana) ? "[]" : string(round.(roots_ana, digits=6))
    
    println(@sprintf("%-8s %20s %20s %12s", tc.name, num_str, ana_str, match ? "✓" : "✗"))
end

# 性能对比
println("\n" * "=" ^ 70)
println("性能对比 (10000次):")
println("=" ^ 70)

tc = test_cases[1]
Emin = tc.m
Emax = OneLoopIntegrals.energy_cutoff(tc.m)

# 预热
for _ in 1:100
    OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    analytical_roots(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
end

n_runs = 10000

t0 = time_ns()
for _ in 1:n_runs
    OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
end
time_num = (time_ns() - t0) / n_runs / 1e6

t0 = time_ns()
for _ in 1:n_runs
    analytical_roots(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
end
time_ana = (time_ns() - t0) / n_runs / 1e6

println(@sprintf("数值求根 (N=2000): %.4f ms", time_num))
println(@sprintf("解析求根:          %.4f ms", time_ana))
println(@sprintf("加速比:            %.1fx", time_num / time_ana))
