"""
测试扫描模块迁移

验证新架构的 TmuScan 和 TrhoScan 是否正常工作。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# 加载模块
include(joinpath(@__DIR__, "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "src", "pnjl", "PNJL.jl"))

using .PNJL
using .PNJL: FixedMu, FixedRho, solve, SolverResult
using .PNJL: run_tmu_scan, run_trho_scan

println("=" ^ 60)
println("测试扫描模块迁移")
println("=" ^ 60)

# 测试 1: 单点求解 (FixedMu)
println("\n[测试 1] 单点求解 (FixedMu)")
T_fm = 100.0 / 197.327  # 100 MeV
μ_fm = 200.0 / 197.327  # 200 MeV

try
    result = solve(FixedMu(), T_fm, μ_fm; xi=0.0)
    println("  收敛: $(result.converged)")
    println("  迭代: $(result.iterations)")
    println("  残差: $(result.residual_norm)")
    println("  压强: $(result.pressure) fm⁻⁴")
    println("  密度: $(result.rho_norm) ρ₀")
    println("  ✓ FixedMu 求解成功")
catch e
    println("  ✗ FixedMu 求解失败: $e")
end

# 测试 2: 单点求解 (FixedRho)
println("\n[测试 2] 单点求解 (FixedRho)")
T_fm = 100.0 / 197.327  # 100 MeV
rho_target = 1.0  # ρ₀

try
    result = solve(FixedRho(rho_target), T_fm; xi=0.0)
    println("  收敛: $(result.converged)")
    println("  迭代: $(result.iterations)")
    println("  残差: $(result.residual_norm)")
    println("  化学势: $(result.mu_vec[1] * 197.327) MeV")
    println("  密度: $(result.rho_norm) ρ₀")
    println("  ✓ FixedRho 求解成功")
catch e
    println("  ✗ FixedRho 求解失败: $e")
end

# 测试 3: 小规模 T-μ 扫描
println("\n[测试 3] 小规模 T-μ 扫描")
output_tmu = joinpath(@__DIR__, "..", "data", "outputs", "results", "pnjl", "test_tmu_scan.csv")

try
    result = run_tmu_scan(
        T_values = [100.0, 150.0],
        mu_values = [0.0, 100.0, 200.0],
        xi_values = [0.0],
        output_path = output_tmu,
        overwrite = true
    )
    println("  总点数: $(result.total)")
    println("  成功: $(result.success)")
    println("  失败: $(result.failure)")
    println("  输出: $(result.output)")
    println("  ✓ T-μ 扫描成功")
catch e
    println("  ✗ T-μ 扫描失败: $e")
    rethrow(e)
end

# 测试 4: 小规模 T-ρ 扫描
println("\n[测试 4] 小规模 T-ρ 扫描")
output_trho = joinpath(@__DIR__, "..", "data", "outputs", "results", "pnjl", "test_trho_scan.csv")

try
    result = run_trho_scan(
        T_values = [100.0, 150.0],
        rho_values = [0.5, 1.0, 1.5],
        xi_values = [0.0],
        output_path = output_trho,
        overwrite = true
    )
    println("  总点数: $(result.total)")
    println("  成功: $(result.success)")
    println("  失败: $(result.failure)")
    println("  输出: $(result.output)")
    println("  ✓ T-ρ 扫描成功")
catch e
    println("  ✗ T-ρ 扫描失败: $e")
    rethrow(e)
end

println("\n" * "=" ^ 60)
println("测试完成")
println("=" ^ 60)
