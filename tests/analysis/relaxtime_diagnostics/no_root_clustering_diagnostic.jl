# 分析脚本：测试“无根”情况下不同聚簇策略的积分误差对比
# 注意：这是诊断/分析用脚本，不属于单元测试套件。

using Printf
using QuadGK: quadgk

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(_ROOT, "src", "relaxtime", "OneLoopIntegralsAniso.jl"))
include(joinpath(_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

λ = -0.5
k = 0.05
m = 0.3
m_prime = 0.3
ξ = -0.1
T = 0.15
μ = 0.0
Φ = 0.0
Φbar = 0.0

Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)

integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
    :quark, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar,
)

ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)

println("无根情况下不同聚簇方法对比 (n=32):")
println("-" ^ 50)

n = 32

# 标准 GL
nodes, weights = OneLoopIntegralsCorrection.GaussLegendre.gauleg(Emin, Emax, n)
val_gl = sum(weights .* integrand.(nodes))
err_gl = abs((val_gl - ref) / ref)
println(@sprintf("标准 GL:        relerr=%.2e", err_gl))

# tanh 聚簇 (对称)
xs, wx = OneLoopIntegralsCorrection.clustered_gl_nodes(Emin, Emax, n; beta=8.0)
val_tanh = sum(wx .* integrand.(xs))
err_tanh = abs((val_tanh - ref) / ref)
println(@sprintf("tanh(β=8):      relerr=%.2e", err_tanh))

# 幂次聚簇于左端（被积函数变化最剧烈处）
xs_left, wx_left = OneLoopIntegralsCorrection.power_left_nodes(Emin, Emax, n; alpha=0.35)
val_left = sum(wx_left .* integrand.(xs_left))
err_left = abs((val_left - ref) / ref)
println(@sprintf("power_left:     relerr=%.2e", err_left))

# DE 变换
xs_de, wx_de = OneLoopIntegralsCorrection.de_nodes(Emin, Emax, n; h=0.15)
val_de = sum(wx_de .* integrand.(xs_de))
err_de = abs((val_de - ref) / ref)
println(@sprintf("DE(h=0.15):     relerr=%.2e", err_de))

println("\n结论：无根情况下，power_left 聚簇于变化剧烈的左端通常更好")
