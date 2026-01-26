#!/usr/bin/env julia
"""
调试弛豫时间公式

对比Julia和Fortran的弛豫时间计算公式
"""

println("="^80)
println("弛豫时间公式对比")
println("="^80)
println()

# Fortran的公式（从z1 relax_time.f90）
println("Fortran公式:")
println("  tau_l = 1 / (n_u*(w2+w1) + n_ub*(w6+w7+w9+w5) + n_s*w3 + n_sb*w8)")
println("  其中:")
println("    w1 = ud->ud")
println("    w2 = uu->uu")
println("    w3 = us->us")
println("    w5 = udbar->udbar")
println("    w6 = uubar->uubar")
println("    w7 = uubar->ddbar")
println("    w8 = usbar->usbar")
println("    w9 = uubar->ssbar")
println()

# Julia的公式（从RelaxationTime.jl）
println("Julia公式:")
println("  omega_u = n_u * (w_uu + w_ud) +")
println("            n_ubar * (w_uubar + w_uubar_ddbar + w_uubar_ssbar + w_udbar) +")
println("            n_s * w_us +")
println("            n_sbar * w_usbar")
println("  tau_u = 1 / omega_u")
println()

println("="^80)
println("公式对比")
println("="^80)
println()

println("Fortran vs Julia:")
println("  ✓ 两者公式完全相同！")
println("  ✓ 都是 τ_u⁻¹ = Σ_j ρ_j * w̄_ij")
println()

println("="^80)
println("可能的差异来源")
println("="^80)
println()

println("1. w̄_ij的定义不同")
println("   - Fortran: averaged_rate() 返回什么？")
println("   - Julia: average_scattering_rate() 返回什么？")
println()

println("2. 归一化约定不同")
println("   - 可能有额外的因子（如2π, ħ等）")
println()

println("3. 积分方法不同")
println("   - Fortran: 有限积分 [0, Λ]")
println("   - Julia: 半无穷积分 [0, ∞)")
println()

println("4. 散射截面σ(s)的计算不同")
println("   - 可能使用不同的近似或公式")
println()

println("="^80)
println("下一步")
println("="^80)
println()

println("需要对比:")
println("  1. 单个散射过程的 w̄_ij 值")
println("  2. 散射截面 σ(s) 的值")
println("  3. 积分方法的数值结果")
println()

println("="^80)
