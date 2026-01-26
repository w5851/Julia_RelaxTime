"""
详细对比 Fortran 和 Julia 的中间量

目标：找出 3 倍差异的真正来源
"""

using Printf

# 添加项目路径
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))

# 加载模块
include(joinpath(@__DIR__, "..", "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "ParameterTypes.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pnjl", "PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "relaxtime", "RelaxationTime.jl"))

using .Constants_PNJL
using .ParameterTypes
using .PNJL
using .RelaxationTime

println("="^70)
println("Fortran vs Julia 详细对比")
println("="^70)

# 参数设置（与 Fortran 一致）
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / Constants_PNJL.ħc_MeV_fm
μ = μ_MeV / Constants_PNJL.ħc_MeV_fm

println("\n参数设置:")
println("  T = $T_MeV MeV = $(round(T, digits=6)) fm⁻¹")
println("  μ = $μ_MeV MeV = $(round(μ, digits=6)) fm⁻¹")
println("  Λ = $(round(Constants_PNJL.Λ_inv_fm, digits=3)) fm⁻¹")

# 求解 PNJL 平衡态
println("\n" * "="^70)
println("步骤 1: 求解 PNJL 平衡态")
println("="^70)

result = PNJL.solve_pnjl_equilibrium(T, μ)
Φ = result.Phi
Φbar = result.Phibar
m_u = result.m_u
m_d = result.m_d
m_s = result.m_s

println("\nPolyakov loops:")
println("  Φ     = $(round(Φ, digits=6))")
println("  Φbar  = $(round(Φbar, digits=6))")

println("\n有效夸克质量 (fm⁻¹):")
println("  m_u   = $(round(m_u, digits=6)) fm⁻¹ = $(round(m_u * Constants_PNJL.ħc_MeV_fm, digits=3)) MeV")
println("  m_d   = $(round(m_d, digits=6)) fm⁻¹ = $(round(m_d * Constants_PNJL.ħc_MeV_fm, digits=3)) MeV")
println("  m_s   = $(round(m_s, digits=6)) fm⁻¹ = $(round(m_s * Constants_PNJL.ħc_MeV_fm, digits=3)) MeV")

# 构建参数
quark_params = QuarkParams(
    m = (u=m_u, d=m_d, s=m_s),
    μ = (u=μ, d=μ, s=μ)
)
thermo_params = ThermoParams(T, Φ, Φbar, 0.0)
K_coeffs = Constants_PNJL.K_COEFFS

# Fortran 结果（半无穷积分）
println("\n" * "="^70)
println("步骤 2: Fortran 结果（半无穷积分）")
println("="^70)

fortran_tau_u = 0.58097
fortran_tau_s = 0.59204
fortran_ratio = fortran_tau_u / fortran_tau_s

println("\nFortran 弛豫时间:")
println("  τ_u   = $(round(fortran_tau_u, digits=5)) fm")
println("  τ_s   = $(round(fortran_tau_s, digits=5)) fm")
println("  τ_u/τ_s = $(round(fortran_ratio, digits=4))")

# Julia 结果
println("\n" * "="^70)
println("步骤 3: Julia 结果（半无穷积分 + Λ截断σ(s)）")
println("="^70)

# 计算弛豫时间
tau_result = compute_relaxation_times(
    quark_params,
    thermo_params,
    K_coeffs;
    p_nodes=20,
    angle_nodes=4,
    phi_nodes=8,
    n_sigma_points=32,
    sigma_cutoff=Constants_PNJL.Λ_inv_fm,  # 使用 Λ 截断
)

julia_tau_u = tau_result.tau_u
julia_tau_s = tau_result.tau_s
julia_ratio = julia_tau_u / julia_tau_s

println("\nJulia 弛豫时间:")
println("  τ_u   = $(round(julia_tau_u, digits=5)) fm")
println("  τ_s   = $(round(julia_tau_s, digits=5)) fm")
println("  τ_u/τ_s = $(round(julia_ratio, digits=4))")

# 对比
println("\n" * "="^70)
println("步骤 4: 对比分析")
println("="^70)

ratio_tau_u = julia_tau_u / fortran_tau_u
ratio_tau_s = julia_tau_s / fortran_tau_s
ratio_ratio = julia_ratio / fortran_ratio

println("\n绝对值对比:")
println("  τ_u (Julia/Fortran)   = $(round(ratio_tau_u, digits=3))")
println("  τ_s (Julia/Fortran)   = $(round(ratio_tau_s, digits=3))")

println("\n比值对比:")
println("  τ_u/τ_s (Julia)       = $(round(julia_ratio, digits=4))")
println("  τ_u/τ_s (Fortran)     = $(round(fortran_ratio, digits=4))")
println("  比值的比值            = $(round(ratio_ratio, digits=4))")

# 详细分析：对比单个散射过程
println("\n" * "="^70)
println("步骤 5: 单个散射过程对比")
println("="^70)

# 选择几个关键过程
processes = [:uu_to_uu, :us_to_us, :ss_to_ss, :ssbar_to_uubar, :uubar_to_ssbar]

println("\n计算平均散射率...")
for process in processes
    rate = RelaxationTime.AverageScatteringRate.average_scattering_rate(
        process,
        quark_params,
        thermo_params,
        K_coeffs;
        p_nodes=20,
        angle_nodes=4,
        phi_nodes=8,
        n_sigma_points=32,
        sigma_cutoff=Constants_PNJL.Λ_inv_fm,
    )
    println("  ⟨Γ⟩($process) = $(round(rate, sigdigits=6)) fm⁻¹")
end

# 对比 σ(s) 的值
println("\n" * "="^70)
println("步骤 6: 截面 σ(s) 对比")
println("="^70)

# 选择一个典型的 s 值
s_test_MeV2 = 500000.0  # 500 GeV²
s_test = s_test_MeV2 / Constants_PNJL.ħc_MeV_fm^2

println("\n测试点: s = $s_test_MeV2 MeV² = $(round(s_test, digits=3)) fm⁻²")

for process in [:uu_to_uu, :us_to_us, :ss_to_ss, :ssbar_to_uubar]
    σ = RelaxationTime.TotalCrossSection.total_cross_section(
        process,
        s_test,
        quark_params,
        thermo_params,
        K_coeffs;
        n_points=32
    )
    println("  σ($process) = $(round(σ, sigdigits=6)) fm²")
end

# 检查数密度
println("\n" * "="^70)
println("步骤 7: 数密度对比")
println("="^70)

# 计算数密度（使用半无穷积分）
ρ_u = RelaxationTime.AverageScatteringRate.number_density(
    :u, m_u, μ, T, Φ, Φbar, 0.0;
    p_nodes=32,
    angle_nodes=2,
    scale=10.0
)

ρ_s = RelaxationTime.AverageScatteringRate.number_density(
    :s, m_s, μ, T, Φ, Φbar, 0.0;
    p_nodes=32,
    angle_nodes=2,
    scale=10.0
)

println("\nJulia 数密度 (半无穷积分):")
println("  ρ_u = $(round(ρ_u, digits=6)) fm⁻³")
println("  ρ_s = $(round(ρ_s, digits=6)) fm⁻³")

println("\nFortran 数密度:")
println("  ρ_u = 1.680901 fm⁻³")
println("  ρ_s = 1.537628 fm⁻³")

println("\n比值:")
println("  ρ_u (Julia/Fortran) = $(round(ρ_u / 1.680901, digits=4))")
println("  ρ_s (Julia/Fortran) = $(round(ρ_s / 1.537628, digits=4))")

# 总结
println("\n" * "="^70)
println("总结")
println("="^70)

println("\n1. 弛豫时间差异:")
println("   - Julia 的 τ 是 Fortran 的 $(round(ratio_tau_u, digits=2))-$(round(ratio_tau_s, digits=2)) 倍")
println("   - 绝对值差异很大")

println("\n2. 比值差异:")
println("   - Julia: τ_u/τ_s = $(round(julia_ratio, digits=3))")
println("   - Fortran: τ_u/τ_s = $(round(fortran_ratio, digits=3))")
println("   - 差异约 $(round((1 - ratio_ratio) * 100, digits=1))%")

println("\n3. 可能的差异来源:")
println("   - 数密度计算？")
println("   - 散射率计算？")
println("   - σ(s) 计算？")
println("   - 积分精度？")

println("\n" * "="^70)
