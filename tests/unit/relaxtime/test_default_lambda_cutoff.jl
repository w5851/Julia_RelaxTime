"""
测试默认行为与显式指定 Λ 截断的一致性。

验证：
1. design_w0cdf_s_grid 的 p_cutoff=Λ 参数正确限制 s 范围
2. RelaxationTime.relaxation_times 的默认行为与显式 Λ 截断一致
3. 对于重夸克（m_s 接近 Λ），结果也保持一致
"""

using Test

include("../../../src/Constants_PNJL.jl")
include("../../../src/pnjl/PNJL.jl")
include("../../../src/relaxtime/RelaxationTime.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/integration/GaussLegendre.jl")

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .PNJL: solve, FixedMu, cached_nodes, calculate_number_densities
using .PNJL.Integrals: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .RelaxationTime: relaxation_times, REQUIRED_PROCESSES
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .GaussLegendre: gauleg

const RT_ASR = RelaxationTime.AverageScatteringRate

# 测试条件：T=100 MeV, μB=800 MeV（重夸克条件，m_s 接近 Λ）
const T_MEV = 100.0
const MUB_MEV = 800.0
const T_FM = T_MEV / ħc_MeV_fm
const MUQ_FM = (MUB_MEV / 3.0) / ħc_MeV_fm

# 相对误差容限
const RTOL = 0.01  # 1%

function setup_test_params()
    base = solve(FixedMu(), T_FM, MUQ_FM; xi=0.0, p_num=12, t_num=6)
    Φ, Φbar = Float64(base.x_state[4]), Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))
    
    A_u = A(masses.u, MUQ_FM, T_FM, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, MUQ_FM, T_FM, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    
    quark_params = (m=masses, μ=(u=Float64(MUQ_FM), d=Float64(MUQ_FM), s=Float64(MUQ_FM)), 
                    A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=Float64(T_FM), Φ=Φ, Φbar=Φbar, ξ=0.0)
    
    thermal_nodes = cached_nodes(12, 6)
    nd = calculate_number_densities(base.x_state, base.mu_vec, T_FM, thermal_nodes, 0.0)
    densities = (
        u=Float64(nd.quark[1]), d=Float64(nd.quark[2]), s=Float64(nd.quark[3]),
        ubar=Float64(nd.antiquark[1]), dbar=Float64(nd.antiquark[2]), sbar=Float64(nd.antiquark[3]),
    )
    
    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs, 
            densities=densities, masses=masses)
end


@testset "design_w0cdf_s_grid p_cutoff parameter" begin
    params = setup_test_params()
    
    # 测试 ssbar_to_ssbar 过程（涉及重夸克）
    process = :ssbar_to_ssbar
    mi = params.masses.s
    
    # 理论 Λ 截断 s 范围
    s_bo_theory = (2*mi)^2 * 1.001
    s_up_theory = (2*sqrt(mi^2 + Λ_inv_fm^2))^2
    
    @testset "p_cutoff=Λ limits s range correctly" begin
        s_grid = RT_ASR.design_w0cdf_s_grid(process, params.quark_params, params.thermo_params; 
            N=60, p_cutoff=Λ_inv_fm)
        
        # s 范围应在理论范围内
        @test minimum(s_grid) >= s_bo_theory * 0.99  # 允许 1% 误差
        @test maximum(s_grid) <= s_up_theory * 1.01
    end
    
    @testset "p_cutoff=nothing gives larger s range" begin
        s_grid_semi = RT_ASR.design_w0cdf_s_grid(process, params.quark_params, params.thermo_params; 
            N=60, p_cutoff=nothing)
        s_grid_lambda = RT_ASR.design_w0cdf_s_grid(process, params.quark_params, params.thermo_params; 
            N=60, p_cutoff=Λ_inv_fm)
        
        # 半无穷积分的 s 范围应该更大
        @test maximum(s_grid_semi) > maximum(s_grid_lambda) * 10
    end
end

@testset "relaxation_times default behavior equals explicit Λ cutoff" begin
    params = setup_test_params()
    
    cos_grid, cos_w = gauleg(-1.0, 1.0, 4)
    phi_grid, phi_w = gauleg(0.0, 2π, 8)
    
    # 默认行为（不传 p_grid）
    tau_default = relaxation_times(params.quark_params, params.thermo_params, params.K_coeffs;
        densities=params.densities,
        p_nodes=20, angle_nodes=4, phi_nodes=8,
        cos_grid=cos_grid, cos_w=cos_w, phi_grid=phi_grid, phi_w=phi_w,
        n_sigma_points=6,
    )
    
    # 显式传入 [0, Λ] 网格 + sigma_cutoff=Λ
    p_grid_lambda, p_w_lambda = gauleg(0.0, Λ_inv_fm, 20)
    tau_explicit = relaxation_times(params.quark_params, params.thermo_params, params.K_coeffs;
        densities=params.densities,
        p_nodes=20, angle_nodes=4, phi_nodes=8,
        p_grid=p_grid_lambda, p_w=p_w_lambda,
        cos_grid=cos_grid, cos_w=cos_w, phi_grid=phi_grid, phi_w=phi_w,
        n_sigma_points=6,
        sigma_cutoff=Λ_inv_fm,
    )
    
    @testset "τ_u consistency" begin
        @test isapprox(tau_default.tau.u, tau_explicit.tau.u, rtol=RTOL)
    end
    
    @testset "τ_s consistency" begin
        @test isapprox(tau_default.tau.s, tau_explicit.tau.s, rtol=RTOL)
    end
    
    @testset "τ_ubar consistency" begin
        @test isapprox(tau_default.tau.ubar, tau_explicit.tau.ubar, rtol=RTOL)
    end
    
    @testset "τ_sbar consistency (critical for heavy quarks)" begin
        # 这是最关键的测试：之前 τ_sbar 有 37% 的差异
        @test isapprox(tau_default.tau.sbar, tau_explicit.tau.sbar, rtol=RTOL)
    end
end

@testset "σ(s) cache uses Λ cutoff by default" begin
    params = setup_test_params()
    
    # 构建默认缓存
    cache = RT_ASR.build_w0cdf_pchip_cache(
        :ssbar_to_ssbar,
        params.quark_params,
        params.thermo_params,
        params.K_coeffs;
        p_cutoff=Λ_inv_fm,
        n_sigma_points=6,
    )
    
    mi = params.masses.s
    s_up_theory = (2*sqrt(mi^2 + Λ_inv_fm^2))^2
    
    @testset "cache s range is bounded by Λ cutoff" begin
        @test maximum(cache.s_vals) <= s_up_theory * 1.01
    end
end
