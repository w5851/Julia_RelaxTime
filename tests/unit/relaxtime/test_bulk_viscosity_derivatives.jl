"""
单元测试：体粘滞系数热力学导数计算

测试 ThermoDerivatives 模块中的以下函数：
- v_n_squared: 热力学速度 v_n²
- dmuB_dT_sigma: 固定 σ 时的导数 ∂μ_B/∂T|_σ
- bulk_viscosity_coefficients: 体粘滞系数所需的所有导数
- compute_B_bracket: 体粘滞公式中的 B 项
"""

using Test
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))

include("../../../src/pnjl/PNJL.jl")
using .PNJL: ThermoDerivatives, AnisoGapSolver

using StaticArrays
using ForwardDiff

# 常数
const hc = 197.327  # MeV·fm

# 测试点：T = 150 MeV, μ_B = 800 MeV
const T_MeV = 150.0
const μB_MeV = 800.0
const T_fm = T_MeV / hc
const μq_fm = μB_MeV / 3.0 / hc

# 辅助函数：提取纯数值（处理可能的 Dual 类型）
function extract_value(x)
    if x isa ForwardDiff.Dual
        return ForwardDiff.value(x)
    else
        return Float64(x)
    end
end

@testset "体粘滞系数热力学导数" begin
    
    @testset "v_n_squared 计算" begin
        v_n_sq = ThermoDerivatives.v_n_squared(T_fm, μq_fm)
        
        # v_n² 应该是正数且在合理范围内
        @test v_n_sq > 0
        @test v_n_sq < 1  # 声速平方应该小于光速
        
        # 预期值（从测试脚本得到）
        @test isapprox(v_n_sq, 0.098640, rtol=1e-4)
    end

    @testset "dmuB_dT_sigma 计算" begin
        dμB_dT_sig = ThermoDerivatives.dmuB_dT_sigma(T_fm, μq_fm)
        
        # ∂μ_B/∂T|_σ 应该是正数（温度升高时化学势也升高以保持 σ 不变）
        @test dμB_dT_sig > 0
        
        # 预期值（从测试脚本得到）
        @test isapprox(dμB_dT_sig, 9.543024, rtol=1e-4)
    end
    
    @testset "bulk_viscosity_coefficients 计算" begin
        coeffs = ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
        
        # 检查返回的 NamedTuple 包含所有必要的字段
        @test haskey(coeffs, :v_n_sq)
        @test haskey(coeffs, :dμB_dT_sigma)
        @test haskey(coeffs, :masses)
        @test haskey(coeffs, :dM_dT)
        @test haskey(coeffs, :dM_dμB)
        @test haskey(coeffs, :s)
        @test haskey(coeffs, :n_B)
        
        # 检查数值
        @test isapprox(coeffs.v_n_sq, 0.098640, rtol=1e-4)
        @test isapprox(coeffs.dμB_dT_sigma, 9.543024, rtol=1e-4)
        
        # 质量应该是正数
        @test all(coeffs.masses .> 0)
        
        # u/d 夸克质量应该相近
        @test isapprox(coeffs.masses[1], coeffs.masses[2], rtol=1e-6)
        
        # s 夸克质量应该更大
        @test coeffs.masses[3] > coeffs.masses[1]
        
        # 质量导数（对重子化学势）
        # ∂M/∂μ_B = (1/3)·∂M/∂μ_q
        @test isapprox(coeffs.dM_dμB[1], -1.053308, rtol=1e-4)
        
        # 熵密度和重子数密度应该是正数
        @test coeffs.s > 0
        @test coeffs.n_B > 0
    end

    @testset "compute_B_bracket 计算" begin
        coeffs = ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
        
        p_test = 1.0  # fm⁻¹
        flavor = 1    # u 夸克
        
        # 夸克的 B 项
        B_q = ThermoDerivatives.compute_B_bracket(
            p_test,
            coeffs.masses[flavor],
            μq_fm,
            T_fm,
            coeffs.v_n_sq,
            coeffs.dμB_dT_sigma,
            coeffs.dM_dT[flavor],
            coeffs.dM_dμB[flavor];
            is_antiquark=false
        )
        
        # 反夸克的 B 项
        B_qbar = ThermoDerivatives.compute_B_bracket(
            p_test,
            coeffs.masses[flavor],
            μq_fm,
            T_fm,
            coeffs.v_n_sq,
            coeffs.dμB_dT_sigma,
            coeffs.dM_dT[flavor],
            coeffs.dM_dμB[flavor];
            is_antiquark=true
        )
        
        # 预期值（从测试脚本得到）
        @test isapprox(B_q, -1.203043, rtol=1e-4)
        @test isapprox(B_qbar, -0.519924, rtol=1e-4)
        
        # 夸克和反夸克的 B 项应该不同（因为化学势符号不同）
        @test B_q != B_qbar
    end
    
    @testset "不同动量下的 B 项" begin
        coeffs = ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
        flavor = 1
        
        # 提取纯数值（避免 Dual 类型问题）
        M_val = extract_value(coeffs.masses[flavor])
        v_n_sq_val = extract_value(coeffs.v_n_sq)
        dμB_dT_sigma_val = extract_value(coeffs.dμB_dT_sigma)
        dM_dT_val = extract_value(coeffs.dM_dT[flavor])
        dM_dμB_val = extract_value(coeffs.dM_dμB[flavor])
        
        # 测试不同动量
        p_values = [0.5, 1.0, 2.0, 3.0]
        B_values = Float64[]
        
        for p in p_values
            B = ThermoDerivatives.compute_B_bracket(
                p, M_val, μq_fm, T_fm,
                v_n_sq_val, dμB_dT_sigma_val,
                dM_dT_val, dM_dμB_val;
                is_antiquark=false
            )
            push!(B_values, extract_value(B))
        end
        
        # B 项应该随动量变化（因为 p² 项）
        @test length(unique(B_values)) == length(B_values)
    end

    @testset "热力学恒等式验证" begin
        # 验证 s = ∂P/∂T 和 n_B = (1/3)·∂P/∂μ_q
        td = ThermoDerivatives.thermo_derivatives(T_fm, μq_fm)
        
        # s = ∂P/∂T
        @test isapprox(td.entropy, td.dP_dT, rtol=1e-6)
        
        # n_B = (1/3)·∂P/∂μ_q
        @test isapprox(td.rho, td.dP_dmu / 3.0, rtol=1e-6)
    end
    
    @testset "链式法则一致性" begin
        coeffs = ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
        
        p_test = 1.0
        M = coeffs.masses[1]
        E = sqrt(p_test^2 + M^2)
        
        # 方法1：使用 Fortran 公式
        # ∂x/∂T|_σ = ∂E/∂T + (∂E/∂μ_B - 1/3)·∂μ_B/∂T|_σ
        dE_dT = (M / E) * coeffs.dM_dT[1]
        dE_dμB = (M / E) * coeffs.dM_dμB[1]
        dx_dT_fortran = dE_dT + (dE_dμB - 1/3) * coeffs.dμB_dT_sigma
        
        # 方法2：链式法则展开
        # ∂x/∂T|_σ = ∂E/∂T + ∂E/∂μ_B·∂μ_B/∂T|_σ - ∂μ_q/∂T|_σ
        # 其中 ∂μ_q/∂T|_σ = (1/3)·∂μ_B/∂T|_σ
        dx_dT_chain = dE_dT + dE_dμB * coeffs.dμB_dT_sigma - (1/3) * coeffs.dμB_dT_sigma
        
        # 两种方法应该给出相同结果
        @test isapprox(dx_dT_fortran, dx_dT_chain, rtol=1e-10)
    end
    
    @testset "不同味道的导数" begin
        coeffs = ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
        
        # u 和 d 夸克的导数应该相同（同位旋对称）
        @test isapprox(coeffs.dM_dT[1], coeffs.dM_dT[2], rtol=1e-6)
        @test isapprox(coeffs.dM_dμB[1], coeffs.dM_dμB[2], rtol=1e-6)
        
        # s 夸克的导数应该不同
        @test coeffs.dM_dT[3] != coeffs.dM_dT[1]
    end
    
end  # @testset "体粘滞系数热力学导数"

println("\n所有单元测试通过！")
