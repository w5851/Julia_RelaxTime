# ScanResultFinalize 模块单元测试

using Test
using StaticArrays

const PROJECT_ROOT_SRF = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_SRF, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ScanResultFinalize 通过 Models 公共接口访问 SolverResult
using Main.Models.ScanResultFinalize: is_success, promote_near_converged

const PNJL_SRF = Models.pnjl_module()
const SolverResult_srf = getproperty(PNJL_SRF, :SolverResult)
const FixedMu_srf = getproperty(PNJL_SRF, :FixedMu)

@testset "ScanResultFinalize" begin
    @testset "is_success 判定" begin
        # SolverResult 位置参数: (mode, converged, solution, x_state, mu_vec, omega, pressure, rho_norm, entropy, energy, masses, iterations, residual_norm, xi)
        result = SolverResult_srf(
            FixedMu_srf(),            # mode
            true,                     # converged
            zeros(5),                 # solution
            SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2),  # x_state
            SVector{3, Float64}(0.0, 0.0, 0.0),                   # mu_vec
            0.0,                      # omega
            0.0,                      # pressure
            0.0,                      # rho_norm
            0.0,                      # entropy
            0.0,                      # energy
            SVector{3, Float64}(0.3, 0.31, 0.5),                  # masses
            5,                        # iterations
            1e-12,                    # residual_norm
            0.0                       # xi
        )
        @test is_success(result; acceptable_residual=1e-6)
    end

    @testset "未收敛结果" begin
        result = SolverResult_srf(
            FixedMu_srf(),
            false,
            zeros(5),
            SVector{5, Float64}(0.0, 0.0, 0.0, 0.0, 0.0),
            SVector{3, Float64}(0.0, 0.0, 0.0),
            0.0, 0.0, 0.0, 0.0, 0.0,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            100,
            1.0,
            0.0
        )
        @test !is_success(result; acceptable_residual=1e-6)
    end
end
