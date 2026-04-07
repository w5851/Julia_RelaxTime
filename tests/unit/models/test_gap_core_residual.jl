using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "gap_core_residual parity" begin
    model = Models.create_model(:PNJL)
    x_state = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
    mu_vec = SVector{3}(0.0, 0.0, 0.0)
    thermal_nodes = Models.cached_nodes(8, 4)
    params = Models.GapParams(0.5, thermal_nodes, 0.0; p_num=8, t_num=4, model_kind=:PNJL)

    F_core = zeros(5)
    Models.gap_core_residual!(F_core, x_state, mu_vec, params)
    F_gap = Models.gap_conditions(x_state, mu_vec, params)

    @test length(F_core) == 5
    for i in 1:5
        @test isapprox(F_core[i], F_gap[i]; rtol=1e-9, atol=1e-9)
    end
end

@testset "build_conditions fixedmu enforces strict x length" begin
    thermal_nodes = Models.cached_nodes(8, 4)
    params = Models.GapParams(0.5, thermal_nodes, 0.0; p_num=8, t_num=4, model_kind=:PNJL)
    cond = Models.build_conditions(Models.FixedMu(), params)
    theta = [0.5, 0.0]

    @test_throws ArgumentError cond(theta, [-1.84, -1.84, -2.23, 0.5, 0.5, 0.0])
end

@testset "build_residual delegates core block" begin
    x_state = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
    mu_vec = SVector{3}(0.0, 0.0, 0.0)
    thermal_nodes = Models.cached_nodes(8, 4)
    params = Models.GapParams(0.5, thermal_nodes, 0.0; p_num=8, t_num=4, model_kind=:PNJL)

    F_core = zeros(5)
    Models.gap_core_residual!(F_core, x_state, mu_vec, params)

    fixed_mu_residual! = Models.build_residual!(Models.FixedMu(), mu_vec, params)
    F_mu = zeros(5)
    fixed_mu_residual!(F_mu, collect(x_state))
    for i in 1:5
        @test isapprox(F_mu[i], F_core[i]; rtol=1e-9, atol=1e-9)
    end

    x_full = vcat(collect(x_state), collect(mu_vec))

    fixed_rho_residual! = Models.build_residual!(Models.FixedRho(0.1), params)
    F_rho = zeros(8)
    fixed_rho_residual!(F_rho, copy(x_full))
    for i in 1:5
        @test isapprox(F_rho[i], F_core[i]; rtol=1e-9, atol=1e-9)
    end

    fixed_entropy_residual! = Models.build_residual!(Models.FixedEntropy(0.1), params)
    F_entropy = zeros(8)
    fixed_entropy_residual!(F_entropy, copy(x_full))
    for i in 1:5
        @test isapprox(F_entropy[i], F_core[i]; rtol=1e-9, atol=1e-9)
    end

    fixed_sigma_residual! = Models.build_residual!(Models.FixedSigma(1.0), params)
    F_sigma = zeros(8)
    fixed_sigma_residual!(F_sigma, copy(x_full))
    for i in 1:5
        @test isapprox(F_sigma[i], F_core[i]; rtol=1e-9, atol=1e-9)
    end

    fixed_asym_residual! = Models.build_residual!(Models.FixedAsymmetricRho(0.1, 1.0, 0.0), params)
    F_asym = zeros(8)
    fixed_asym_residual!(F_asym, copy(x_full))
    for i in 1:5
        @test isapprox(F_asym[i], F_core[i]; rtol=1e-9, atol=1e-9)
    end
end
