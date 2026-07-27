using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

function _s_shape_curve()
    rho = Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    mu = Float64[250.0, 255.0, 260.0, 265.0, 268.0, 266.0, 262.0, 258.0, 260.0, 265.0, 272.0, 280.0, 288.0, 296.0, 305.0]
    return mu, rho
end

function _monotonic_curve()
    rho = Float64[0.1, 0.4, 0.7, 1.0, 1.3, 1.6]
    mu = Float64[250.0, 260.0, 270.0, 280.0, 290.0, 300.0]
    return mu, rho
end

@testset "Phase core algorithms smoke" begin
    mu_s, rho_s = _s_shape_curve()
    sres = Models.detect_s_shape(mu_s, rho_s)
    @test sres.has_s_shape

    mres = Models.maxwell_construction(mu_s, rho_s; min_samples=8)
    @test mres.converged
    @test mres.mu_transition !== nothing
    @test mres.rho_hadron !== nothing
    @test mres.rho_quark !== nothing

    mu_m, rho_m = _monotonic_curve()
    mres2 = Models.maxwell_construction(mu_m, rho_m; min_samples=6)
    @test !mres2.converged

    curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
        140.0 => (mu_s, rho_s),
        160.0 => (mu_m, rho_m),
    )
    cep = Models.find_cep(curves; tol=0.5, max_bisect_iter=8)
    @test !cep.found
    @test cep.result_status == :ambiguous
    @test cep.T_last_first_order_MeV == 140.0

    curves_nonuniform = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
        140.0 => (mu_s, rho_s),
        160.0 => (mu_s[1:end-1], rho_s[1:end-1]),
    )
    cep2 = Models.find_cep(curves_nonuniform; tol=0.5, max_bisect_iter=8)
    @test !cep2.found
    @test cep2.result_status == :ambiguous

    evaluator = function (T_mid::Float64, level::Int)
        if T_mid <= 150.0
            return mu_s, rho_s
        end
        rho_refined = if level <= 0
            rho_m
        else
            collect(range(first(rho_m), stop=last(rho_m), length=length(rho_m) + level * 2))
        end
        mu_refined = collect(range(first(mu_m), stop=last(mu_m), length=length(rho_refined)))
        return mu_refined, rho_refined
    end

    cep3 = Models.find_cep(curves_nonuniform;
        strategy=:direct,
        evaluate_at_T=evaluator,
        tol=0.5,
        max_bisect_iter=8,
        direct_bracket_mode=:directional,
        direct_start=:mid,
        direct_initial_step=5.0,
        direct_expand_factor=2.0,
        direct_max_expand_steps=6,
        direct_fallback_scan=false,
        max_refine_level=2)
    @test !cep3.found
    @test cep3.result_status == :ambiguous
    @test cep3.eval_count > 0
end
