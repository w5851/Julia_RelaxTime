# Models phase-0 contract smoke test
#
# Validates the new unified entrypoints:
# - Models.solve_gap(model, T, mu_vec)
# - Models.omega(model, x_state, T, mu_vec)
# - stable state + mu_vec normalization helpers

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Load models entry into Main to match how core bridge modules use it.
const _MODELS_ENTRY = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _MODELS_ENTRY)
end

@testset "Models phase 0 contract" begin
    model = Main.Models.create_model(:PNJL)

    # Keep this intentionally cheap; solver smoke already covers robustness.
    T_fm = 0.15
    mu = 0.0

    x_state = Main.Models.solve_gap(model, T_fm, mu; p_num=12, t_num=4, xi=0.0, residual_norm_max=1e-5)
    @test x_state isa Main.Models.MeanFieldState

    # Accept scalar μ by contract (treated as (μ, μ, μ)).
    ω = Main.Models.omega(model, x_state, T_fm, mu; p_num=12, t_num=4, xi=0.0)
    @test isfinite(ω)

    comps = Main.Models.omega_components(model, x_state, T_fm, mu; p_num=12, t_num=4, xi=0.0)
    @test all(isfinite, (comps.chi, comps.poly, comps.vac, comps.therm, comps.omega))

    # Contract helpers
    v = Main.Models.state_vector(x_state)
    @test length(v) == 5
    μv = Main.Models.normalize_mu_vec(mu)
    @test length(μv) == 3
end
