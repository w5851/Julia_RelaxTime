# Models legacy adapter smoke test
#
# Validates phase-1 direction: legacy core can be wrapped as a Models.*Model
# and used through the same unified entrypoints.

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@testset "Models legacy adapter (PNJL)" begin
    model = Main.Models.create_model(:LegacyPNJL)

    # Keep this cheap; PNJL solver robustness is covered elsewhere.
    T_fm = 0.15
    mu = 0.0

    st = Main.Models.solve_gap(model, T_fm, mu; p_num=12, t_num=4, xi=0.0, residual_norm_max=1e-5)
    @test st isa Main.Models.MeanFieldState

    ω = Main.Models.omega(model, st, T_fm, mu; p_num=12, t_num=4, xi=0.0)
    @test isfinite(ω)

    comps = Main.Models.omega_components(model, st, T_fm, mu; p_num=12, t_num=4, xi=0.0)
    @test all(isfinite, (comps.chi, comps.poly, comps.vac, comps.therm, comps.omega))
end
