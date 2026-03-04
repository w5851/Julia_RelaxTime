# Models legacy adapter smoke test (NJL)
#
# Goal: ensure NJL also has a legacy/models compare loop through the same entrypoints.
# The legacy adapter is currently a thin wrapper (see src/models/legacy/LegacyNJLModel.jl),
# but it locks in the API surface: solve_gap -> omega -> number_densities.

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@testset "Models legacy adapter (NJL)" begin
    # Keep this cheap.
    T_fm = 0.15
    mu = 0.0
    p_num = 24
    t_num = 6

    for kind in (:NJL, :LegacyNJL)
        model = Main.Models.create_model(kind)

        st = Main.Models.solve_gap(model, T_fm, mu; p_num=p_num, t_num=t_num, xi=0.0, residual_norm_max=1e-5)
        @test st isa Main.Models.MeanFieldState

        ω = Main.Models.omega(model, st, T_fm, mu; p_num=p_num, t_num=t_num, xi=0.0)
        @test isfinite(ω)

        dens = Main.Models.number_densities(model, st, T_fm, mu; p_num=p_num, t_num=t_num, xi=0.0)
        @test haskey(dens, :quark)
        @test haskey(dens, :antiquark)
        @test length(dens.quark) == 3
        @test length(dens.antiquark) == 3
        @test all(isfinite, dens.quark)
        @test all(isfinite, dens.antiquark)
    end
end
