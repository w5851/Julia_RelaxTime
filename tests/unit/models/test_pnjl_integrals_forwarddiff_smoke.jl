using Test
using StaticArrays
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@testset "Models PNJLIntegrals ForwardDiff smoke" begin
    # Use cached nodes (Float64) but propagate Dual through mu/T into calculate_log_sum.
    p_mesh, cosθ_mesh, coeffs = Main.Models.PNJLIntegrals.cached_nodes(10, 2)

    masses = SVector{3, Float64}(0.30, 0.40, 0.50)
    Φ = 0.12
    Φbar = 0.34
    xi = 0.3

    f_mu(x) = Main.Models.PNJLIntegrals.calculate_log_sum(
        masses,
        p_mesh,
        cosθ_mesh,
        coeffs,
        Φ,
        Φbar,
        SVector{3, typeof(x)}(x, x, x),
        0.20,
        xi,
    )

    f_T(x) = Main.Models.PNJLIntegrals.calculate_log_sum(
        masses,
        p_mesh,
        cosθ_mesh,
        coeffs,
        Φ,
        Φbar,
        SVector{3, Float64}(0.1, 0.1, 0.1),
        x,
        xi,
    )

    dmu = ForwardDiff.derivative(f_mu, 0.10)
    dT = ForwardDiff.derivative(f_T, 0.20)

    @test isfinite(dmu)
    @test isfinite(dT)
end
