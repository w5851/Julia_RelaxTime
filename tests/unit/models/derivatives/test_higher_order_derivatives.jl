using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "higher-order derivatives utilities" begin
    f(x) = x^5

    @test Models.nth_derivative(f, 2.0, 4) ≈ 240.0 atol = 1e-10
    @test Models.nth_derivative(f, 2.0, 5) ≈ 120.0 atol = 1e-10

    @test_throws ArgumentError Models.susceptibility_scale(0.0, 2)
    @test_throws ArgumentError Models.susceptibility_scale(-1.0, 2)

    struct DummyChiProvider <: Models.AbstractChiProvider end
    @test Models.chi_provider_name(DummyChiProvider()) == :unknown
end
