using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const MODEL_LAYOUT = Dict(
    :njl => joinpath(PROJECT_ROOT, "src", "models", "njl"),
    :njl2 => joinpath(PROJECT_ROOT, "src", "models", "njl2"),
    :pnjl => joinpath(PROJECT_ROOT, "src", "models", "pnjl"),
    :rpnjl => joinpath(PROJECT_ROOT, "src", "models", "rpnjl"),
    :pnjl_magnetic => joinpath(PROJECT_ROOT, "src", "models", "pnjl_magnetic"),
    :rotation => joinpath(PROJECT_ROOT, "src", "models", "rotation"),
    :gas_liquid => joinpath(PROJECT_ROOT, "src", "models", "gas_liquid"),
)

@testset "models structure homomorphism" begin
    for (kind, base) in MODEL_LAYOUT
        @test isdir(base)
        @test isfile(joinpath(base, "api.jl"))
        @test isfile(joinpath(base, "capabilities.jl"))
        @test isdir(joinpath(base, "core"))
        @test isdir(joinpath(base, "adapters"))
        @test isdir(joinpath(base, "workflows"))
    end
end
