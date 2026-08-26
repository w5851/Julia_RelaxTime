using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const MODEL_LAYOUT = Dict(
    :njl => (base=joinpath(PROJECT_ROOT, "src", "models", "njl"), core_required=true),
    :njl2 => (base=joinpath(PROJECT_ROOT, "src", "models", "njl2"), core_required=true),
    :pnjl => (base=joinpath(PROJECT_ROOT, "src", "models", "pnjl"), core_required=true),
    :rpnjl => (base=joinpath(PROJECT_ROOT, "src", "models", "rpnjl"), core_required=true),
    # Magnetic physics is owned by pnjl_physics; this directory is an adapter facade.
    :pnjl_magnetic => (base=joinpath(PROJECT_ROOT, "src", "models", "pnjl_magnetic"), core_required=false),
    :rotation => (base=joinpath(PROJECT_ROOT, "src", "models", "rotation"), core_required=true),
    :gas_liquid => (base=joinpath(PROJECT_ROOT, "src", "models", "gas_liquid"), core_required=true),
)

@testset "models structure homomorphism" begin
    for (kind, layout) in MODEL_LAYOUT
        base = layout.base
        @test isdir(base)
        @test isfile(joinpath(base, "api.jl"))
        @test isfile(joinpath(base, "capabilities.jl"))
        layout.core_required && @test isdir(joinpath(base, "core"))
        @test isdir(joinpath(base, "adapters"))
        @test isdir(joinpath(base, "workflows"))
    end
end
