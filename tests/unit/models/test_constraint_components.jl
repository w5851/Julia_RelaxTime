using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "constraint components contract" begin
    @test isdefined(Models, :AbstractConstraintComponent)
    @test isdefined(Models, :constraint_dim)
    @test isdefined(Models, :constraint_name)
    @test isdefined(Models, :build_constraint_components)
    @test isdefined(Models, :constraint_total_dim)

    fixedrho_components = Models.build_constraint_components(Models.FixedRho(0.6))
    fixedasym_components = Models.build_constraint_components(Models.FixedAsymmetricRho(0.6, 1.0, 0.0))

    @test !isempty(fixedrho_components)
    @test !isempty(fixedasym_components)

    @test first(fixedrho_components) isa Models.AbstractConstraintComponent
    @test first(fixedasym_components) isa Models.AbstractConstraintComponent

    names_fixedrho = Models.constraint_name.(fixedrho_components)
    names_fixedasym = Models.constraint_name.(fixedasym_components)

    @test :stationarity in names_fixedrho
    @test :equal_mu in names_fixedrho
    @test :fixed_baryon_density in names_fixedrho

    @test :stationarity in names_fixedasym
    @test :asymmetric_density in names_fixedasym

    @test Models.constraint_total_dim(fixedrho_components) == Models.state_dim(Models.FixedRho(0.6))
    @test Models.constraint_total_dim(fixedasym_components) == Models.state_dim(Models.FixedAsymmetricRho(0.6, 1.0, 0.0))
end
