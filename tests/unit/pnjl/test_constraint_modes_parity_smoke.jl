using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

@testset "ConstraintModes parity smoke" begin
    p_modes = (
        PNJL.FixedMu(),
        PNJL.FixedRho(0.6),
        PNJL.FixedAsymmetricRho(0.8, 1.2, 0.0),
        PNJL.FixedEntropy(0.4),
        PNJL.FixedSigma(8.0),
    )

    m_modes = (
        Models.FixedMu(),
        Models.FixedRho(0.6),
        Models.FixedAsymmetricRho(0.8, 1.2, 0.0),
        Models.FixedEntropy(0.4),
        Models.FixedSigma(8.0),
    )

    @test length(p_modes) == length(m_modes)

    for i in eachindex(p_modes)
        pm = p_modes[i]
        mm = m_modes[i]

        @test PNJL.state_dim(pm) == Models.state_dim(mm)
        @test PNJL.param_dim(pm) == Models.param_dim(mm)
        @test PNJL.constraint_description(pm) == Models.constraint_description(mm)
        @test sprint(show, pm) == sprint(show, mm)
    end
end
