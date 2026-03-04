using Test

include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "njl", "NJLCore.jl"))

@testset "Config profile selection (NJL)" begin
    withenv("PHYSICS_PARAM_PROFILE" => "default") do
        # New scheme preferred over old scheme when both exist:
        # - config/models/njl/<profile>.toml has priority over config/njl/<profile>.toml
        p_new = NJLCore.njl_params(profile="unittest")
        @test p_new.label == "njl-unittest-new"

        # Fallback to legacy path when only old config exists.
        p_old = NJLCore.njl_params(profile="unittest_oldonly")
        @test p_old.label == "njl-unittest-oldonly"
    end

    # Physics profile is selected via PHYSICS_PARAM_PROFILE.
    withenv("PHYSICS_PARAM_PROFILE" => "unittest") do
        p_phys = NJLCore.njl_params(profile="default")
        @test p_phys.hbarc_MeV_fm == 199.0
    end
end
