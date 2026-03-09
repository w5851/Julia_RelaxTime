using Test

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
end

@testset "PNJL config dynamic constants (no load-time cache dependency)" begin
    # Model profiles are loaded only from config/models/pnjl/<profile>.toml.
    c_new = Constants_PNJL.pnjl_constants(profile="unittest", physics_profile="default")
    expected = 333.0 / c_new.hbarc_MeV_fm
    @test isapprox(c_new.T0_inv_fm, expected; rtol=0.0, atol=0.0)

    # Legacy-only profile has been migrated into config/models/pnjl/.
    c_old = Constants_PNJL.pnjl_constants(profile="unittest_oldonly", physics_profile="default")
    expected_old = 111.0 / c_old.hbarc_MeV_fm
    @test isapprox(c_old.T0_inv_fm, expected_old; rtol=0.0, atol=0.0)

    # Physics profile is selected explicitly and should not depend on module load time.
    c_phys = Constants_PNJL.pnjl_constants(profile="default", physics_profile="unittest")
    @test c_phys.hbarc_MeV_fm == 199.0
end
