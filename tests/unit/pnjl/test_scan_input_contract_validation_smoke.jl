using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.legacy_pnjl_module()
const PNJL = Models.legacy_pnjl_module()

@testset "Scan input contract validation smoke" begin
    @test_throws ArgumentError run_tmu_scan(
        T_values=1.0,
        mu_values=[10.0],
        xi_values=[0.0],
        output_path=joinpath(mktempdir(), "x.csv"),
        overwrite=true,
        resume=false,
    )

    @test_throws ArgumentError run_tmu_scan(
        T_values=[90.0],
        mu_values=[10.0],
        xi_values=Float64[],
        output_path=joinpath(mktempdir(), "x.csv"),
        overwrite=true,
        resume=false,
    )

    @test_throws ArgumentError run_trho_scan(
        T_values=[90.0],
        rho_values=[0.2],
        xi_values=[0.0],
        seed_policy=:invalid_seed,
        output_path=joinpath(mktempdir(), "x.csv"),
        overwrite=true,
        resume=false,
    )

    @test_throws ArgumentError run_trho_scan(
        T_values=[90.0],
        rho_values=[0.2],
        xi_values=[NaN],
        output_path=joinpath(mktempdir(), "x.csv"),
        overwrite=true,
        resume=false,
    )
end
