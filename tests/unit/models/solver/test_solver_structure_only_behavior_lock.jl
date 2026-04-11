using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const SNAPSHOT_SOURCE = "migration-pre-main smoke sample"
const SNAPSHOT_COMMAND = "julia --project=. -e 'include(\"src/models/Models.jl\"); model=Main.Models.create_model(:PNJL); mode=Main.Models.FixedMu(); T_fm=100.0/197.327; mu_fm=0.0; r=Main.Models.solve(model, mode, T_fm, mu_fm; p_num=8, t_num=4, residual_norm_max=1e-6); println((converged=r.converged, solution_length=length(r.solution), residual_norm=r.residual_norm, omega=r.omega, pressure=r.pressure, rho_norm=r.rho_norm, mu_vec=collect(r.mu_vec), masses=collect(r.masses)))'"

const STRUCTURE_LOCK_BASELINE = (
    converged=true,
    solution_length=5,
    residual_norm=1.1087008365216668e-15,
    omega=-21.60808181242041,
    pressure=21.60808181242041,
    rho_norm=3.127689543869961e-19,
    mu_vec=[0.0, 0.0, 0.0],
    masses=[1.8629520813762488, 1.8629520813762488, 2.7845080443471675],
)

function assert_nan_aware_close(actual::Real, expected::Real; rtol::Float64=1e-6, atol::Float64=1e-8)
    if isnan(expected)
        @test isnan(actual)
    else
        @test isapprox(actual, expected; rtol=rtol, atol=atol)
    end
end

function assert_nan_aware_close(actual::AbstractVector{<:Real}, expected::AbstractVector{<:Real}; rtol::Float64=1e-6, atol::Float64=1e-8)
    @test length(actual) == length(expected)
    for i in eachindex(actual, expected)
        if isnan(expected[i])
            @test isnan(actual[i])
        else
            @test isapprox(actual[i], expected[i]; rtol=rtol, atol=atol)
        end
    end
end

@testset "solver structure-only behavior lock" begin
    @testset "smoke snapshot lock for key fields" begin
        _ = (SNAPSHOT_SOURCE, SNAPSHOT_COMMAND)

        model = Models.create_model(:PNJL)
        mode = Models.FixedMu()
        T_fm = 100.0 / 197.327
        μ_fm = 0.0

        result = Models.solve(model, mode, T_fm, μ_fm; p_num=8, t_num=4, residual_norm_max=1e-6)

        @test result.converged == STRUCTURE_LOCK_BASELINE.converged
        @test length(result.solution) == STRUCTURE_LOCK_BASELINE.solution_length

        assert_nan_aware_close(result.residual_norm, STRUCTURE_LOCK_BASELINE.residual_norm)
        assert_nan_aware_close(result.omega, STRUCTURE_LOCK_BASELINE.omega)
        assert_nan_aware_close(result.pressure, STRUCTURE_LOCK_BASELINE.pressure)
        assert_nan_aware_close(result.rho_norm, STRUCTURE_LOCK_BASELINE.rho_norm)
        assert_nan_aware_close(collect(result.mu_vec), STRUCTURE_LOCK_BASELINE.mu_vec)
        assert_nan_aware_close(collect(result.masses), STRUCTURE_LOCK_BASELINE.masses)
    end
end
