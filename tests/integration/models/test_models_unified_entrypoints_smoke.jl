using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Models unified entrypoints smoke" begin
    @testset "scan entrypoints" begin
        tmp = mktempdir()
        out_tmu = joinpath(tmp, "models_entry_tmu.csv")
        stats_tmu = Models.run_tmu_scan(
            T_values=[150.0],
            mu_values=[0.0],
            xi_values=[0.0],
            output_path=out_tmu,
            overwrite=true,
            resume=false,
            use_phase_aware=false,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=60,
        )
        @test stats_tmu.total == 1
        @test stats_tmu.success == 1

        out_trho = joinpath(tmp, "models_entry_trho.csv")
        stats_trho = Models.run_trho_scan(
            T_values=[150.0],
            rho_values=[0.2],
            xi_values=[0.0],
            output_path=out_trho,
            overwrite=true,
            resume=false,
            reverse_rho=false,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=60,
        )
        @test stats_trho.total == 1
        @test stats_trho.success == 1
    end

    @testset "workflow entrypoint" begin
        T = 150.0 / 197.327
        mu = 0.0
        res = Models.solve_gap_and_transport(
            T,
            mu;
            xi=0.0,
            solver_backend=:legacy,
            compute_tau=false,
            compute_bulk=false,
            tau=(u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0),
            p_num=12,
            t_num=4,
            solver_kwargs=(iterations=60,),
        )
        @test isfinite(res.transport.eta)
        @test (isnan(res.transport.zeta) || isfinite(res.transport.zeta))
        @test isfinite(res.transport.sigma)
    end
end
