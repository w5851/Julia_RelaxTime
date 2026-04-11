using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const ETA_RTOL = 1e-6
const ETA_ATOL = 1e-8

function _write_minimal_cross_section_config(path::String)
    open(path, "w") do io
        write(io, """
schema_version = "v1"

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0]
xi_list = [0.0]
processes = ["ud_to_ud"]

[scan.cross_section.energy]
mode = "list"
sqrt_s_list_MeV = [500.0]
""")
    end
    return path
end

@testset "workflow/scan/orchestrator consistency regression (task-234)" begin
    @testset "workflow eta consistency: solve_gap_and_transport vs run_workflow_pipeline" begin
        tau0 = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
        T_fm = 0.15
        mu_fm = 0.0

        direct = Models.solve_gap_and_transport(
            T_fm,
            mu_fm;
            xi=0.0,
            compute_tau=false,
            tau=tau0,
            compute_bulk=false,
            p_num=12,
            t_num=4,
        )

        tmp = mktempdir()
        via_pipeline = Models.run_workflow_pipeline(
            :transport;
            T_fm=T_fm,
            mu_fm=mu_fm,
            xi=0.0,
            output_dir=tmp,
            compute_tau=false,
            tau=tau0,
            compute_bulk=false,
            p_num=12,
            t_num=4,
        )

        @test isapprox(
            Float64(direct.transport.eta),
            Float64(via_pipeline.transport.eta);
            rtol=ETA_RTOL,
            atol=ETA_ATOL,
        )
        @test isfile(String(via_pipeline.run_manifest))
    end

    @testset "scan stats consistency: run_tmu_scan vs run_scan_pipeline(:tmu)" begin
        tmp = mktempdir()
        output_direct = joinpath(tmp, "scan_tmu_direct.csv")
        output_pipeline = joinpath(tmp, "scan_tmu_pipeline.csv")

        common_kwargs = (
            model_kind=:PNJL,
            T_values=[150.0],
            mu_values=[0.0],
            xi_values=[0.0],
            overwrite=true,
            resume=false,
            use_phase_aware=false,
            solver_backend=:models,
            p_num=12,
            t_num=4,
            iterations=40,
        )

        direct = Models.run_tmu_scan(; common_kwargs..., output_path=output_direct)
        via_pipeline = Models.run_scan_pipeline(:tmu; common_kwargs..., output_path=output_pipeline)

        @test direct.total == via_pipeline.total
        @test direct.success == via_pipeline.success
        @test direct.failure == via_pipeline.failure
        @test direct.skipped == via_pipeline.skipped
        @test isfile(String(via_pipeline.manifest_path))
    end

    @testset "orchestrator manifest/path consistency: run_relaxtime_orchestrator_pipeline" begin
        outdir = mktempdir()
        cfg_dir = mktempdir()
        cfg = _write_minimal_cross_section_config(joinpath(cfg_dir, "regression_cfg.toml"))

        result = Models.run_relaxtime_orchestrator_pipeline(
            :cross_section;
            config_path=cfg,
            outdir=outdir,
            processes=["ud_to_ud"],
            overwrite=true,
            resume=true,
            fail_on_fallback=true,
        )

        @test isfile(String(result.manifest_path))
        @test isfile(String(result.effective_config_path))
        @test isfile(String(result.consumption_report_path))
        @test isfile(String(result.cross_section_path))

        manifest = JSON3.read(read(String(result.manifest_path), String))
        @test haskey(manifest, :pipeline)
        @test String(manifest.pipeline.pipeline_family) == "relaxtime_orchestrator"
    end
end
