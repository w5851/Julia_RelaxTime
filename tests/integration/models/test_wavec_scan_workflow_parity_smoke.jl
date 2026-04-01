using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Wave-C scan/workflow parity smoke" begin
    @test isdefined(Main.Models, :scan_workflow_migration_status)

    migration = Main.Models.scan_workflow_migration_status("scripts/pnjl/run_tmu_scan.jl")
    @test migration.status in (:hard_deprecated, :removed, :archived)
    @test migration.route in (:compat_adapter, :unified_cli)
    @test migration.removal_wave in (:D, :E)

    tmu_script_path = joinpath(PROJECT_ROOT, "scripts", "pnjl", "run_tmu_scan.jl")
    if !isfile(tmu_script_path)
        @test migration.status in (:removed, :archived)
        @test occursin("scripts/models/run_unified_scan.jl", migration.unified_entry)
        return
    end

    tmu_script = Module(:WaveCTmuScript)
    Base.include(tmu_script, tmu_script_path)

    @test isdefined(tmu_script, :ScriptTmuScanConfig)
    @test hasfield(getfield(tmu_script, :ScriptTmuScanConfig), :model_kind)

    cfg = getfield(tmu_script, :parse_args)([
        "--mode=point",
        "--T_mev=150",
        "--mu_mev=0",
        "--xi=0.0",
        "--no_phase_aware",
        "--solver_backend=models",
        "--model_kind=PNJL",
        "--p_num=8",
        "--t_num=4",
        "--verbose",
    ])
    @test getfield(cfg, :model_kind) == :PNJL

    mktempdir() do outdir
        script_out = joinpath(outdir, "script_tmu.csv")
        unified_out = joinpath(outdir, "unified_tmu.csv")

        via_script = getfield(tmu_script, :main)([
            "--mode=point",
            "--T_mev=150",
            "--mu_mev=0",
            "--xi=0.0",
            "--output=$(script_out)",
            "--overwrite",
            "--no_phase_aware",
            "--solver_backend=models",
            "--model_kind=PNJL",
            "--p_num=8",
            "--t_num=4",
        ])

        via_unified = Main.Models.run_tmu_scan(
            T_values=[150.0],
            mu_values=[0.0],
            xi_values=[0.0],
            output_path=unified_out,
            overwrite=true,
            resume=false,
            use_phase_aware=false,
            solver_backend=:models,
            model_kind=:PNJL,
            p_num=8,
            t_num=4,
        )

        @test via_script.total == via_unified.total == 1
        @test via_script.success == via_unified.success
        @test via_script.failure == via_unified.failure
    end

    dense_script_path = joinpath(PROJECT_ROOT, "scripts", "pnjl", "run_dense_trho_scan.jl")
    dense_script = Module(:WaveCDenseScript)
    Base.include(dense_script, dense_script_path)
    @test isdefined(dense_script, :DenseScanOptions)
    @test hasfield(getfield(dense_script, :DenseScanOptions), :model_kind)

    dense_cfg = getfield(dense_script, :parse_args)([
        "--tmin", "130.0",
        "--tmax", "130.0",
        "--tstep", "1.0",
        "--rho-max", "0.2",
        "--coarse-step", "0.2",
        "--model-kind", "rpnjl",
    ])
    @test getfield(dense_cfg, :model_kind) == :RPNJL

    dense_err = try
        getfield(dense_script, :parse_args)(["--model-kind", "bad_model"])
        nothing
    catch exc
        exc
    end
    @test dense_err isa ArgumentError
    @test occursin("accepted: PNJL|RPNJL", sprint(showerror, dense_err))

    adaptive_script_path = joinpath(PROJECT_ROOT, "scripts", "pnjl", "run_adaptive_trho_scan.jl")
    adaptive_script = Module(:WaveCAdaptiveScript)
    Base.include(adaptive_script, adaptive_script_path)
    @test isdefined(adaptive_script, :AdaptiveCLIOptions)
    @test hasfield(getfield(adaptive_script, :AdaptiveCLIOptions), :model_kind)

    adaptive_cfg = getfield(adaptive_script, :parse_args)([
        "--source", joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl", "scan", "trho", "trho_scan.csv"),
        "--model-kind", "pnjl",
    ])
    @test getfield(adaptive_cfg, :model_kind) == :PNJL

    adaptive_err = try
        getfield(adaptive_script, :parse_args)([
            "--source", joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl", "scan", "trho", "trho_scan.csv"),
            "--model-kind", "bad_model",
        ])
        nothing
    catch exc
        exc
    end
    @test adaptive_err isa ArgumentError
    @test occursin("accepted: PNJL|RPNJL", sprint(showerror, adaptive_err))
end
