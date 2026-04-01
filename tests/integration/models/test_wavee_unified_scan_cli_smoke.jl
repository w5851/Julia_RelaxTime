using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "models", "run_unified_scan.jl")

@testset "Wave-E unified scan CLI smoke" begin
    @test isfile(CLI_SCRIPT)

    tmpdir = mktempdir()

    tmu_output = joinpath(tmpdir, "wavee_tmu.csv")
    tmu_cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) scan tmu --model_kind=PNJL --T_values=150 --mu_values=0 --xi_values=0 --output_path=$(tmu_output) --overwrite=true --resume=false --use_phase_aware=false --solver_backend=legacy --p_num=12 --t_num=4 --iterations=40`
    run(tmu_cmd)
    @test isfile(tmu_output)
    @test length(readlines(tmu_output)) >= 2

    trho_output = joinpath(tmpdir, "wavee_trho.csv")
    trho_cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) scan trho --model_kind=PNJL --T_values=150 --rho_values=0.2 --xi_values=0 --output_path=$(trho_output) --overwrite=true --resume=false --reverse_rho=false --seed_policy=hybrid_continuity --constraint_mode=fixed_rho --solver_backend=legacy --p_num=12 --t_num=4 --iterations=40`
    run(trho_cmd)
    @test isfile(trho_output)
    @test length(readlines(trho_output)) >= 2

    phase_output_dir = joinpath(tmpdir, "wavee_phase")
    phase_cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) workflow phase --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --profile=smoke --solver_backend=legacy --iterations=10 --output_dir=$(phase_output_dir)`
    run(phase_cmd)

    @test isfile(joinpath(phase_output_dir, "phase_summary.json"))
    @test isfile(joinpath(phase_output_dir, "run_manifest.json"))

    err = try
        bad_cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) scan tmu --model_kind=pnjl_aniso --T_values=150 --mu_values=0 --xi_values=0.1 --output_path=$(joinpath(tmpdir, "bad_aniso.csv")) --overwrite=true --resume=false --use_phase_aware=false --solver_backend=models --p_num=8 --t_num=4 --iterations=20`
        run(bad_cmd)
        nothing
    catch exc
        exc
    end

    @test err !== nothing
    msg = lowercase(sprint(showerror, err))
    @test occursin("pnjl_aniso", msg)

    err_trho = try
        bad_trho_cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) scan trho --model_kind=pnjl_aniso --T_values=150 --rho_values=0.2 --xi_values=0.1 --output_path=$(joinpath(tmpdir, "bad_aniso_trho.csv")) --overwrite=true --resume=false --reverse_rho=false --seed_policy=hybrid_continuity --constraint_mode=fixed_rho --solver_backend=models --p_num=8 --t_num=4 --iterations=20`
        run(bad_trho_cmd)
        nothing
    catch exc
        exc
    end

    @test err_trho !== nothing
    msg_trho = lowercase(sprint(showerror, err_trho))
    @test occursin("pnjl_aniso", msg_trho)
end
