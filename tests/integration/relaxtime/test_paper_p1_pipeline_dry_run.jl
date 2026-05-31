using Test

const REPO_ROOT_P1_PIPELINE = normpath(joinpath(@__DIR__, "..", "..", ".."))
const P1_PIPELINE_SCRIPT = joinpath(REPO_ROOT_P1_PIPELINE, "scripts", "relaxtime", "run_paper_p1_pipeline.jl")

@testset "P1 paper pipeline dry-run contract" begin
    @test isfile(P1_PIPELINE_SCRIPT)

    tmp = mktempdir()
    result_dir = joinpath(tmp, "results")
    slice_plan = joinpath(tmp, "mott_slice_plan.csv")
    write(slice_plan, join([
        "muB_MeV,T_min_MeV,T_max_MeV,T_step_MeV",
        "0,100,300,5",
        "1100,30,300,5",
    ], "\n") * "\n")

    mott_output = read(`$(Base.julia_cmd()) --project=$REPO_ROOT_P1_PIPELINE $P1_PIPELINE_SCRIPT --stage=mott --dry-run --result-dir $result_dir --mott-slice-plan $slice_plan --xi-list -0.3,0,0.3`, String)
    @test occursin("mott slice plan:", mott_output)
    @test occursin("mott_muB0000", mott_output)
    @test occursin("mott_muB1100", mott_output)
    @test occursin("write_mott_combined_manifest", mott_output)

    high_mu_cfg = read(joinpath(result_dir, "configs", "mott_muB1100.toml"), String)
    @test occursin("muB_MeV = 1100.0", high_mu_cfg)
    @test occursin("T_min_MeV = 30.0", high_mu_cfg)
    @test occursin("equilibrium_branch_mode = \"stable\"", high_mu_cfg)

    phase_output = read(`$(Base.julia_cmd()) --project=$REPO_ROOT_P1_PIPELINE $P1_PIPELINE_SCRIPT --stage=phase --dry-run --result-dir $result_dir --phase-output-root $(joinpath(tmp, "phase")) --phase-reference-root $(joinpath(tmp, "reference")) --phase-tag dry --xi-list 0`, String)
    @test occursin("--T-min 30.0", phase_output)
    @test occursin("--p-num 24", phase_output)
    @test occursin("--t-num 8", phase_output)
    @test occursin("--tag dry", phase_output)
end
