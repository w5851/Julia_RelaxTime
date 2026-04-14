using Test
using JSON3

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Workflow adapter pipeline smoke" begin
    tmp = mktempdir()

    result = Models.run_workflow_pipeline(
        :transport;
        T_fm=0.15,
        mu_fm=0.0,
        xi=0.0,
        output_dir=tmp,
        compute_tau=false,
        tau=(u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0),
        compute_bulk=false,
    )

    @test hasproperty(result, :run_manifest)
    manifest_path = result.run_manifest
    @test isfile(manifest_path)

    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, :pipeline)
    @test String(manifest.pipeline.pipeline_family) == "workflow"
end

@testset "Workflow adapter diagnostics unavailable is auditable" begin
    tmp = mktempdir()
    diag_dir = joinpath(tmp, "diag")
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    lib_path = joinpath(repo_root, "scripts", "analysis", "relaxtime", "t190_sigma_chain_decomposition_lib.jl")
    backup_path = lib_path * ".bak_unavailable"

    mv(lib_path, backup_path; force=true)
    try
        result = Models.run_workflow_pipeline(
            :transport;
            T_fm=0.15,
            mu_fm=0.0,
            xi=0.0,
            output_dir=tmp,
            compute_tau=false,
            tau=(u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0),
            compute_bulk=false,
            diagnostics_mode=:t190_chain,
            diagnostics_output_dir=diag_dir,
            diagnostics_strict=false,
        )

        @test result.diagnostics.status == :unavailable
        @test isfile(String(result.diagnostics.index_path))

        index_json = JSON3.read(read(String(result.diagnostics.index_path), String))
        @test String(index_json.status) == "unavailable"
        @test haskey(index_json, :error)

        manifest = JSON3.read(read(String(result.run_manifest), String))
        @test String(manifest.pipeline.diagnostics_status) == "unavailable"
    finally
        mv(backup_path, lib_path; force=true)
    end
end

@testset "Workflow adapter diagnostics index smoke" begin
    tmp = mktempdir()
    diag_dir = joinpath(tmp, "diag")

    result = Models.run_workflow_pipeline(
        :transport;
        T_fm=0.15,
        mu_fm=0.0,
        xi=0.0,
        output_dir=tmp,
        compute_tau=false,
        tau=(u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0),
        compute_bulk=false,
        diagnostics_mode=:t190_chain,
        diagnostics_output_dir=diag_dir,
        diagnostics_strict=true,
    )

    @test hasproperty(result, :diagnostics)
    @test hasproperty(result.diagnostics, :mode)
    @test result.diagnostics.mode == :t190_chain
    @test hasproperty(result.diagnostics, :index_path)
    @test isfile(String(result.diagnostics.index_path))
    @test hasproperty(result.diagnostics, :artifacts)
    @test haskey(result.diagnostics.artifacts, "t190_mixed_p_chain_csv")
    @test haskey(result.diagnostics.artifacts, "t190_mixed_p_chain_summary_csv")
    @test isfile(String(result.diagnostics.artifacts["t190_mixed_p_chain_csv"]))
    @test isfile(String(result.diagnostics.artifacts["t190_mixed_p_chain_summary_csv"]))

    index_json = JSON3.read(read(String(result.diagnostics.index_path), String))
    @test String(index_json.mode) == "t190_chain"
    @test haskey(index_json, :status)
    @test haskey(index_json, :run_context)
    @test haskey(index_json, :artifacts)
    @test haskey(index_json.artifacts, :t190_mixed_p_chain_csv)
    @test haskey(index_json.artifacts, :t190_mixed_p_chain_summary_csv)

    manifest_path = result.run_manifest
    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, :pipeline)
    @test haskey(manifest.pipeline, :diagnostics_mode)
    @test haskey(manifest.pipeline, :diagnostics_status)
    @test haskey(manifest.pipeline, :diagnostics_index_path)
    @test haskey(manifest.pipeline, :diagnostics_t190_ratio_detM_area_B_over_A)
    @test haskey(manifest.pipeline, :diagnostics_t190_ratio_Dmixed_area_B_over_A)
    @test haskey(manifest.pipeline, :diagnostics_t190_ratio_detM_point_B_over_A)
    @test haskey(manifest.pipeline, :diagnostics_t190_ratio_Dmixed_point_B_over_A)
    @test String(manifest.pipeline.diagnostics_mode) == "t190_chain"
    @test String(manifest.pipeline.diagnostics_status) == "success"
    @test endswith(String(manifest.pipeline.diagnostics_index_path), "diagnostics_index.json")
    @test Float64(manifest.pipeline.diagnostics_t190_ratio_detM_area_B_over_A) > 1.0
    @test Float64(manifest.pipeline.diagnostics_t190_ratio_Dmixed_area_B_over_A) < 1.0
    @test Float64(manifest.pipeline.diagnostics_t190_ratio_detM_point_B_over_A) > 1.0
    @test Float64(manifest.pipeline.diagnostics_t190_ratio_Dmixed_point_B_over_A) < 1.0
end
