using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const SCALAR_RTOL = 1e-6
const SCALAR_ATOL = 1e-8

function _assert_scalar_close(a::Real, b::Real)
    @test isapprox(Float64(a), Float64(b); rtol=SCALAR_RTOL, atol=SCALAR_ATOL, nans=true)
end

function _assert_vector_close(a::AbstractVector, b::AbstractVector)
    @test length(a) == length(b)
    for (x, y) in zip(a, b)
        _assert_scalar_close(x, y)
    end
end

function _run_phase_result(runfn; output_dir::String)
    return runfn(
        :PNJL;
        mode=:research,
        T_grid=[150.0],
        rho_grid=[0.1, 0.2, 0.3],
        xi=0.0,
        output_dir=output_dir,
        profile=:smoke,
        solver_backend=:models,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=10,
        compute_crossover=true,
        crossover_method=:inflection,
        crossover_variable=:phi_u,
        crossover_n_mu=6,
        cep_strategy=:interpolate,
        cep_max_bisect_iter=1,
        cep_max_refine_level=0,
        promote_reference=false,
    )
end

function _read_manifest(path::String)
    return JSON3.read(read(path, String))
end

function _is_utc_iso8601(value::AbstractString)
    return occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$", value)
end

@testset "Phase pipeline old/new consistency regression" begin
    old_output = mktempdir()
    new_output = mktempdir()

    @test isdefined(Models, :_run_phase_pipeline_core)
    old_runner = getproperty(Models, :_run_phase_pipeline_core)

    old_result = _run_phase_result(old_runner; output_dir=old_output)
    new_result = _run_phase_result(Models.run_phase_pipeline; output_dir=new_output)

    @test old_result.model_kind == new_result.model_kind
    @test old_result.model_variant == new_result.model_variant
    _assert_scalar_close(old_result.xi, new_result.xi)

    @test old_result.cep.found == new_result.cep.found
    _assert_scalar_close(old_result.cep.T_cep_MeV, new_result.cep.T_cep_MeV)
    _assert_scalar_close(old_result.cep.mu_cep_MeV, new_result.cep.mu_cep_MeV)

    old_boundary = collect(old_result.first_order_boundary)
    new_boundary = collect(new_result.first_order_boundary)
    @test length(old_boundary) == length(new_boundary)
    old_mu = [row.mu_transition_MeV for row in old_boundary]
    new_mu = [row.mu_transition_MeV for row in new_boundary]
    _assert_vector_close(old_mu, new_mu)

    old_crossover = collect(old_result.crossover_line)
    new_crossover = collect(new_result.crossover_line)
    @test length(old_crossover) == length(new_crossover)
    old_tc = [row.T_crossover_MeV for row in old_crossover]
    new_tc = [row.T_crossover_MeV for row in new_crossover]
    _assert_vector_close(old_tc, new_tc)
end

@testset "Phase CLI manifest v1 compatibility regression" begin
    @test isfile(CLI_SCRIPT)
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --preset=smoke --iterations=8 --p_num=8 --t_num=4 --cep_max_bisect_iter=1 --cep_max_refine_level=0 --output_dir=$(output_dir)`
    run(cmd)

    manifest_path = joinpath(output_dir, "run_manifest.json")
    @test isfile(manifest_path)
    manifest = _read_manifest(manifest_path)

    # legacy keys
    @test haskey(manifest, "generated_at")
    @test haskey(manifest, "argv")
    @test haskey(manifest, "config_hash")
    @test haskey(manifest, "run_id")
    @test haskey(manifest, "mode")
    @test haskey(manifest, "model_kind")
    @test haskey(manifest, "artifact_paths")
    @test haskey(manifest, "effective_config")

    # manifest_v1 required checks
    @test haskey(manifest, "pipeline")
    pipe = manifest["pipeline"]
    @test haskey(pipe, "name")
    @test haskey(pipe, "version")
    @test haskey(pipe, "model_kind")
    @test haskey(pipe, "run_id")
    @test haskey(pipe, "git_commit")
    @test haskey(pipe, "manifest_schema_version")
    @test haskey(pipe, "timestamp")
    @test haskey(pipe, "config_hash")
    @test haskey(pipe, "artifact_hash")

    @test _is_utc_iso8601(String(pipe["timestamp"]))
    @test !isempty(strip(String(pipe["config_hash"])))
    @test !isempty(strip(String(pipe["artifact_hash"])))

    # old/new key-field consistency
    @test String(manifest["model_kind"]) == String(pipe["model_kind"])
    @test !isempty(strip(String(manifest["run_id"])))
    @test !isempty(strip(String(pipe["run_id"])))
    @test !isempty(strip(String(manifest["config_hash"])))
    @test !isempty(strip(String(pipe["config_hash"])))
end
