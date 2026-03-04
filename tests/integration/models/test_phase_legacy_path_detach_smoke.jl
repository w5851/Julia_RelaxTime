using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

function _read_text(path::String)
    return read(path, String)
end

function _all_julia_files(root::String)
    files = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            endswith(lowercase(name), ".jl") || continue
            push!(files, joinpath(dir, name))
        end
    end
    return files
end

@testset "Phase legacy path detach smoke" begin
    cli_file = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")
    pipeline_file = joinpath(PROJECT_ROOT, "src", "models", "phase", "PhasePipeline.jl")
    entrypoints_file = joinpath(PROJECT_ROOT, "src", "models", "entrypoints.jl")

    cli_text = _read_text(cli_file)
    pipeline_text = _read_text(pipeline_file)
    entrypoints_text = _read_text(entrypoints_file)

    @test occursin("Models.run_phase_pipeline", cli_text)
    @test !occursin("PhaseTransition", cli_text)
    @test !occursin("PhaseTransition", pipeline_text)
    @test occursin("export run_phase_pipeline", entrypoints_text)

    script_files = _all_julia_files(joinpath(PROJECT_ROOT, "scripts"))
    legacy_script_refs = String[]
    for file in script_files
        text = _read_text(file)
        occursin("PhaseTransition", text) && push!(legacy_script_refs, file)
    end
    @test isempty(legacy_script_refs)
end
