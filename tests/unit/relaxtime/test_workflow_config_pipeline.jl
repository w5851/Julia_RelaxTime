using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const WORKFLOW_CONFIG_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "config", "WorkflowConfig.jl")
@test isfile(WORKFLOW_CONFIG_PATH)
include(WORKFLOW_CONFIG_PATH)

using .WorkflowConfig: normalize_merge_validate

@testset "workflow config pipeline order" begin
    default_cfg = Dict(
        "schema_version" => "v1",
        "scan" => Dict("transport" => Dict("muB_MeV" => 0.0, "xi_list" => [0.0]))
    )

    toml_cfg = Dict(
        "scan" => Dict("transport" => Dict("xi_list" => [-0.3, 0.0, 0.3]))
    )

    cli_cfg = Dict(
        "scan" => Dict("transport" => Dict("muB_MeV" => 100.0))
    )

    aliases = Dict("process_aliases" => Dict("mappings" => Dict("udtoud" => "ud_to_ud")))

    out = normalize_merge_validate(default_cfg, toml_cfg, cli_cfg, aliases)

    @test out.trace == [
        :alias_normalization,
        :source_validation,
        :precedence_merge,
        :effective_validation,
    ]
    @test out.effective["scan"]["transport"]["muB_MeV"] == 100.0
    @test out.effective["scan"]["transport"]["xi_list"] == [-0.3, 0.0, 0.3]
end
