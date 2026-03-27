using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CHECKER_PATH = joinpath(REPO_ROOT, "scripts", "dev", "check_relaxtime_script_governance.jl")
const MANIFEST_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "workflow", "classification_manifest.toml")

@testset "relaxtime script governance checker" begin
    @test isfile(CHECKER_PATH)
    @test isfile(MANIFEST_PATH)

    code = """
        include(raw\"$CHECKER_PATH\")
        ok = Main.RelaxTimeScriptGovernance.check_governance()
        ok || error(\"governance failed\")
    """
    run(`julia --project=. -e $code`)
    @test true
end
