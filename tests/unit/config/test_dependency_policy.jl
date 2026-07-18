using Test

const DEPENDENCY_POLICY_CHECKER = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_dependency_policy.jl"))
if !isdefined(Main, :DependencyPolicyGovernance)
    include(DEPENDENCY_POLICY_CHECKER)
end

function _write_dependency_fixture(root; root_quadgk=false, runtime_quadgk=false, benchmark_quadgk=true)
    mkpath(joinpath(root, "config", "ci"))
    mkpath(joinpath(root, "benchmark"))
    mkpath(joinpath(root, "src"))
    write(joinpath(root, "Project.toml"), root_quadgk ? "[deps]\nQuadGK = \"uuid\"\n" : "[deps]\n")
    write(joinpath(root, "benchmark", "Project.toml"), benchmark_quadgk ? "[deps]\nQuadGK = \"uuid\"\n" : "[deps]\n")
    forbidden_fixture = "using " * "QuadGK: quad" * "gk\nquad" * "gk(identity, 0, 1)\n"
    write(joinpath(root, "src", "fixture.jl"), runtime_quadgk ? forbidden_fixture : "identity(1)\n")
    write(joinpath(root, "config", "ci", "dependency_policy.toml"), """
schema_version = 1
root_project = "Project.toml"
runtime_scan_roots = ["src", "scripts", "tests"]
forbidden_root_direct_dependencies = ["QuadGK"]

[isolated_dependencies.QuadGK]
project = "benchmark/Project.toml"
forbidden_symbols = ["quadgk"]
""")
    return root
end

@testset "dependency policy governance" begin
    mktempdir() do temp
        valid = _write_dependency_fixture(joinpath(temp, "valid"))
        @test isempty(Main.DependencyPolicyGovernance.validate_repository(valid))

        root_dependency = _write_dependency_fixture(joinpath(temp, "root_dependency"); root_quadgk=true)
        violations = Main.DependencyPolicyGovernance.validate_repository(root_dependency)
        @test any(item -> occursin("root Project.toml [deps]", item), violations)

        runtime_import = _write_dependency_fixture(joinpath(temp, "runtime_import"); runtime_quadgk=true)
        violations = Main.DependencyPolicyGovernance.validate_repository(runtime_import)
        @test any(item -> occursin("src/fixture.jl:1", item), violations)
        @test any(item -> occursin("src/fixture.jl:2", item), violations)

        missing_isolated = _write_dependency_fixture(joinpath(temp, "missing_isolated"); benchmark_quadgk=false)
        violations = Main.DependencyPolicyGovernance.validate_repository(missing_isolated)
        @test any(item -> occursin("must be declared in isolated project", item), violations)
    end
end
