using Test

module ScriptEntrypointsTest
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_script_entrypoints.jl"))
end

@testset "Script entrypoint syntax and reference checks" begin
    mktempdir() do dir
        script = joinpath(dir, "entry.jl")
        violations = NamedTuple[]
        write(script, "error(\"syntax checks must not execute the script\")\n")
        ScriptEntrypointsTest.scan_file(script, violations)
        @test isempty(violations)

        write(script, "abspath(PROGRAM_FILE) == @__FILE__ && main(ARGS)\n")
        ScriptEntrypointsTest.scan_file(script, violations)
        @test length(violations) == 1
        @test only(violations).kind == "syntax"

        empty!(violations)
        write(script, "function incomplete(\n")
        ScriptEntrypointsTest.scan_file(script, violations)
        @test length(violations) == 1
        @test only(violations).kind == "syntax"

        empty!(violations)
        write(script, "include(\"missing.jl\")\n")
        ScriptEntrypointsTest.scan_file(script, violations)
        @test length(violations) == 1
        @test only(violations).kind == "include"
        write(joinpath(dir, "missing.jl"), "nothing\n")
        empty!(violations)
        ScriptEntrypointsTest.scan_file(script, violations)
        @test isempty(violations)
    end
end
