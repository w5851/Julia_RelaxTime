using Test

const AGENT_INSTRUCTION_CHECKER = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_agent_instructions.jl"))
if !isdefined(Main, :AgentInstructionGovernance)
    include(AGENT_INSTRUCTION_CHECKER)
end

const _AGENT_REQUIRED = """
# AGENTS.md
## Working Context
## Codex Collaboration Rules
## Command Reference
See docs/dev/agent_command_reference.md.
## Repository Layout
## Commit Message Governance (Mandatory)
"""

const _COMMAND_REQUIRED = """
# Command reference
## Environment and setup
## Focused test profiles
## Governance and audit commands
## Benchmark commands
"""

function _write_agent_fixture(root; agents=_AGENT_REQUIRED, command_reference=_COMMAND_REQUIRED)
    mkpath(joinpath(root, "docs", "dev"))
    write(joinpath(root, "AGENTS.md"), agents)
    command_reference === nothing || write(joinpath(root, "docs", "dev", "agent_command_reference.md"), command_reference)
    return root
end

@testset "agent instruction governance" begin
    mktempdir() do temp
        valid = _write_agent_fixture(joinpath(temp, "valid"))
        @test isempty(Main.AgentInstructionGovernance.validate_repository(valid; max_root_lines=20))

        embedded = _write_agent_fixture(joinpath(temp, "embedded"); agents=_AGENT_REQUIRED * "\n## Test Commands\n")
        violations = Main.AgentInstructionGovernance.validate_repository(embedded; max_root_lines=20)
        @test any(item -> occursin("embedding ## Test Commands", item), violations)

        too_long = _write_agent_fixture(joinpath(temp, "too_long"); agents=_AGENT_REQUIRED * repeat("extra\n", 20))
        violations = Main.AgentInstructionGovernance.validate_repository(too_long; max_root_lines=10)
        @test any(item -> occursin("maximum is 10", item), violations)

        missing_reference = _write_agent_fixture(joinpath(temp, "missing_reference"); command_reference=nothing)
        violations = Main.AgentInstructionGovernance.validate_repository(missing_reference; max_root_lines=20)
        @test any(item -> occursin("missing command reference", item), violations)
    end
end
