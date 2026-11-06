using Test

const AGENT_INSTRUCTION_CHECKER = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_agent_instructions.jl"))
if !isdefined(Main, :AgentInstructionGovernance)
    include(AGENT_INSTRUCTION_CHECKER)
end

const _AGENT_REQUIRED = "# Project guidance\nSee docs/dev/agent_command_reference.md.\n"
const _COMMAND_REQUIRED = "# Commands\n"

function _write_agent_fixture(root; agents=_AGENT_REQUIRED, command_reference=_COMMAND_REQUIRED)
    mkpath(joinpath(root, "docs", "dev"))
    write(joinpath(root, "AGENTS.md"), agents)
    command_reference === nothing || write(joinpath(root, "docs", "dev", "agent_command_reference.md"), command_reference)
    return root
end

@testset "agent instruction governance" begin
    mktempdir() do temp
        valid = _write_agent_fixture(joinpath(temp, "valid"))
        @test isempty(Main.AgentInstructionGovernance.validate_repository(valid))

        unlinked = _write_agent_fixture(joinpath(temp, "unlinked"); agents="# Project guidance\n")
        @test any(item -> occursin("must link", item), Main.AgentInstructionGovernance.validate_repository(unlinked))

        missing_reference = _write_agent_fixture(joinpath(temp, "missing_reference"); command_reference=nothing)
        @test any(item -> occursin("missing command reference", item), Main.AgentInstructionGovernance.validate_repository(missing_reference))

        rm(joinpath(valid, "AGENTS.md"))
        @test Main.AgentInstructionGovernance.validate_repository(valid) == ["missing root AGENTS.md"]
    end
end
