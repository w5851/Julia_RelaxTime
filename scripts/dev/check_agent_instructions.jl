#!/usr/bin/env julia

module AgentInstructionGovernance

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const COMMAND_REFERENCE = "docs/dev/agent_command_reference.md"

normalize_text(text::AbstractString) = replace(text, "\r\n" => "\n", "\r" => "\n")

function validate_repository(root::AbstractString)
    violations = String[]
    agents_path = joinpath(root, "AGENTS.md")
    command_path = joinpath(root, split(COMMAND_REFERENCE, '/')...)

    isfile(agents_path) || return ["missing root AGENTS.md"]
    content = normalize_text(read(agents_path, String))
    occursin(COMMAND_REFERENCE, content) || push!(violations,
        "root AGENTS.md must link $(COMMAND_REFERENCE)")

    if !isfile(command_path)
        push!(violations, "missing command reference $(COMMAND_REFERENCE)")
    end
    return violations
end

function main()
    violations = validate_repository(REPO_ROOT)
    if !isempty(violations)
        println("[agent-instructions] FAILED: $(length(violations)) violation(s)")
        foreach(item -> println(" - " * item), violations)
        return 1
    end
    println("[agent-instructions] OK")
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(AgentInstructionGovernance.main())
end
