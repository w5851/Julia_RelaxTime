#!/usr/bin/env julia

module AgentInstructionGovernance

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MAX_ROOT_LINES = 220
const REQUIRED_HEADINGS = [
    "## Working Context",
    "## Codex Collaboration Rules",
    "## Command Reference",
    "## Repository Layout",
    "## Commit Messages",
]
const FORBIDDEN_CATALOG_HEADINGS = [
    "## Setup Commands",
    "## Test Commands",
    "## Governance And Audit Commands",
    "## Benchmark Commands",
]
const COMMAND_REFERENCE = "docs/dev/agent_command_reference.md"
const REQUIRED_COMMAND_SECTIONS = [
    "## Environment and setup",
    "## Focused test profiles",
    "## Governance and audit commands",
    "## Benchmark commands",
]

normalize_text(text::AbstractString) = replace(text, "\r\n" => "\n", "\r" => "\n")

function validate_repository(root::AbstractString; max_root_lines::Int=MAX_ROOT_LINES, advisories::Vector{String}=String[])
    violations = String[]
    agents_path = joinpath(root, "AGENTS.md")
    command_path = joinpath(root, split(COMMAND_REFERENCE, '/')...)

    isfile(agents_path) || return ["missing root AGENTS.md"]
    content = normalize_text(read(agents_path, String))
    line_count = length(split(content, '\n'; keepempty=true))
    line_count <= max_root_lines || push!(advisories,
        "root AGENTS.md has $(line_count) lines; consider shortening beyond $(max_root_lines)")

    for heading in REQUIRED_HEADINGS
        occursin(heading, content) || push!(advisories, "root AGENTS.md does not use suggested heading $(heading)")
    end
    for heading in FORBIDDEN_CATALOG_HEADINGS
        occursin(heading, content) && push!(advisories,
            "consider moving $(heading) to the command reference")
    end
    occursin(COMMAND_REFERENCE, content) || push!(violations,
        "root AGENTS.md must link $(COMMAND_REFERENCE)")

    if !isfile(command_path)
        push!(violations, "missing command reference $(COMMAND_REFERENCE)")
    else
        command_content = normalize_text(read(command_path, String))
        for heading in REQUIRED_COMMAND_SECTIONS
            occursin(heading, command_content) || push!(advisories,
                "command reference does not use suggested heading $(heading)")
        end
    end
    return violations
end

function main()
    advisories = String[]
    violations = validate_repository(REPO_ROOT; advisories)
    foreach(item -> println("[agent-instructions] advisory: " * item), advisories)
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
