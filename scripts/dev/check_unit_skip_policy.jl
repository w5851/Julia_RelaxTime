#!/usr/bin/env julia

using TOML

const ROOT = pwd()
const RUNTESTS_PATH = joinpath(ROOT, "tests", "unit", "runtests.jl")
const POLICY_PATH = joinpath(ROOT, "config", "ci", "unit_skip_policy.toml")

function read_default_skip_entries(path::String)
    content = read(path, String)
    m = match(r"const\s+DEFAULT_SKIP\s*=\s*Set\(\[(.*?)\]\)"s, content)
    m === nothing && error("Failed to locate DEFAULT_SKIP in $(path)")

    body = m.captures[1]
    entries = String[]
    for line in split(body, '\n')
        stripped = strip(line)
        startswith(stripped, "#") && continue
        isempty(stripped) && continue
        hit = match(r"\"([^\"]+)\"", stripped)
        hit === nothing && continue
        push!(entries, hit.captures[1])
    end
    return sort(entries)
end

function main()
    policy = TOML.parsefile(POLICY_PATH)
    max_skip = Int(policy["max_skip"])
    required_entries = sort(String.(policy["targets"]["required_entries"]))

    current_entries = read_default_skip_entries(RUNTESTS_PATH)

    violations = String[]
    if length(current_entries) > max_skip
        push!(violations, "DEFAULT_SKIP size $(length(current_entries)) exceeds policy max_skip=$(max_skip)")
    end

    if current_entries != required_entries
        push!(violations, "DEFAULT_SKIP entries differ from config/ci/unit_skip_policy.toml")
        only_in_current = setdiff(current_entries, required_entries)
        only_in_policy = setdiff(required_entries, current_entries)
        if !isempty(only_in_current)
            push!(violations, "  extra in runtests: " * join(only_in_current, ", "))
        end
        if !isempty(only_in_policy)
            push!(violations, "  missing in runtests: " * join(only_in_policy, ", "))
        end
    end

    if !isempty(violations)
        println("[unit-skip-policy] FAILED")
        for item in violations
            println(" - " * item)
        end
        exit(1)
    end

    println("[unit-skip-policy] OK")
    println("  phase = " * string(policy["phase"]["name"]))
    println("  max_skip = " * string(max_skip))
    println("  entries = " * join(current_entries, ", "))
end

main()
