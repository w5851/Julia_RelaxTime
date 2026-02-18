#!/usr/bin/env julia

using Dates

const ROOT = pwd()
const ACTIVE_DIR = joinpath(ROOT, "docs", "dev", "active")

const NAME_RE = r"^\d{4}-\d{2}-\d{2}_.+\.md$"
const MAX_ACTIVE_AGE_DAYS = 60

function is_stale(date_str)
    created = Date(String(date_str), dateformat"yyyy-mm-dd")
    return Dates.value(today() - created) > MAX_ACTIVE_AGE_DAYS
end

function main()
    isdir(ACTIVE_DIR) || error("Active docs dir not found: $(ACTIVE_DIR)")

    violations = String[]
    files = sort(readdir(ACTIVE_DIR))
    for file in files
        path = joinpath(ACTIVE_DIR, file)
        isfile(path) || continue
        endswith(file, ".md") || continue

        if !occursin(NAME_RE, file)
            push!(violations, "invalid name format: $(file) (expected YYYY-MM-DD_*.md)")
            continue
        end

        date_str = first(split(file, '_'))
        if is_stale(date_str)
            push!(violations, "stale active doc (>$(MAX_ACTIVE_AGE_DAYS)d): $(file)")
        end
    end

    if !isempty(violations)
        println("[active-docs-governance] FAILED")
        for v in violations
            println(" - " * v)
        end
        println("hint: completed tasks should be archived via scripts/dev/archive_docs.jl")
        exit(1)
    end

    println("[active-docs-governance] OK")
    println("  rule: filename=YYYY-MM-DD_*.md, max_age_days=$(MAX_ACTIVE_AGE_DAYS)")
end

main()
