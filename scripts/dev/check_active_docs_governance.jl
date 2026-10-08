#!/usr/bin/env julia

module ActiveDocsGovernance

using Dates

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const NAME_RE = r"^\d{4}-\d{2}-\d{2}_.+\.md$"
const REVIEW_AGE_DAYS = 60

function review_documents(root::AbstractString=ROOT; current_date::Date=today())
    active_dir = joinpath(root, "docs", "dev", "active")
    violations = String[]
    advisories = String[]
    isdir(active_dir) || return (; violations=["active docs directory not found: $(active_dir)"], advisories)
    for file in sort(readdir(active_dir))
        isfile(joinpath(active_dir, file)) && endswith(file, ".md") || continue
        if !occursin(NAME_RE, file)
            push!(violations, "invalid name format: $(file) (expected YYYY-MM-DD_*.md)")
            continue
        end
        created = try
            Date(first(split(file, '_')), dateformat"yyyy-mm-dd")
        catch
            push!(violations, "invalid date in active document filename: $(file)")
            continue
        end
        if Dates.value(current_date - created) > REVIEW_AGE_DAYS
            push!(advisories, "review active doc (>$(REVIEW_AGE_DAYS)d): $(file); archive by task status, not age")
        end
    end
    return (; violations, advisories)
end

function main()
    result = review_documents()
    foreach(item -> println("[active-docs-governance] advisory: " * item), result.advisories)
    if !isempty(result.violations)
        println("[active-docs-governance] FAILED")
        foreach(item -> println(" - " * item), result.violations)
        return 1
    end
    println("[active-docs-governance] OK")
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(ActiveDocsGovernance.main())
end
