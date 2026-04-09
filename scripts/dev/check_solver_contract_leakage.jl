#!/usr/bin/env julia

"""
Solver contract leakage guard.

Purpose:
- Keep `phase/scans/workflows` consumption restricted to stable solver contracts.
- Block direct reads of solver-private attempt origin fields and selector internals.

Policy:
- Allowed in `src/models/solver/**` and tests (internal implementation/testing scope).
- Disallowed in runtime consumers under `src/models/phase/**`, `src/models/scans/**`, `src/models/workflows/**`.
"""

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const RUNTIME_TARGET_DIRS = (
    "src/models/phase",
    "src/models/scans",
    "src/models/workflows",
)

const PRIVATE_PATTERNS = (
    :governed_attempt_origin => r"\bgoverned_attempt_origin\b",
    :fixedrho_joint_attempt_origin => r"\bfixedrho_joint_attempt_origin\b",
    :entropy_attempt_origin => r"\bentropy_attempt_origin\b",
    :sigma_attempt_origin => r"\bsigma_attempt_origin\b",
    :asym_attempt_origin => r"\basym_attempt_origin\b",
    :selection_reason_source => r"\bselection_reason_source\b",
)

const RUNTIME_ALLOWLIST = Set([
    "src/models/phase/PMPhaseDiagnostic.jl",
])

normalize_rel(path::AbstractString) = replace(path, '\\' => '/')

function collect_runtime_files(root::String)
    files = String[]
    for rel_dir in RUNTIME_TARGET_DIRS
        abs_dir = joinpath(root, split(rel_dir, '/')...)
        isdir(abs_dir) || continue
        for (dir, _, names) in walkdir(abs_dir)
            for name in names
                endswith(name, ".jl") || continue
                abs_path = joinpath(dir, name)
                push!(files, normalize_rel(relpath(abs_path, root)))
            end
        end
    end
    sort!(files)
    return files
end

function detect_private_hits(root::String, rel_path::String)
    abs_path = joinpath(root, split(rel_path, '/')...)
    hits = NamedTuple{(:line, :field, :code), Tuple{Int, Symbol, String}}[]
    line_no = 0
    for raw in eachline(abs_path)
        line_no += 1
        for (name, pat) in PRIVATE_PATTERNS
            occursin(pat, raw) || continue
            push!(hits, (line=line_no, field=name, code=strip(raw)))
        end
    end
    return hits
end

function main()
    files = collect_runtime_files(ROOT)
    violations = String[]
    allowlist_hits = Dict{String, Int}()

    for rel_path in files
        hits = detect_private_hits(ROOT, rel_path)
        isempty(hits) && continue

        if rel_path in RUNTIME_ALLOWLIST
            allowlist_hits[rel_path] = length(hits)
            continue
        end

        for hit in hits
            push!(violations, "$(rel_path):$(hit.line): solver private field $(hit.field) leakage => $(hit.code)")
        end
    end

    if !isempty(violations)
        println("[solver-contract-leakage] FAILED")
        println("  runtime_allowlist: " * join(sort!(collect(RUNTIME_ALLOWLIST)), ", "))
        for item in violations
            println(" - " * item)
        end
        println("hint: runtime consumers must use stable diagnostic/result contracts only")
        exit(1)
    end

    println("[solver-contract-leakage] OK")
    println("  runtime_allowlist_count=" * string(length(RUNTIME_ALLOWLIST)))
    for rel_path in sort!(collect(keys(allowlist_hits)))
        println("  hit_count $(rel_path) => $(allowlist_hits[rel_path])")
    end
end

main()
