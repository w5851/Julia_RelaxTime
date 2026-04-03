#!/usr/bin/env julia

"""
Legacy solver switch leakage guard.

Purpose:
- keep transitional legacy switches scoped to explicit compatibility boundaries,
- block accidental spread of `use_problem_spec` / `allow_legacy_path` /
  `warn_on_legacy_path` in `src/**/*.jl` and `tests/**/*.jl`.

Allowlist (current transitional boundary):
- `src/models/solver/Solver.jl`
- `src/models/scans/TrhoScan.jl`
- `src/models/scans/TmuScan.jl`
- `src/models/scans/DualBranchScan.jl`
- `tests/unit/models/test_problem_spec_contract.jl`
- `tests/unit/models/test_trho_scan.jl`
- `tests/unit/models/test_tmu_scan.jl`
- `tests/unit/models/test_dual_branch_scan.jl`
"""

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const TARGET_DIRS = (
    "src",
    "tests",
)

const SWITCH_PATTERNS = (
    :use_problem_spec => r"\buse_problem_spec\b",
    :allow_legacy_path => r"\ballow_legacy_path\b",
    :warn_on_legacy_path => r"\bwarn_on_legacy_path\b",
)

const ALLOWLIST = Set([
    "src/models/solver/Solver.jl",
    "src/models/scans/TrhoScan.jl",
    "src/models/scans/TmuScan.jl",
    "src/models/scans/DualBranchScan.jl",
    "tests/unit/models/test_problem_spec_contract.jl",
    "tests/unit/models/test_trho_scan.jl",
    "tests/unit/models/test_tmu_scan.jl",
    "tests/unit/models/test_dual_branch_scan.jl",
])

normalize_rel(path::AbstractString) = replace(path, '\\' => '/')

function collect_jl_files(root::String)
    files = String[]
    for rel_dir in TARGET_DIRS
        abs_dir = joinpath(root, rel_dir)
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

function detect_switch_hits(root::String, rel_path::String)
    abs_path = joinpath(root, split(rel_path, '/')...)
    hits = NamedTuple{(:line, :switch, :code), Tuple{Int, Symbol, String}}[]
    line_no = 0
    for raw in eachline(abs_path)
        line_no += 1
        for (name, pat) in SWITCH_PATTERNS
            occursin(pat, raw) || continue
            push!(hits, (line=line_no, switch=name, code=strip(raw)))
        end
    end
    return hits
end

function main()
    files = collect_jl_files(ROOT)
    violations = String[]
    allowlist_hits = Dict{String, Int}()

    for rel_path in files
        hits = detect_switch_hits(ROOT, rel_path)
        isempty(hits) && continue

        if rel_path in ALLOWLIST
            allowlist_hits[rel_path] = length(hits)
            continue
        end

        for hit in hits
            push!(violations, "$(rel_path):$(hit.line): detected $(hit.switch) outside allowlist => $(hit.code)")
        end
    end

    if !isempty(violations)
        println("[legacy-switch-governance] FAILED")
        println("  allowlist: " * join(sort!(collect(ALLOWLIST)), ", "))
        for item in violations
            println(" - " * item)
        end
        println("hint: keep transitional legacy switches only in solver boundary and contract test")
        exit(1)
    end

    println("[legacy-switch-governance] OK")
    println("  allowlist_count=" * string(length(ALLOWLIST)))
    for rel_path in sort!(collect(keys(allowlist_hits)))
        println("  hit_count $(rel_path) => $(allowlist_hits[rel_path])")
    end
end

main()
