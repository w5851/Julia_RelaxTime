#!/usr/bin/env julia

"""
Models 统一入口契约检查（轻量骨架）

目的：
1) 确保 `src/models/entrypoints.jl` 公开入口已在 `src/models/Models.jl` 导出。
2) 确保核心 workflow/scan/phase 入口未丢失。
3) 给出文档与调用边界的软告警（不阻断）。

用法：
  julia --project=. scripts/dev/check_models_entry_contract.jl
"""

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS_FILE = joinpath(ROOT, "src", "models", "Models.jl")
const ENTRYPOINTS_FILE = joinpath(ROOT, "src", "models", "entrypoints.jl")
const ENTRY_DOC_FILE = joinpath(ROOT, "docs", "api", "generated", "models", "EntryPointsExportIndex.md")

const REQUIRED_PUBLIC_ENTRYPOINTS = (
    "run_tmu_scan",
    "run_trho_scan",
    "build_default_rho_grid",
    "solve_gap_and_transport",
    "solve_transport_from_equilibrium",
    "solve_gap_and_meson_point",
    "solve_gas_liquid_point",
    "solve_rotation_point",
    "run_phase_pipeline",
    "run_production_phase_pipeline",
    "find_cep",
    "build_phase_artifacts",
    "resolve_phase_output_target",
    "promote_phase_artifacts",
)

function parse_export_symbols(line::AbstractString)
    s = strip(line)
    startswith(s, "export ") || return String[]
    body = strip(s[8:end])
    isempty(body) && return String[]

    out = String[]
    for part in split(body, ',')
        token = strip(part)
        isempty(token) && continue
        name = split(token)[1]
        isempty(name) && continue
        push!(out, name)
    end
    return out
end

function collect_exports(path::String)
    isfile(path) || error("missing file: $(relpath(path, ROOT))")
    syms = Set{String}()
    for line in eachline(path)
        for sym in parse_export_symbols(line)
            push!(syms, sym)
        end
    end
    return syms
end

function collect_julia_files_under(paths::Vector{String})
    files = String[]
    for base in paths
        isdir(base) || continue
        for (root, _, names) in walkdir(base)
            for name in names
                endswith(name, ".jl") || continue
                push!(files, joinpath(root, name))
            end
        end
    end
    return sort(files)
end

function check_boundary_warnings()
    target_dirs = [
        joinpath(ROOT, "src"),
        joinpath(ROOT, "scripts"),
        joinpath(ROOT, "tests"),
    ]
    files = collect_julia_files_under(target_dirs)
    warnings = String[]

    direct_workflow_include_pattern = r"\binclude\(.*(src[\\/]+models[\\/]+workflows[\\/]+(TransportWorkflow|MesonMassWorkflow|GasLiquidWorkflow|RotationWorkflow)\.jl|models[\\/]+workflows[\\/]+(TransportWorkflow|MesonMassWorkflow|GasLiquidWorkflow|RotationWorkflow)\.jl)"
    for file in files
        rel = replace(relpath(file, ROOT), '\\' => '/')
        startswith(rel, "src/models/") && continue
        line_no = 0
        for line in eachline(file)
            line_no += 1
            occursin(direct_workflow_include_pattern, line) || continue
            push!(warnings, "$(rel):$(line_no): direct include of internal workflow file detected")
        end
    end

    return warnings
end

function main()
    violations = String[]
    soft_warnings = String[]

    models_exports = collect_exports(MODELS_FILE)
    entry_exports = collect_exports(ENTRYPOINTS_FILE)

    # Soft check: entrypoints 导出与 Models 导出差异（历史兼容期仅告警，不阻断）
    for sym in sort(collect(entry_exports))
        sym in models_exports || push!(soft_warnings, "entrypoints export not re-exported by Models: $(sym)")
    end

    # Hard check 2: 核心公开入口不得缺失
    for sym in REQUIRED_PUBLIC_ENTRYPOINTS
        sym in entry_exports || push!(violations, "missing required entrypoint export in entrypoints.jl: $(sym)")
        sym in models_exports || push!(violations, "missing required export in Models.jl: $(sym)")
    end

    # Soft check: 入口索引文档是否覆盖关键入口
    if isfile(ENTRY_DOC_FILE)
        doc_content = read(ENTRY_DOC_FILE, String)
        for sym in REQUIRED_PUBLIC_ENTRYPOINTS
            occursin("`$(sym)`", doc_content) || push!(soft_warnings, "EntryPointsExportIndex missing symbol: $(sym)")
        end
    else
        push!(soft_warnings, "missing docs file: $(replace(relpath(ENTRY_DOC_FILE, ROOT), '\\' => '/'))")
    end

    append!(soft_warnings, check_boundary_warnings())

    if !isempty(violations)
        println("[models-entry-contract] FAILED: $(length(violations)) violation(s)")
        for item in violations
            println(" - " * item)
        end
        if !isempty(soft_warnings)
            println("[models-entry-contract] WARNINGS: $(length(soft_warnings))")
            for item in Iterators.take(soft_warnings, 10)
                println(" - " * item)
            end
            length(soft_warnings) > 10 && println(" - ... (truncated)")
        end
        exit(1)
    end

    println("[models-entry-contract] OK")
    println("  entrypoints_exports=$(length(entry_exports)) models_exports=$(length(models_exports))")

    if !isempty(soft_warnings)
        println("[models-entry-contract] WARNINGS: $(length(soft_warnings))")
        for item in Iterators.take(soft_warnings, 10)
            println(" - " * item)
        end
        length(soft_warnings) > 10 && println(" - ... (truncated)")
    end
end

main()
