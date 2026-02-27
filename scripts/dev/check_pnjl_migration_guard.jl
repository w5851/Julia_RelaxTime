#!/usr/bin/env julia

using TOML

"""
PNJL 迁移门禁检查（阶段 A）

目的：在“全量迁移到 src/models”期间，限制 `src/pnjl` 下新增核心实现，防止技术债继续增长。

默认检查范围：工作区相对 HEAD 的变更（含未暂存 + 已暂存）。

用法：
  julia --project=. scripts/dev/check_pnjl_migration_guard.jl
  julia --project=. scripts/dev/check_pnjl_migration_guard.jl --staged
  julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base origin/main --head HEAD
"""

const PROTECTED_PREFIXES = (
    "src/pnjl/core/",
    "src/pnjl/solver/",
    "src/pnjl/derivatives/",
)

const ALLOWED_FACADE_SUFFIXES = (
    "Facade.jl",
)

const FORBIDDEN_MODEL_KIND_HARDCODE = r"model_kind\s*(::\s*Symbol)?\s*=\s*:PNJL"
const MODELS_INVOKELATEST_ALLOWLIST_FILE = joinpath("config", "ci", "models_invokelatest_allowlist.toml")
const FORBIDDEN_PNJL_SCAN_RUNTIME_DEP = r"\bPNJL\.run_(tmu|trho)_scan\b"
const FORBIDDEN_PNJL_SCAN_SUBMODULE_DEP = r"\bPNJL\.(TmuScan|TrhoScan)\b|using\s+\.PNJL\.(TmuScan|TrhoScan)\b|using\s+PNJL\.(TmuScan|TrhoScan)\b"
const FORBIDDEN_PNJL_ENTRY_INCLUDE_DEP = r"src/pnjl/PNJL\.jl|joinpath\(.*\"pnjl\"\s*,\s*\"PNJL\.jl\"\)|include\(.*PNJL\.jl\)"
const FORBIDDEN_PNJL_SCAN_DEFAULT_INCLUDE = r"^include\(joinpath\(@__DIR__,\s*\"scans\",\s*\".+\.jl\"\)\)$"
const PNJL_SCAN_DEFAULT_INCLUDE_ALLOWLIST = (
    "ScanConfig.jl",
)

struct GuardConfig
    mode::Symbol  # :working | :staged | :range
    base::Union{Nothing, String}
    head::Union{Nothing, String}
end

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function load_models_invokelatest_allowlist(root::String)
    allowlist_path = joinpath(root, MODELS_INVOKELATEST_ALLOWLIST_FILE)
    isfile(allowlist_path) || error("missing allowlist file: $(MODELS_INVOKELATEST_ALLOWLIST_FILE)")

    raw = TOML.parsefile(allowlist_path)
    table = get(raw, "models_invokelatest_allowlist", nothing)
    table isa AbstractDict || error("invalid allowlist format: [models_invokelatest_allowlist] table is required")

    allowlist = Dict{String, Vector{Regex}}()
    for (path, patterns_any) in table
        patterns_any isa AbstractVector || error("invalid allowlist entry for $(path): value must be an array of regex strings")
        patterns = Regex[]
        for p in patterns_any
            p isa AbstractString || error("invalid allowlist pattern for $(path): expected string")
            push!(patterns, Regex(String(p)))
        end
        allowlist[normalize_path(String(path))] = patterns
    end
    return allowlist
end

normalize_path(path::AbstractString) = replace(path, '\\' => '/')

function parse_args(args::Vector{String})
    staged = false
    base = nothing
    head = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--staged"
            staged = true
            i += 1
        elseif arg == "--base"
            i < length(args) || error("--base requires a value")
            base = args[i + 1]
            i += 2
        elseif arg == "--head"
            i < length(args) || error("--head requires a value")
            head = args[i + 1]
            i += 2
        else
            error("unknown arg: $(arg)")
        end
    end

    if staged && (base !== nothing || head !== nothing)
        error("--staged cannot be used with --base/--head")
    end

    if (base === nothing) ⊻ (head === nothing)
        error("--base and --head must be provided together")
    end

    mode = staged ? :staged : ((base !== nothing && head !== nothing) ? :range : :working)
    return GuardConfig(mode, base, head)
end

function git_read_lines(root::String, args::Vector{String})
    cmd = Cmd(["git", "-C", root, args...])
    out = read(cmd, String)
    return isempty(out) ? String[] : split(chomp(out), '\n')
end

function target_name_status(root::String, cfg::GuardConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--name-status", "--", "src/pnjl"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--name-status", range, "--", "src/pnjl"])
    else
        return git_read_lines(root, ["diff", "--name-status", "HEAD", "--", "src/pnjl"])
    end
end

function target_patch(root::String, cfg::GuardConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", "src/pnjl"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", "src/pnjl"])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", "src/pnjl"])
    end
end

function target_patch_models(root::String, cfg::GuardConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", "src/models"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", "src/models"])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", "src/models"])
    end
end

function target_patch_src(root::String, cfg::GuardConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", "src"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", "src"])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", "src"])
    end
end

function target_patch_workspace(root::String, cfg::GuardConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", "src", "tests", "scripts"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", "src", "tests", "scripts"])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", "src", "tests", "scripts"])
    end
end

function target_patch_allowlist(root::String, cfg::GuardConfig)
    allowlist_rel = MODELS_INVOKELATEST_ALLOWLIST_FILE
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", allowlist_rel])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", allowlist_rel])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", allowlist_rel])
    end
end

function allowlist_change_summary(root::String, cfg::GuardConfig)
    added = String[]
    removed = String[]

    for line in target_patch_allowlist(root, cfg)
        startswith(line, "+++") && continue
        startswith(line, "---") && continue
        startswith(line, "@@") && continue

        if startswith(line, "+")
            s = strip(line[2:end])
            isempty(s) && continue
            push!(added, s)
        elseif startswith(line, "-")
            s = strip(line[2:end])
            isempty(s) && continue
            push!(removed, s)
        end
    end

    return (added=added, removed=removed)
end

function models_invokelatest_audit(root::String, allowlist::Dict{String, Vector{Regex}})
    models_root = joinpath(root, "src", "models")
    allowlist_baseline = sum(length(v) for v in values(allowlist))

    observed = 0
    allowlisted = 0
    non_allowlisted_hits = String[]

    if !isdir(models_root)
        return (observed=observed, allowlist_baseline=allowlist_baseline, allowlisted=allowlisted, non_allowlisted_hits=non_allowlisted_hits)
    end

    for (dir, _, files) in walkdir(models_root)
        for file in files
            endswith(file, ".jl") || continue
            abs_path = joinpath(dir, file)
            rel_path = normalize_path(relpath(abs_path, root))
            allowed_patterns = get(allowlist, rel_path, Regex[])

            line_no = 0
            for line in eachline(abs_path)
                line_no += 1
                if occursin(r"Base\.invokelatest\(", line)
                    observed += 1
                    if any(occursin(re, line) for re in allowed_patterns)
                        allowlisted += 1
                    else
                        push!(non_allowlisted_hits, "$(rel_path):$(line_no): $(strip(line))")
                    end
                end
            end
        end
    end

    return (
        observed=observed,
        allowlist_baseline=allowlist_baseline,
        allowlisted=allowlisted,
        non_allowlisted_hits=non_allowlisted_hits,
    )
end

function pnjl_scan_runtime_dependency_audit(root::String)
    src_root = joinpath(root, "src")
    hits = String[]
    isdir(src_root) || return hits

    for (dir, _, files) in walkdir(src_root)
        for file in files
            endswith(file, ".jl") || continue
            abs_path = joinpath(dir, file)
            rel_path = normalize_path(relpath(abs_path, root))
            startswith(rel_path, "src/pnjl/") && continue

            line_no = 0
            for line in eachline(abs_path)
                line_no += 1
                occursin(FORBIDDEN_PNJL_SCAN_RUNTIME_DEP, line) || continue
                push!(hits, "$(rel_path):$(line_no): $(strip(line))")
            end
        end
    end
    return hits
end

function pnjl_scan_default_include_audit(root::String)
    pnjl_path = joinpath(root, "src", "pnjl", "PNJL.jl")
    isfile(pnjl_path) || return String[]

    hits = String[]
    line_no = 0
    for line in eachline(pnjl_path)
        line_no += 1
        s = strip(line)
        startswith(line, "include(") || continue
        occursin(FORBIDDEN_PNJL_SCAN_DEFAULT_INCLUDE, s) || continue
        if any(occursin("\"$(allowed)\"", s) for allowed in PNJL_SCAN_DEFAULT_INCLUDE_ALLOWLIST)
            continue
        end
        push!(hits, "src/pnjl/PNJL.jl:$(line_no): $(s)")
    end
    return hits
end

function is_protected_file(path::String)
    return any(startswith(path, p) for p in PROTECTED_PREFIXES)
end

function is_allowed_facade_file(path::String)
    return any(endswith(path, s) for s in ALLOWED_FACADE_SUFFIXES)
end

function parse_export_symbols(line::AbstractString)
    s = strip(line)
    startswith(s, "export ") || return String[]
    body = strip(s[8:end])
    isempty(body) && return String[]

    symbols = String[]
    for part in split(body, ',')
        token = strip(part)
        isempty(token) && continue
        name = split(token)[1]
        isempty(name) && continue
        push!(symbols, name)
    end
    return symbols
end

function collect_models_exports(root::String)
    models_path = joinpath(root, "src", "models", "Models.jl")
    isfile(models_path) || return Set{String}()

    exports = Set{String}()
    for line in eachline(models_path)
        for sym in parse_export_symbols(line)
            push!(exports, sym)
        end
    end
    return exports
end

function main()
    cfg = parse_args(ARGS)
    root = project_root()
    allowlist = load_models_invokelatest_allowlist(root)

    violations = String[]
    models_exports = collect_models_exports(root)
    added_pnjl_exports = String[]

    # 规则1：禁止在 src/pnjl 新增文件（迁移期只允许收敛/替换）。
    for line in target_name_status(root, cfg)
        isempty(line) && continue
        fields = split(line, '\t')
        length(fields) >= 2 || continue
        status = fields[1]
        path = normalize_path(fields[end])
        if startswith(status, "A")
            push!(violations, "forbid adding new file under src/pnjl: $(path)")
        end
    end

    # 规则2：受保护目录仅允许 façade 薄层变更；禁止新增 function 定义。
    current_file = ""
    for line in target_patch(root, cfg)
        if startswith(line, "+++ b/")
            current_file = normalize_path(line[7:end])
            continue
        end

        if isempty(current_file)
            continue
        end
        startswith(current_file, "src/pnjl/") || continue

        if startswith(line, "+") && !startswith(line, "+++")
            added = strip(line[2:end])
            isempty(added) && continue
            startswith(added, "#") && continue

            if current_file == "src/pnjl/PNJL.jl"
                append!(added_pnjl_exports, parse_export_symbols(added))
            end

            if is_protected_file(current_file) && !is_allowed_facade_file(current_file)
                if occursin(r"^function\s+", added)
                    push!(violations, "forbid new function in protected legacy module: $(current_file) => $(added)")
                end
            end

            # 规则4：禁止在 src/pnjl 新增 model_kind=:PNJL 硬编码，避免后端语义漂移回归。
            if occursin(FORBIDDEN_MODEL_KIND_HARDCODE, added)
                push!(violations, "forbid hardcoded model_kind=:PNJL under src/pnjl: $(current_file) => $(added)")
            end

            # 规则7：禁止在 PNJL 兼容入口新增 scans 默认 include（避免默认耦合回流）。
            added_raw = line[2:end]
            if current_file == "src/pnjl/PNJL.jl" && startswith(added_raw, "include(") && occursin(FORBIDDEN_PNJL_SCAN_DEFAULT_INCLUDE, added)
                is_allowlisted = any(occursin("\"$(allowed)\"", added) for allowed in PNJL_SCAN_DEFAULT_INCLUDE_ALLOWLIST)
                if !is_allowlisted
                    push!(violations, "forbid new default include under src/pnjl/PNJL.jl scans/*: $(added)")
                end
            end
        end
    end

    # 规则3：新增 PNJL compat export 必须已在 Models 导出（先 models 后 compat）。
    for sym in unique(added_pnjl_exports)
        if !(sym in models_exports)
            push!(violations, "forbid PNJL-only new export without Models export: $(sym) (add to src/models/Models.jl first)")
        end
    end

    # 规则5：src/models 中新增 invokelatest 必须命中白名单（防止新增分散动态调用点）。
    current_file = ""
    for line in target_patch_models(root, cfg)
        if startswith(line, "+++ b/")
            current_file = normalize_path(line[7:end])
            continue
        end

        if isempty(current_file)
            continue
        end
        startswith(current_file, "src/models/") || continue

        if startswith(line, "+") && !startswith(line, "+++")
            added = strip(line[2:end])
            isempty(added) && continue
            startswith(added, "#") && continue

            if occursin(r"Base\.invokelatest\(", added)
                allowed_patterns = get(allowlist, current_file, Regex[])
                is_allowed = any(occursin(re, added) for re in allowed_patterns)
                if !is_allowed
                    push!(violations, "forbid non-allowlisted invokelatest under src/models: $(current_file) => $(added)")
                end
            end
        end
    end

    # 规则6：禁止新增对 PNJL.run_tmu_scan/run_trho_scan 的运行时依赖（收敛到 Models 统一入口）。
    current_file = ""
    for line in target_patch_src(root, cfg)
        if startswith(line, "+++ b/")
            current_file = normalize_path(line[7:end])
            continue
        end

        if isempty(current_file)
            continue
        end
        startswith(current_file, "src/") || continue

        if startswith(line, "+") && !startswith(line, "+++")
            added = strip(line[2:end])
            isempty(added) && continue
            startswith(added, "#") && continue

            if occursin(FORBIDDEN_PNJL_SCAN_RUNTIME_DEP, added)
                push!(violations, "forbid new runtime dependency on PNJL.run_* scan entrypoints: $(current_file) => $(added)")
            end
        end
    end

    # 规则8：禁止新增对 PNJL.TmuScan/PNJL.TrhoScan 子模块的直接依赖（统一走 PNJL 顶层入口）。
    current_file = ""
    for line in target_patch_workspace(root, cfg)
        if startswith(line, "+++ b/")
            current_file = normalize_path(line[7:end])
            continue
        end

        if isempty(current_file)
            continue
        end
        (startswith(current_file, "src/") || startswith(current_file, "tests/") || startswith(current_file, "scripts/")) || continue

        if startswith(line, "+") && !startswith(line, "+++")
            added = strip(line[2:end])
            isempty(added) && continue
            startswith(added, "#") && continue

            if occursin(FORBIDDEN_PNJL_SCAN_SUBMODULE_DEP, added)
                push!(violations, "forbid direct dependency on PNJL.TmuScan/PNJL.TrhoScan submodules: $(current_file) => $(added)")
            end
        end
    end

    # 规则9：禁止新增对 src/pnjl/PNJL.jl 文件路径的直接 include 依赖（统一走 Models 入口）。
    current_file = ""
    for line in target_patch_workspace(root, cfg)
        if startswith(line, "+++ b/")
            current_file = normalize_path(line[7:end])
            continue
        end

        if isempty(current_file)
            continue
        end
        (startswith(current_file, "src/") || startswith(current_file, "tests/") || startswith(current_file, "scripts/")) || continue

        if startswith(line, "+") && !startswith(line, "+++")
            added = strip(line[2:end])
            isempty(added) && continue
            startswith(added, "#") && continue

            if occursin(FORBIDDEN_PNJL_ENTRY_INCLUDE_DEP, added)
                push!(violations, "forbid direct include dependency on src/pnjl/PNJL.jl: $(current_file) => $(added)")
            end
        end
    end

    scan_default_include_hits = pnjl_scan_default_include_audit(root)
    if !isempty(scan_default_include_hits)
        push!(violations, "forbid default include under src/pnjl/PNJL.jl scans/* beyond allowlist (observed=$(length(scan_default_include_hits)))")
    end

    if !isempty(violations)
        println("[pnjl-migration-guard] FAILED")
        for item in violations
            println(" - " * item)
        end
        if !isempty(scan_default_include_hits)
            println("  pnjl-scan-default-include-audit observed=", length(scan_default_include_hits), " scope=src/pnjl/PNJL.jl scans/* (excluding allowlist)")
            for item in Iterators.take(scan_default_include_hits, 5)
                println("   - " * item)
            end
        end
        println("hint: move new core logic to src/models, keep src/pnjl as compatibility thin layer")
        exit(1)
    end

    println("[pnjl-migration-guard] OK")
    if cfg.mode === :staged
        println("  mode=staged")
    elseif cfg.mode === :range
        println("  mode=range ($(cfg.base)...$(cfg.head))")
    else
        println("  mode=working-tree-vs-HEAD")
    end

    audit = models_invokelatest_audit(root, allowlist)
    drift = audit.observed != audit.allowlist_baseline
    status = drift ? "drift" : "stable"
    println("  models-invokelatest-audit=$(status) observed=$(audit.observed) allowlist_baseline=$(audit.allowlist_baseline) allowlisted=$(audit.allowlisted)")

    if !isempty(audit.non_allowlisted_hits)
        println("  models-invokelatest-audit-warning: found non-allowlisted hits (showing up to 5)")
        for item in Iterators.take(audit.non_allowlisted_hits, 5)
            println("   - " * item)
        end
    end

    if drift
        println("  hint: if this drift is intentional, update $(MODELS_INVOKELATEST_ALLOWLIST_FILE) and migration task evidence")
    end

    allowlist_delta = allowlist_change_summary(root, cfg)
    if !isempty(allowlist_delta.added) || !isempty(allowlist_delta.removed)
        println("  models-invokelatest-allowlist-diff added=$(length(allowlist_delta.added)) removed=$(length(allowlist_delta.removed))")
        for item in Iterators.take(allowlist_delta.added, 5)
            println("   + " * item)
        end
        for item in Iterators.take(allowlist_delta.removed, 5)
            println("   - " * item)
        end
    end

    scan_runtime_hits = pnjl_scan_runtime_dependency_audit(root)
    println("  pnjl-scan-runtime-dependency-audit observed=", length(scan_runtime_hits), " scope=src/** (excluding src/pnjl/**)")
    if !isempty(scan_runtime_hits)
        println("  pnjl-scan-runtime-dependency-warning: found PNJL.run_* runtime callsites (showing up to 5)")
        for item in Iterators.take(scan_runtime_hits, 5)
            println("   - " * item)
        end
    end

    println("  pnjl-scan-default-include-audit observed=", length(scan_default_include_hits), " scope=src/pnjl/PNJL.jl scans/* (excluding allowlist)")
end

main()
