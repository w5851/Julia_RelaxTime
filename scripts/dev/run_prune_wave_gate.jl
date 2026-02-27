#!/usr/bin/env julia

"""
PNJL 裁剪波次门禁（删除前快照 + 删除后 smoke + 失败回滚）

用法：
  julia --project=. scripts/dev/run_prune_wave_gate.jl
  julia --project=. scripts/dev/run_prune_wave_gate.jl --staged
  julia --project=. scripts/dev/run_prune_wave_gate.jl --base origin/main --head HEAD
  julia --project=. scripts/dev/run_prune_wave_gate.jl --auto-rollback
  julia --project=. scripts/dev/run_prune_wave_gate.jl --unit-profile smoke

说明：
- 默认模式：working tree vs HEAD（用于本地裁剪波次验收）。
- 快照产物：data/outputs/results/pnjl_prune_wave_snapshot_<timestamp>.txt
- 路径门禁：新增脚本默认输出不得回流到根级 outputs/results（由 check_data_output_path_guard.jl 校验）
- smoke 命令：`UNIT_PROFILE=<profile> julia --project=. tests/unit/runtests.jl`
- 可选 `--auto-rollback` 仅在 working 模式生效，失败时执行：
  `git restore --worktree --staged --source=HEAD -- src/pnjl`
"""

using Dates

const MODELS_INVOKELATEST_ALLOWLIST_FILE = joinpath("config", "ci", "models_invokelatest_allowlist.toml")
const PNJL_SCAN_DEFAULT_INCLUDE_ALLOWLIST = Set(["ScanConfig.jl"])

struct GateConfig
    mode::Symbol
    base::Union{Nothing, String}
    head::Union{Nothing, String}
    unit_profile::String
    auto_rollback::Bool
end

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

normalize_path(path::AbstractString) = replace(path, '\\' => '/')

function parse_args(args::Vector{String})
    staged = false
    base = nothing
    head = nothing
    unit_profile = "smoke"
    auto_rollback = false

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
        elseif arg == "--unit-profile"
            i < length(args) || error("--unit-profile requires a value")
            unit_profile = args[i + 1]
            i += 2
        elseif arg == "--auto-rollback"
            auto_rollback = true
            i += 1
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
    return GateConfig(mode, base, head, unit_profile, auto_rollback)
end

function git_read_lines(root::String, args::Vector{String})
    cmd = Cmd(["git", "-C", root, args...])
    out = read(cmd, String)
    return isempty(out) ? String[] : split(chomp(out), '\n')
end

function target_name_status(root::String, cfg::GateConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--name-status", "--", "src/pnjl"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--name-status", range, "--", "src/pnjl"])
    else
        return git_read_lines(root, ["diff", "--name-status", "HEAD", "--", "src/pnjl"])
    end
end

function target_patch_allowlist(root::String, cfg::GateConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", MODELS_INVOKELATEST_ALLOWLIST_FILE])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", MODELS_INVOKELATEST_ALLOWLIST_FILE])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", MODELS_INVOKELATEST_ALLOWLIST_FILE])
    end
end

function collect_allowlist_delta(root::String, cfg::GateConfig)
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

function collect_deleted_pnjl_files(root::String, cfg::GateConfig)
    deleted = String[]
    for line in target_name_status(root, cfg)
        isempty(line) && continue
        fields = split(line, '\t')
        length(fields) >= 2 || continue
        status = fields[1]
        path = normalize_path(fields[end])
        if startswith(status, "D")
            push!(deleted, path)
        end
    end
    return sort!(unique(deleted))
end

function result_timestamp()
    return Dates.format(now(), "yyyymmdd_HHMMSS")
end

function snapshot_path(root::String, ts::String)
    out_dir = joinpath(root, "data", "outputs", "results")
    mkpath(out_dir)
    return joinpath(out_dir, "pnjl_prune_wave_snapshot_$(ts).txt")
end

function allowlist_delta_path(root::String, ts::String)
    out_dir = joinpath(root, "data", "outputs", "results")
    mkpath(out_dir)
    return joinpath(out_dir, "models_invokelatest_allowlist_delta_$(ts).txt")
end

function scan_default_include_audit_path(root::String, ts::String)
    out_dir = joinpath(root, "data", "outputs", "results")
    mkpath(out_dir)
    return joinpath(out_dir, "pnjl_scan_default_include_audit_$(ts).txt")
end

function collect_scan_default_includes(root::String)
    pnjl_path = joinpath(root, "src", "pnjl", "PNJL.jl")
    isfile(pnjl_path) || return String[]

    hits = String[]
    line_no = 0
    for line in eachline(pnjl_path)
        line_no += 1
        startswith(line, "include(") || continue
        s = strip(line)
        m = match(r"^include\(joinpath\(@__DIR__,\s*\"scans\",\s*\"([^\"]+\.jl)\"\)\)$", s)
        m === nothing && continue
        name = m.captures[1]
        name in PNJL_SCAN_DEFAULT_INCLUDE_ALLOWLIST && continue
        push!(hits, "src/pnjl/PNJL.jl:$(line_no): $(s)")
    end
    return sort!(unique(hits))
end

function write_snapshot(path::String, cfg::GateConfig, deleted::Vector{String}, allowlist_delta)
    open(path, "w") do io
        println(io, "# PNJL prune-wave snapshot")
        println(io, "timestamp = ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "mode = ", cfg.mode)
        if cfg.mode === :range
            println(io, "range = ", cfg.base, "...", cfg.head)
        end
        println(io, "deleted_count = ", length(deleted))
        println(io, "allowlist_added_count = ", length(allowlist_delta.added))
        println(io, "allowlist_removed_count = ", length(allowlist_delta.removed))
        println(io)
        for p in deleted
            println(io, p)
        end

        if !isempty(allowlist_delta.added)
            println(io)
            println(io, "[allowlist_added]")
            for item in allowlist_delta.added
                println(io, item)
            end
        end

        if !isempty(allowlist_delta.removed)
            println(io)
            println(io, "[allowlist_removed]")
            for item in allowlist_delta.removed
                println(io, item)
            end
        end
    end
end

function write_allowlist_delta(path::String, cfg::GateConfig, allowlist_delta)
    open(path, "w") do io
        println(io, "# models invokelatest allowlist delta")
        println(io, "timestamp = ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "mode = ", cfg.mode)
        if cfg.mode === :range
            println(io, "range = ", cfg.base, "...", cfg.head)
        end
        println(io, "allowlist_file = ", MODELS_INVOKELATEST_ALLOWLIST_FILE)
        println(io, "added_count = ", length(allowlist_delta.added))
        println(io, "removed_count = ", length(allowlist_delta.removed))

        println(io)
        println(io, "[added]")
        for item in allowlist_delta.added
            println(io, item)
        end

        println(io)
        println(io, "[removed]")
        for item in allowlist_delta.removed
            println(io, item)
        end
    end
end

function write_scan_default_include_audit(path::String, cfg::GateConfig, hits::Vector{String})
    open(path, "w") do io
        println(io, "# pnjl scan default include audit")
        println(io, "timestamp = ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "mode = ", cfg.mode)
        if cfg.mode === :range
            println(io, "range = ", cfg.base, "...", cfg.head)
        end
        println(io, "scope = src/pnjl/PNJL.jl scans/* (excluding allowlist)")
        println(io, "allowlist = ", join(sort!(collect(PNJL_SCAN_DEFAULT_INCLUDE_ALLOWLIST)), ","))
        println(io, "observed = ", length(hits))
        println(io)
        for item in hits
            println(io, item)
        end
    end
end

function guard_cmd(root::String, cfg::GateConfig)
    if cfg.mode === :staged
        return Cmd(Cmd(["julia", "--project=.", joinpath("scripts", "dev", "check_pnjl_migration_guard.jl"), "--staged"]); dir=root)
    elseif cfg.mode === :range
        return Cmd(Cmd(["julia", "--project=.", joinpath("scripts", "dev", "check_pnjl_migration_guard.jl"), "--base", cfg.base, "--head", cfg.head]); dir=root)
    else
        return Cmd(Cmd(["julia", "--project=.", joinpath("scripts", "dev", "check_pnjl_migration_guard.jl")]); dir=root)
    end
end

function output_path_guard_cmd(root::String, cfg::GateConfig)
    if cfg.mode === :staged
        return Cmd(Cmd(["julia", "--project=.", joinpath("scripts", "dev", "check_data_output_path_guard.jl"), "--staged"]); dir=root)
    elseif cfg.mode === :range
        return Cmd(Cmd(["julia", "--project=.", joinpath("scripts", "dev", "check_data_output_path_guard.jl"), "--base", cfg.base, "--head", cfg.head]); dir=root)
    else
        return Cmd(Cmd(["julia", "--project=.", joinpath("scripts", "dev", "check_data_output_path_guard.jl")]); dir=root)
    end
end

function smoke_cmd(root::String)
    return Cmd(Cmd(["julia", "--project=.", joinpath("tests", "unit", "runtests.jl")]); dir=root)
end

function rollback_cmd(root::String)
    return Cmd(["git", "-C", root, "restore", "--worktree", "--staged", "--source=HEAD", "--", "src/pnjl"])
end

function main()
    cfg = parse_args(ARGS)
    root = project_root()

    ts = result_timestamp()
    deleted = collect_deleted_pnjl_files(root, cfg)
    allowlist_delta = collect_allowlist_delta(root, cfg)
    scan_default_include_hits = collect_scan_default_includes(root)
    snap = snapshot_path(root, ts)
    allowlist_audit = allowlist_delta_path(root, ts)
    scan_default_include_audit = scan_default_include_audit_path(root, ts)
    write_snapshot(snap, cfg, deleted, allowlist_delta)
    write_allowlist_delta(allowlist_audit, cfg, allowlist_delta)
    write_scan_default_include_audit(scan_default_include_audit, cfg, scan_default_include_hits)

    println("[pnjl-prune-wave-gate] snapshot: ", normalize_path(relpath(snap, root)))
    println("[pnjl-prune-wave-gate] allowlist audit: ", normalize_path(relpath(allowlist_audit, root)))
    println("[pnjl-prune-wave-gate] scan default include audit: ", normalize_path(relpath(scan_default_include_audit, root)))
    println("[pnjl-prune-wave-gate] scan default include observed: ", length(scan_default_include_hits))
    println("[pnjl-prune-wave-gate] deleted files: ", length(deleted))
    if !isempty(allowlist_delta.added) || !isempty(allowlist_delta.removed)
        println("[pnjl-prune-wave-gate] allowlist diff: added=", length(allowlist_delta.added), " removed=", length(allowlist_delta.removed))
        for item in Iterators.take(allowlist_delta.added, 3)
            println("  + ", item)
        end
        for item in Iterators.take(allowlist_delta.removed, 3)
            println("  - ", item)
        end
    end

    println("[pnjl-prune-wave-gate] run migration guard ...")
    run(guard_cmd(root, cfg))

    println("[pnjl-prune-wave-gate] run data output path guard ...")
    run(output_path_guard_cmd(root, cfg))

    println("[pnjl-prune-wave-gate] run unit smoke ...")
    withenv("UNIT_PROFILE" => cfg.unit_profile) do
        run(smoke_cmd(root))
    end

    println("[pnjl-prune-wave-gate] OK (profile=$(cfg.unit_profile))")
    return
end

try
    main()
catch err
    println("[pnjl-prune-wave-gate] FAILED: ", sprint(showerror, err))

    cfg = try
        parse_args(ARGS)
    catch
        nothing
    end

    if cfg !== nothing && cfg.mode === :working && cfg.auto_rollback
        root = project_root()
        println("[pnjl-prune-wave-gate] auto rollback: git restore --worktree --staged --source=HEAD -- src/pnjl")
        try
            run(rollback_cmd(root))
            println("[pnjl-prune-wave-gate] rollback done")
        catch rb_err
            println("[pnjl-prune-wave-gate] rollback failed: ", sprint(showerror, rb_err))
        end
    else
        println("[pnjl-prune-wave-gate] rollback hint: git restore --worktree --staged --source=HEAD -- src/pnjl")
    end

    exit(1)
end
