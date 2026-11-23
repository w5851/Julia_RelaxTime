#!/usr/bin/env julia

"""
数据输出路径门禁（scripts/*.jl）

目的：
1) 禁止新增默认输出路径回流到根目录 outputs/results，统一要求 data/outputs/results 分层。
2) 禁止新增正式图像资产写入 data/outputs/results；图像应写入 data/outputs/figures。
3) 禁止新增非图像/非 plot_manifest 资产写入 data/outputs/figures。

检查范围：
- 默认 working tree vs HEAD
- --staged
- --base <rev> --head <rev>

判定规则（仅检查新增行）：
1) 出现 `joinpath(..., "outputs", "results", ...)` 且同一行未包含 `"data"`。
2) 出现文本 `outputs/results` 且不包含 `data/outputs/results`。
3) 新增/复制/重命名 tracked 图像文件位于 `data/outputs/results/**.{png,svg,pdf}`。
   已验收的 immutable charged GBU v2 bundle 是唯一保留的哈希绑定历史包例外；
   新生产运行仍必须把图像写入 `data/outputs/figures/`。
4) 新增/复制/重命名 tracked 非图像文件位于 `data/outputs/figures/**`，且不是 `plot_manifest.json`。
5) 归档策略（config/governance/data_archive_policy.toml）：新增/改名的 `data/outputs/**` 下 .zip
   必须登记在 retained_archives 或 pending_retirement；删除必须登记在 retired_archives；
   任何未被登记的 tracked .zip 都视为 violation（hard fail）；
   已登记条目的 sha256 必须与工作区文件一致，retired 条目不得仍在工作区存在，
   且其 git_blob 必须仍在对象库中且哈希回墓碑记录的 sha256。
"""

import SHA
import TOML

struct GuardConfig
    mode::Symbol
    base::Union{Nothing, String}
    head::Union{Nothing, String}
end

const IMMUTABLE_RESULT_BUNDLE_PREFIXES = (
    "data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/freezeout_20260907_v2/",
)

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
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

function target_patch_scripts(root::String, cfg::GuardConfig)
    if cfg.mode === :staged
        return git_read_lines(root, ["diff", "--cached", "--unified=0", "--", "scripts"])
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        return git_read_lines(root, ["diff", "--unified=0", range, "--", "scripts"])
    else
        return git_read_lines(root, ["diff", "--unified=0", "HEAD", "--", "scripts"])
    end
end

function target_changed_result_files(root::String, cfg::GuardConfig)
    args = if cfg.mode === :staged
        ["diff", "--cached", "--name-status", "--diff-filter=ACMR", "--", "data/outputs/results"]
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        ["diff", "--name-status", "--diff-filter=ACMR", range, "--", "data/outputs/results"]
    else
        ["diff", "--name-status", "--diff-filter=ACMR", "HEAD", "--", "data/outputs/results"]
    end
    return git_read_lines(root, args)
end

function target_changed_figure_files(root::String, cfg::GuardConfig)
    args = if cfg.mode === :staged
        ["diff", "--cached", "--name-status", "--diff-filter=ACMR", "--", "data/outputs/figures"]
    elseif cfg.mode === :range
        range = "$(cfg.base)...$(cfg.head)"
        ["diff", "--name-status", "--diff-filter=ACMR", range, "--", "data/outputs/figures"]
    else
        ["diff", "--name-status", "--diff-filter=ACMR", "HEAD", "--", "data/outputs/figures"]
    end
    return git_read_lines(root, args)
end

function is_forbidden_output_line(line::AbstractString)
    s = String(strip(line))

    joinpath_outputs_results = occursin(r"joinpath\([^\n]*\"outputs\"\s*,\s*\"results\"", s)
    if joinpath_outputs_results && !occursin("\"data\"", s)
        return true
    end

    textual_outputs_results = occursin("outputs/results", s)
    if textual_outputs_results && !occursin("data/outputs/results", s)
        return true
    end

    return false
end

function changed_path_from_name_status(line::AbstractString)
    parts = split(String(line), '\t')
    length(parts) >= 2 || return nothing
    status = parts[1]
    if startswith(status, "R") || startswith(status, "C")
        length(parts) >= 3 || return nothing
        return normalize_path(parts[end])
    end
    return normalize_path(parts[2])
end

function is_result_figure_path(path::AbstractString)
    normalized = normalize_path(path)
    startswith(normalized, "data/outputs/results/") || return false
    any(startswith(normalized, prefix) for prefix in IMMUTABLE_RESULT_BUNDLE_PREFIXES) && return false
    ext = lowercase(splitext(normalized)[2])
    return ext in (".png", ".svg", ".pdf")
end

function is_allowed_figure_asset_path(path::AbstractString)
    normalized = normalize_path(path)
    startswith(normalized, "data/outputs/figures/") || return true
    basename(normalized) == "plot_manifest.json" && return true
    ext = lowercase(splitext(normalized)[2])
    return ext in (".png", ".svg", ".pdf")
end

function collect_violations(root::String, cfg::GuardConfig)
    violations = String[]

    current_file = ""
    new_lineno = 0

    for line in target_patch_scripts(root, cfg)
        if startswith(line, "+++ b/")
            current_file = normalize_path(line[7:end])
            new_lineno = 0
            continue
        elseif startswith(line, "@@")
            m = match(r"\+(\d+)", line)
            if m !== nothing
                new_lineno = parse(Int, m.captures[1]) - 1
            end
            continue
        end

        isempty(current_file) && continue
        endswith(current_file, ".jl") || continue

        if startswith(line, "+") && !startswith(line, "+++")
            new_lineno += 1
            added = line[2:end]
            startswith(strip(added), "#") && continue
            if is_forbidden_output_line(added)
                push!(violations, "$(current_file):$(new_lineno): $(strip(added))")
            end
        elseif startswith(line, "-") && !startswith(line, "---")
            continue
        else
            new_lineno += 1
        end
    end

    return violations
end

function collect_result_figure_violations(root::String, cfg::GuardConfig)
    violations = String[]
    for line in target_changed_result_files(root, cfg)
        path = changed_path_from_name_status(line)
        path === nothing && continue
        if is_result_figure_path(path)
            push!(violations, path)
        end
    end
    return violations
end

function collect_figure_asset_violations(root::String, cfg::GuardConfig)
    violations = String[]
    for line in target_changed_figure_files(root, cfg)
        path = changed_path_from_name_status(line)
        path === nothing && continue
        if !is_allowed_figure_asset_path(path)
            push!(violations, path)
        end
    end
    return violations
end

const ARCHIVE_POLICY_RELATIVE_PATH = "config/governance/data_archive_policy.toml"
const ARCHIVE_POLICY_SCHEMA = "data_archive_policy_v1"
const ARCHIVE_OUTPUT_PREFIX = "data/outputs/"

const RETAINED_REQUIRED_FIELDS = ("path", "sha256", "entry_count", "rationale", "approved_at")
const RETIRED_REQUIRED_FIELDS = (
    "path", "sha256", "entry_count", "declared_commit", "uncommitted_source",
    "snapshot_mtime", "replacement", "rationale", "retired_at", "git_blob", "retired_before_commit",
)
const PENDING_REQUIRED_FIELDS = ("path", "sha256", "entry_count", "reason", "planned_phase")
const ALLOWED_RETAINED_STATUS = ("active", "historical_superseded", "historical_evidence")

struct ArchivePolicy
    retained::Vector{Dict{String, Any}}
    retired::Vector{Dict{String, Any}}
    pending::Vector{Dict{String, Any}}
end

function load_archive_policy(root::AbstractString)
    path = joinpath(root, ARCHIVE_POLICY_RELATIVE_PATH)
    isfile(path) || return nothing
    data = TOML.parsefile(path)
    schema = get(data, "schema", "")
    schema == ARCHIVE_POLICY_SCHEMA ||
        error("archive policy schema mismatch in $(ARCHIVE_POLICY_RELATIVE_PATH): got $(repr(schema))")
    to_entries(key) = [
        Dict{String, Any}(String(k) => v for (k, v) in entry)
        for entry in get(data, key, [])
    ]
    return ArchivePolicy(to_entries("retained_archives"), to_entries("retired_archives"), to_entries("pending_retirement"))
end

is_sha256_hex(text) = text isa AbstractString && occursin(r"^[0-9a-f]{64}$", text)
is_full_git_id(text) = text isa AbstractString && occursin(r"^[0-9a-f]{40}$", text)
policy_path(entry) = normalize_path(string(get(entry, "path", "")))

function archive_policy_entry_violations(entry::Dict{String, Any}, required, label::AbstractString)
    violations = String[]
    path = policy_path(entry)
    if isempty(path)
        push!(violations, "$(label): entry without path")
        return violations
    end
    startswith(path, ARCHIVE_OUTPUT_PREFIX) ||
        push!(violations, "$(path): archive path must live under $(ARCHIVE_OUTPUT_PREFIX)")
    endswith(lowercase(path), ".zip") ||
        push!(violations, "$(path): archive entries must be .zip files")
    for field in required
        haskey(entry, field) || push!(violations, "$(path): missing field $(field)")
    end
    if haskey(entry, "sha256") && !is_sha256_hex(entry["sha256"])
        push!(violations, "$(path): sha256 must be 64 lowercase hex characters")
    end
    if haskey(entry, "entry_count")
        count = entry["entry_count"]
        (count isa Integer && count > 0) ||
            push!(violations, "$(path): entry_count must be a positive integer")
    end
    if haskey(entry, "status") && !(entry["status"] in ALLOWED_RETAINED_STATUS)
        push!(
            violations,
            "$(path): status must be one of $(join(ALLOWED_RETAINED_STATUS, " / ")), got $(repr(entry["status"]))",
        )
    end
    return violations
end

function changed_archive_paths(root::AbstractString, cfg::GuardConfig, filter::AbstractString)
    args = if cfg.mode === :staged
        ["diff", "--cached", "--name-status", "--diff-filter=$(filter)", "--", "data/outputs"]
    elseif cfg.mode === :range
        ["diff", "--name-status", "--diff-filter=$(filter)", "$(cfg.base)...$(cfg.head)", "--", "data/outputs"]
    else
        ["diff", "--name-status", "--diff-filter=$(filter)", "HEAD", "--", "data/outputs"]
    end
    paths = String[]
    for line in git_read_lines(root, args)
        path = changed_path_from_name_status(line)
        path === nothing && continue
        endswith(lowercase(path), ".zip") || continue
        push!(paths, path)
    end
    return paths
end

function tracked_archive_paths(root::AbstractString)
    files = git_read_lines(root, ["ls-files", "--", "data/outputs"])
    return [normalize_path(file) for file in files if endswith(lowercase(file), ".zip")]
end

file_sha256(path::AbstractString) = bytes2hex(SHA.sha256(read(path)))

function is_shallow_repository(root::AbstractString)
    return strip(read(Cmd(["git", "-C", root, "rev-parse", "--is-shallow-repository"]), String)) == "true"
end

function git_blob_bytes(root::AbstractString, blob::AbstractString)
    buffer = IOBuffer()
    run(pipeline(Cmd(["git", "-C", root, "cat-file", "blob", blob]); stdout = buffer, stderr = devnull))
    return take!(buffer)
end

"""
Verify that every tombstone's `git_blob` really is in the object database and that its bytes
still hash to the recorded `sha256`. Skipped on shallow clones, where historical blobs are absent.
"""
function collect_retired_blob_violations(root::AbstractString, policy::ArchivePolicy)
    violations = String[]
    isempty(policy.retired) && return violations
    is_shallow_repository(root) && return violations

    for entry in policy.retired
        blob = string(get(entry, "git_blob", ""))
        is_full_git_id(blob) || continue  # format violations are reported elsewhere
        path = policy_path(entry)
        try
            actual = bytes2hex(SHA.sha256(git_blob_bytes(root, blob)))
            expected = string(get(entry, "sha256", ""))
            actual == expected ||
                push!(violations, "$(path): recorded git_blob $(blob) hashes to $(actual), not the tombstone sha256 $(expected)")
        catch err
            push!(violations, "$(path): recorded git_blob $(blob) is missing from the local object database")
        end
    end
    return violations
end

function collect_retired_blob_advisories(root::AbstractString, policy::ArchivePolicy)
    isempty(policy.retired) && return String[]
    is_shallow_repository(root) || return String[]
    return ["shallow clone detected: $(length(policy.retired)) tombstones were not verified against the local object database"]
end

function registered_path_set(entries)
    return Set(policy_path(entry) for entry in entries)
end

function collect_archive_violations(root::AbstractString, cfg::GuardConfig, policy::ArchivePolicy)
    violations = String[]

    retained_paths = registered_path_set(policy.retained)
    pending_paths = registered_path_set(policy.pending)
    retired_paths = registered_path_set(policy.retired)

    for path in changed_archive_paths(root, cfg, "ACR")
        if !(path in retained_paths) && !(path in pending_paths)
            push!(violations, "$(path): new or renamed archive is not registered in $(ARCHIVE_POLICY_RELATIVE_PATH)")
        end
    end

    for path in changed_archive_paths(root, cfg, "D")
        if !(path in retired_paths)
            push!(violations, "$(path): deleted archive requires a retired_archives tombstone in $(ARCHIVE_POLICY_RELATIVE_PATH)")
        end
    end

    for (label, entries, required) in (
        ("retained_archives", policy.retained, RETAINED_REQUIRED_FIELDS),
        ("pending_retirement", policy.pending, PENDING_REQUIRED_FIELDS),
    )
        for entry in entries
            append!(violations, archive_policy_entry_violations(entry, required, label))
            path = policy_path(entry)
            isempty(path) && continue
            full = joinpath(root, path)
            if isfile(full)
                expected = get(entry, "sha256", "")
                actual = file_sha256(full)
                actual == expected ||
                    push!(violations, "$(path): sha256 mismatch (recorded $(expected), actual $(actual))")
            else
                push!(violations, "$(path): registered archive is missing from the working tree")
            end
        end
    end

    for entry in policy.retired
        append!(violations, archive_policy_entry_violations(entry, RETIRED_REQUIRED_FIELDS, "retired_archives"))
        path = policy_path(entry)
        isempty(path) && continue
        isfile(joinpath(root, path)) &&
            push!(violations, "$(path): retired archive is still present in the working tree")
        if haskey(entry, "git_blob") && !is_full_git_id(entry["git_blob"])
            push!(violations, "$(path): git_blob must be a 40 character git object id")
        end
        if haskey(entry, "retired_before_commit") && !is_full_git_id(entry["retired_before_commit"])
            push!(violations, "$(path): retired_before_commit must be a 40 character git object id")
        end
    end

    append!(violations, collect_retired_blob_violations(root, policy))

    for path in tracked_archive_paths(root)
        if !(path in retained_paths) && !(path in pending_paths) && !(path in retired_paths)
            push!(violations, "$(path): tracked archive under data/outputs is not registered in $(ARCHIVE_POLICY_RELATIVE_PATH)")
        end
    end

    return violations
end

function collect_archive_advisories(root::AbstractString, policy::ArchivePolicy)
    return collect_retired_blob_advisories(root, policy)
end

function main()
    cfg = parse_args(ARGS)
    root = project_root()

    violations = collect_violations(root, cfg)
    figure_violations = collect_result_figure_violations(root, cfg)
    figure_asset_violations = collect_figure_asset_violations(root, cfg)

    policy = load_archive_policy(root)
    archive_violations = String[]
    archive_advisories = String[]
    if policy !== nothing
        archive_violations = collect_archive_violations(root, cfg, policy)
        archive_advisories = collect_archive_advisories(root, policy)
    end

    if !isempty(violations) || !isempty(figure_violations) || !isempty(figure_asset_violations) || !isempty(archive_violations)
        println("[data-output-path-guard] FAILED")
        if !isempty(violations)
            println("  Detected new outputs/results paths outside data/outputs/results:")
        end
        for item in violations
            println("   - " * item)
        end
        if !isempty(figure_violations)
            println("  Detected figure artifacts under data/outputs/results:")
            for item in figure_violations
                println("   - " * item)
            end
        end
        if !isempty(figure_asset_violations)
            println("  Detected non-figure artifacts under data/outputs/figures:")
            for item in figure_asset_violations
                println("   - " * item)
            end
        end
        if !isempty(archive_violations)
            println("  Detected archive policy violations ($(ARCHIVE_POLICY_RELATIVE_PATH)):")
            for item in archive_violations
                println("   - " * item)
            end
        end
        println("  hint: use data/outputs/results/<domain>/<category>/... as default output path")
        println("  hint: use data/outputs/figures/<domain>/<category>/... for PNG/SVG/PDF artifacts and plot_manifest.json only")
        println("  hint: register archives in $(ARCHIVE_POLICY_RELATIVE_PATH); retired archives need a tombstone with git_blob")
        exit(1)
    end

    println("[data-output-path-guard] OK")
    for item in archive_advisories
        println("  [ADVISORY] " * item)
    end
    if cfg.mode === :staged
        println("  mode=staged")
    elseif cfg.mode === :range
        println("  mode=range ($(cfg.base)...$(cfg.head))")
    else
        println("  mode=working-tree-vs-HEAD")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
