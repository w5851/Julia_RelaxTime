#!/usr/bin/env julia

"""
数据输出路径门禁（scripts/*.jl）

目的：禁止新增默认输出路径回流到根目录 outputs/results，统一要求 data/outputs/results 分层。

检查范围：
- 默认 working tree vs HEAD
- --staged
- --base <rev> --head <rev>

判定规则（仅检查新增行）：
1) 出现 `joinpath(..., "outputs", "results", ...)` 且同一行未包含 `"data"`。
2) 出现文本 `outputs/results` 且不包含 `data/outputs/results`。
"""

struct GuardConfig
    mode::Symbol
    base::Union{Nothing, String}
    head::Union{Nothing, String}
end

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

function is_forbidden_output_line(line::String)
    s = strip(line)

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

function main()
    cfg = parse_args(ARGS)
    root = project_root()

    violations = collect_violations(root, cfg)

    if !isempty(violations)
        println("[data-output-path-guard] FAILED")
        println("  Detected new outputs/results paths outside data/outputs/results:")
        for item in violations
            println("   - " * item)
        end
        println("  hint: use data/outputs/results/<domain>/<category>/... as default output path")
        exit(1)
    end

    println("[data-output-path-guard] OK")
    if cfg.mode === :staged
        println("  mode=staged")
    elseif cfg.mode === :range
        println("  mode=range ($(cfg.base)...$(cfg.head))")
    else
        println("  mode=working-tree-vs-HEAD")
    end
end

main()
