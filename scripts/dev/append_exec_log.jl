#!/usr/bin/env julia

using Dates

const ROOT = pwd()
const ACTIVE_DIR = joinpath(ROOT, "docs", "dev", "active")

function print_usage()
    println("""
Usage: julia --project=. scripts/dev/append_exec_log.jl [OPTIONS]

Append one execution record block to an execution log file (append-only).
By default this script does NOT read log history body.

Options:
  -h, --help                    Show this help message
  --log-file PATH               Target log file path (absolute or relative)
  --task-file PATH              Task file path, used to auto-derive log file if --log-file is absent
  --batch TEXT                  Batch id, e.g. \"Batch N2\" (required)
  --goal TEXT                   Goal of this batch (required)
  --code-change TEXT            Code change summary (default: empty)
  --cmd TEXT                    Validation command (repeatable)
  --artifact TEXT               Artifact path (repeatable)
  --result TEXT                 Result summary (required)
  --mainline TEXT               Mainline mapping, e.g. N1/N2/N3/N4 (required)
  --notes TEXT                  Optional notes
  --dry-run                     Print append block without writing file

Examples:
  julia --project=. scripts/dev/append_exec_log.jl \\
    --task-file docs/dev/active/2026-02-26_多重派发重构与PNJL迁移下线开发任务单.md \\
    --batch \"Batch N2\" \\
    --goal \"移除旧入口并完成冒烟校验\" \\
    --code-change \"删除 src/pnjl 兼容路径调用点\" \\
    --cmd \"julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD\" \\
    --cmd \"UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl\" \\
    --artifact \"data/outputs/results/pnjl_prune_wave_snapshot_20260226_101500.txt\" \\
    --result \"通过\" \\
    --mainline \"N2\"
""")
end

normalize_sep(path::String) = replace(path, '\\' => '/')

function abs_path(path::String)
    if isabspath(path)
        return normpath(path)
    end
    return normpath(joinpath(ROOT, path))
end

function rel_to_root(path::String)
    ap = abs_path(path)
    rootp = abs_path(ROOT)
    ap_norm = normalize_sep(ap)
    root_norm = normalize_sep(rootp)
    if startswith(ap_norm, root_norm * "/")
        return ap_norm[length(root_norm)+2:end]
    elseif ap_norm == root_norm
        return "."
    else
        return ap_norm
    end
end

function discover_latest_task_file()
    if !isdir(ACTIVE_DIR)
        return nothing
    end
    candidates = String[]
    for file in readdir(ACTIVE_DIR)
        if endswith(file, ".md") && occursin("开发任务单", file)
            push!(candidates, joinpath(ACTIVE_DIR, file))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates)
    return candidates[end]
end

function derive_log_file_from_task(task_file::String)
    filename = basename(task_file)
    if occursin("开发任务单", filename)
        log_filename = replace(filename, "开发任务单" => "执行台账")
    elseif endswith(filename, ".md")
        log_filename = replace(filename, ".md" => "_执行台账.md")
    else
        log_filename = filename * "_执行台账.md"
    end
    return joinpath(dirname(task_file), log_filename)
end

function infer_title_from_filename(log_file::String)
    filename = basename(log_file)
    base = replace(filename, r"\.md$" => "")
    base = replace(base, r"^\d{4}[-_]\d{2}[-_]\d{2}_" => "")
    base = replace(base, "_" => " ")
    return isempty(strip(base)) ? "执行台账" : strip(base)
end

function bootstrap_log_file(log_file::String, task_file::Union{Nothing, String})
    mkpath(dirname(log_file))
    title = infer_title_from_filename(log_file)
    today_str = Dates.format(today(), "yyyy-mm-dd")

    task_line = if task_file === nothing
        "- 本阶段开发任务定义统一维护在：\n  - （未指定）"
    else
        "- 本阶段开发任务定义统一维护在：\n  - `$(rel_to_root(task_file))`"
    end

    content = """
# $(title)

更新日期：$(today_str)

---

## 0. 台账定位

- 本文档仅记录“执行事实”：变更点、命令、产物、结果。
$(task_line)

---

## 1. 记录规范（强约束）

- [x] 执行记录与开发任务分开保存：
  - 任务目标/范围/DoD 只写“开发任务单”；
  - 执行过程/命令/结果只写“执行台账”。
- [x] 追加记录时采用“直接追加”策略：
  - 不要求回读本台账历史上下文；
  - 仅按统一模板追加到文档末尾；
  - 每条记录必须自包含（含目标、命令、产物、结果）。
- [x] 每条记录必须可追溯到输出产物（`data/outputs/results/*`）。

---

## 2. 执行记录
"""

    open(log_file, "w") do io
        write(io, content)
    end
end

function render_entry(
    batch::String,
    goal::String,
    code_change::String,
    commands::Vector{String},
    artifacts::Vector{String},
    result::String,
    mainline::String,
    notes::String,
)
    cmds = isempty(commands) ? ["（无）"] : commands
    arts = isempty(artifacts) ? ["（无）"] : artifacts

    cmd_lines = join(["    - `$(c)`" for c in cmds], "\n")
    art_lines = join(["    - `$(a)`" for a in arts], "\n")

    notes_line = isempty(strip(notes)) ? "" : "\n  - 备注：$(notes)"

    return """

- [ ] 批次号：$(batch)
  - 目标：$(goal)
  - 代码变更：$(isempty(strip(code_change)) ? "（无）" : code_change)
  - 验证命令：
$(cmd_lines)
  - 产物：
$(art_lines)
  - 结果：$(result)
  - 主线映射（N1/N2/N3/N4）：$(mainline)$(notes_line)
"""
end

function parse_args(args::Vector{String})
    opts = Dict{String, Any}(
        "log_file" => nothing,
        "task_file" => nothing,
        "batch" => "",
        "goal" => "",
        "code_change" => "",
        "commands" => String[],
        "artifacts" => String[],
        "result" => "",
        "mainline" => "",
        "notes" => "",
        "dry_run" => false,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ["-h", "--help"]
            print_usage()
            return nothing
        elseif arg == "--dry-run"
            opts["dry_run"] = true
            i += 1
        elseif arg in ["--log-file", "--task-file", "--batch", "--goal", "--code-change", "--cmd", "--artifact", "--result", "--mainline", "--notes"]
            if i == length(args)
                error("Missing value for option: $(arg)")
            end
            val = args[i + 1]
            if arg == "--log-file"
                opts["log_file"] = val
            elseif arg == "--task-file"
                opts["task_file"] = val
            elseif arg == "--batch"
                opts["batch"] = val
            elseif arg == "--goal"
                opts["goal"] = val
            elseif arg == "--code-change"
                opts["code_change"] = val
            elseif arg == "--cmd"
                push!(opts["commands"], val)
            elseif arg == "--artifact"
                push!(opts["artifacts"], val)
            elseif arg == "--result"
                opts["result"] = val
            elseif arg == "--mainline"
                opts["mainline"] = val
            elseif arg == "--notes"
                opts["notes"] = val
            end
            i += 2
        else
            error("Unknown option: $(arg). Use --help for usage.")
        end
    end

    return opts
end

function resolve_log_and_task(log_file_opt, task_file_opt)
    task_file = task_file_opt === nothing ? nothing : abs_path(task_file_opt)

    if log_file_opt !== nothing
        return abs_path(log_file_opt), task_file
    end

    if task_file === nothing
        task_file = discover_latest_task_file()
    end

    if task_file !== nothing
        return derive_log_file_from_task(task_file), task_file
    end

    fallback = joinpath(ACTIVE_DIR, "$(Dates.format(today(), "yyyy-mm-dd"))_执行台账.md")
    return fallback, nothing
end

function validate_required(opts)
    required = ["batch", "goal", "result", "mainline"]
    missing = String[]
    for k in required
        if isempty(strip(String(opts[k])))
            push!(missing, k)
        end
    end
    if !isempty(missing)
        error("Missing required options: " * join(missing, ", "))
    end
end

function main()
    if isempty(ARGS)
        print_usage()
        exit(1)
    end

    opts = parse_args(ARGS)
    opts === nothing && return

    validate_required(opts)

    log_file, task_file = resolve_log_and_task(opts["log_file"], opts["task_file"])

    entry = render_entry(
        String(opts["batch"]),
        String(opts["goal"]),
        String(opts["code_change"]),
        Vector{String}(opts["commands"]),
        Vector{String}(opts["artifacts"]),
        String(opts["result"]),
        String(opts["mainline"]),
        String(opts["notes"]),
    )

    if Bool(opts["dry_run"])
        println("[append-exec-log] DRY RUN")
        println("  log_file: $(rel_to_root(log_file))")
        println("  task_file: $(task_file === nothing ? "(auto: none)" : rel_to_root(task_file))")
        println(entry)
        return
    end

    if !isfile(log_file)
        bootstrap_log_file(log_file, task_file)
        println("[append-exec-log] created log file: $(rel_to_root(log_file))")
    end

    open(log_file, "a") do io
        write(io, entry)
    end

    println("[append-exec-log] appended")
    println("  log_file: $(rel_to_root(log_file))")
    println("  batch: $(opts["batch"])")
    println("  commands: $(length(opts["commands"]))")
    println("  artifacts: $(length(opts["artifacts"]))")
end

main()
