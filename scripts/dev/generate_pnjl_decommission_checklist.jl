#!/usr/bin/env julia

using Dates

struct ChecklistConfig
    base::String
    head::String
end

const SCAN_DEFAULT_INCLUDE_ALLOWLIST = Set(["ScanConfig.jl"])

normalize_path(path::AbstractString) = replace(path, '\\' => '/')

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function parse_args(args::Vector{String})
    base = "HEAD"
    head = "HEAD"

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--base"
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

    return ChecklistConfig(base, head)
end

function git_lines(root::String, args::Vector{String})
    cmd = Cmd(["git", "-C", root, args...])
    out = read(cmd, String)
    return isempty(out) ? String[] : split(chomp(out), '\n')
end

function grep_lines(root::String, pattern::String, paths::Vector{String})
    cmd = ["grep", "-n", "-E", pattern, "--"]
    append!(cmd, paths)
    try
        return git_lines(root, cmd)
    catch
        return String[]
    end
end

function collect_scan_files(root::String)
    scans_dir = joinpath(root, "src", "models", "scans")
    if !isdir(scans_dir)
        return String[]
    end
    files = String[]
    for name in readdir(scans_dir)
        endswith(name, ".jl") || continue
        push!(files, normalize_path(joinpath("src", "models", "scans", name)))
    end
    return sort(files)
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
        name in SCAN_DEFAULT_INCLUDE_ALLOWLIST && continue
        push!(hits, "src/pnjl/PNJL.jl:$(line_no): $(s)")
    end
    return unique_sorted(hits)
end

function collect_scan_loaders(root::String)
    pnjl_path = joinpath(root, "src", "pnjl", "PNJL.jl")
    isfile(pnjl_path) || return String[]

    hits = String[]
    line_no = 0
    for line in eachline(pnjl_path)
        line_no += 1
        s = strip(line)
        startswith(s, "function load_") || continue
        occursin("scan", lowercase(s)) || continue
        push!(hits, "src/pnjl/PNJL.jl:$(line_no): $(s)")
    end
    return unique_sorted(hits)
end

function unique_sorted(lines::AbstractVector{<:AbstractString})
    return sort!(unique(String.(lines)))
end

function filter_non_definition(lines::AbstractVector{<:AbstractString})
    kept = String[]
    for line in lines
        occursin(r"function\s+run_(tmu|trho)_scan", line) && continue
        occursin(r"export\s+run_(tmu|trho)_scan", line) && continue
        push!(kept, line)
    end
    return unique_sorted(kept)
end

function filter_julia_callsites(lines::AbstractVector{<:AbstractString})
    kept = String[]
    for line in lines
        parts = split(String(line), ":"; limit=2)
        isempty(parts) && continue
        path = parts[1]
        endswith(path, ".jl") || continue
        push!(kept, String(line))
    end
    return unique_sorted(kept)
end

function result_timestamp()
    return Dates.format(now(), "yyyymmdd_HHMMSS")
end

function output_path(root::String, ts::String)
    out_dir = joinpath(root, "data", "outputs", "results")
    mkpath(out_dir)
    return joinpath(out_dir, "pnjl_decommission_checklist_$(ts).md")
end

function write_markdown(path::String, cfg::ChecklistConfig, scan_files, scan_includes, pnjl_export_callers, models_callers, pnjl_solver_callers, scan_default_includes, scan_loaders)
    generated = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    open(path, "w") do io
        println(io, "# PNJL 下线前检查表（scans/* + PNJL.jl）")
        println(io)
        println(io, "- generated_at: ", generated)
        println(io, "- diff_range: ", cfg.base, "...", cfg.head)
        println(io, "- scope: src/models/scans/*, src/pnjl/PNJL.jl")
        println(io)

        println(io, "## 1) 调用方盘点")
        println(io)
        println(io, "### 1.1 scans 文件清单")
        for item in scan_files
            println(io, "- [ ] ", item)
        end
        isempty(scan_files) && println(io, "- [ ] (未发现 scans 文件)")
        println(io)

        println(io, "### 1.2 对 scans 的 include/路径引用")
        for item in scan_includes
            println(io, "- [ ] ", item)
        end
        isempty(scan_includes) && println(io, "- [x] 未发现额外 include/路径引用")
        println(io)

        println(io, "### 1.2b PNJL 顶层默认 include 审计（排除白名单）")
        for item in scan_default_includes
            println(io, "- [ ] ", item)
        end
        isempty(scan_default_includes) && println(io, "- [x] 顶层默认 include 为 0（白名单除外）")
        println(io)

        println(io, "### 1.2c PNJL scan loader 可见性")
        for item in scan_loaders
            println(io, "- [ ] ", item)
        end
        isempty(scan_loaders) && println(io, "- [ ] 未发现 scan loader（需人工确认）")
        println(io)

        println(io, "### 1.3 `PNJL.run_tmu_scan/run_trho_scan` 调用方")
        for item in pnjl_export_callers
            println(io, "- [ ] ", item)
        end
        isempty(pnjl_export_callers) && println(io, "- [x] 代码库内未发现运行时调用")
        println(io)

        println(io, "### 1.4 `Models.run_tmu_scan/run_trho_scan` 调用方")
        for item in models_callers
            println(io, "- [ ] ", item)
        end
        isempty(models_callers) && println(io, "- [ ] 未发现（需人工确认是否存在外部调用）")
        println(io)

        println(io, "### 1.5 `PNJL.solve/solve_multi` 外部调用方（排除 src/pnjl/**）")
        for item in pnjl_solver_callers
            println(io, "- [ ] ", item)
        end
        isempty(pnjl_solver_callers) && println(io, "- [x] 未发现外部运行时调用")
        println(io)

        println(io, "## 2) 风险检查")
        println(io)
        println(io, "- [ ] 风险 R1：外部脚本或服务仍直接依赖 `PNJL.run_*`，下线后出现入口丢失。")
        println(io, "- [ ] 风险 R2：`src/models/scans/*` 迁移后 include 链或旧测试加载失败。")
        println(io, "- [ ] 风险 R3：扫描结果口径漂移（默认 backend/seed 差异）引发回归不一致。")
        println(io, "- [ ] 风险 R4：回退开关未覆盖扫描入口，故障时无法秒级止损。")
        println(io)

        println(io, "## 3) 回退策略")
        println(io)
        println(io, "- [ ] 回退 S1：保留 `PNJL.jl` 对 `Models.run_*` 的 compat 转发一个稳定周期。")
        println(io, "- [ ] 回退 S2：保留 `scripts/dev/run_prune_wave_gate.jl` 快照产物并对比删除前后调用差异。")
        println(io, "- [ ] 回退 S3：定向回归固定执行 `models_unified_entrypoints + pnjl scan backend smoke`。")
        println(io, "- [ ] 回退 S4：如出现生产故障，执行 `git revert <wave-commit>` 恢复 `src/models/scans/*` 与 `PNJL` 导出。")
    end
end

function main()
    cfg = parse_args(ARGS)
    root = project_root()

    scan_files = collect_scan_files(root)
    scan_includes = filter_julia_callsites(grep_lines(root, "scans/|src/models/scans|src/pnjl/scans", ["src", "tests", "scripts"]))
    scan_default_includes = collect_scan_default_includes(root)
    scan_loaders = collect_scan_loaders(root)

    pnjl_callers_raw = grep_lines(root, "PNJL\\.run_tmu_scan|PNJL\\.run_trho_scan", ["src", "tests", "scripts"])
    pnjl_export_callers = filter_non_definition(filter_julia_callsites(pnjl_callers_raw))

    models_callers_raw = grep_lines(root, "Models\\.run_tmu_scan|Models\\.run_trho_scan", ["src", "tests", "scripts"])
    models_callers = filter_non_definition(filter_julia_callsites(models_callers_raw))

    pnjl_solver_callers_raw = grep_lines(root, "PNJL\\.solve\\(|PNJL\\.solve_multi\\(", ["src", "tests", "scripts"])
    pnjl_solver_callers_all = filter_non_definition(filter_julia_callsites(pnjl_solver_callers_raw))
    pnjl_solver_callers = String[]
    for line in pnjl_solver_callers_all
        path = split(line, ':'; limit=2)[1]
        startswith(path, "src/pnjl/") && continue
        push!(pnjl_solver_callers, line)
    end
    pnjl_solver_callers = unique_sorted(pnjl_solver_callers)

    ts = result_timestamp()
    out = output_path(root, ts)
    write_markdown(out, cfg, scan_files, scan_includes, pnjl_export_callers, models_callers, pnjl_solver_callers, scan_default_includes, scan_loaders)

    println("[pnjl-decommission-checklist] output: ", normalize_path(relpath(out, root)))
    println("[pnjl-decommission-checklist] scans files: ", length(scan_files))
    println("[pnjl-decommission-checklist] default includes (excluding allowlist): ", length(scan_default_includes))
    println("[pnjl-decommission-checklist] scan loaders: ", length(scan_loaders))
    println("[pnjl-decommission-checklist] pnjl callers: ", length(pnjl_export_callers))
    println("[pnjl-decommission-checklist] models callers: ", length(models_callers))
    println("[pnjl-decommission-checklist] pnjl solver callers (external): ", length(pnjl_solver_callers))
end

main()
