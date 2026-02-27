#!/usr/bin/env julia

using Dates

normalize_path(path::AbstractString) = replace(path, '\\' => '/')

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function latest_file(pattern::Regex, dir::String)
    isdir(dir) || return nothing
    files = filter(f -> occursin(pattern, basename(f)), readdir(dir; join=true))
    isempty(files) && return nothing
    sort!(files)
    return files[end]
end

function read_text(path::AbstractString)
    isfile(path) || return ""
    return read(path, String)
end

function has_line(path::String, rx::Regex)
    isfile(path) || return false
    for line in eachline(path)
        occursin(rx, line) && return true
    end
    return false
end

function scan_solver_direct_refs(root::String)
    scans_dir = joinpath(root, "src", "models", "scans")
    isdir(scans_dir) || return String[]
    hits = String[]
    for name in sort(readdir(scans_dir))
        endswith(name, ".jl") || continue
        path = joinpath(scans_dir, name)
        line_no = 0
        for line in eachline(path)
            line_no += 1
            s = strip(line)
            if occursin("ImplicitSolver", s) || occursin("using ..ImplicitSolver", s)
                push!(hits, "src/models/scans/$(name):$(line_no): $(s)")
            end
        end
    end
    return hits
end

function solver_core_impl_markers(root::String)
    path = joinpath(root, "src", "pnjl", "solver", "ImplicitSolver.jl")
    markers = Dict{String,Bool}(
        "function solve_with_derivatives(" => false,
        "function create_implicit_solver(" => false,
        "thermo_backend::Symbol=:models [FixedMu default]" => false,
        "thermo_backend::Symbol=:models [Derivative default]" => false,
    )
    isfile(path) || return markers
    txt = read(path, String)

    markers["function solve_with_derivatives("] = occursin("function solve_with_derivatives(", txt)
    markers["function create_implicit_solver("] = occursin("function create_implicit_solver(", txt)
    markers["thermo_backend::Symbol=:models [FixedMu default]"] = occursin(r"function solve\(::FixedMu[\s\S]*?thermo_backend::Symbol=:models", txt)
    markers["thermo_backend::Symbol=:models [Derivative default]"] = occursin(r"create_implicit_solver\s*=\s*function\s*\([\s\S]*?thermo_backend::Symbol=:models", txt) && occursin(r"solve_with_derivatives\s*=\s*function\s*\([\s\S]*?thermo_backend::Symbol=:models", txt)
    return markers
end

function bool_mark(b::Bool)
    return b ? "[x]" : "[ ]"
end

function main()
    root = project_root()
    out_dir = joinpath(root, "data", "outputs", "results")
    mkpath(out_dir)

    v2_path = joinpath(root, "docs", "dev", "active", "2026-02-25_新架构主线迁移开发任务单_v2.md")

    checklist = latest_file(r"^pnjl_decommission_checklist_\d{8}_\d{6}\.md$", out_dir)
    prune = latest_file(r"^pnjl_prune_wave_snapshot_\d{8}_\d{6}\.txt$", out_dir)
    allow = latest_file(r"^models_invokelatest_allowlist_delta_\d{8}_\d{6}\.txt$", out_dir)
    include_audit = latest_file(r"^pnjl_scan_default_include_audit_\d{8}_\d{6}\.txt$", out_dir)

    checklist_txt = checklist === nothing ? "" : read_text(checklist)
    prune_txt = prune === nothing ? "" : read_text(prune)

    scans_solver_refs = scan_solver_direct_refs(root)
    solver_markers = solver_core_impl_markers(root)

        dod1 = isempty(scans_solver_refs) &&
            !solver_markers["function solve_with_derivatives("] &&
            !solver_markers["function create_implicit_solver("] &&
            solver_markers["thermo_backend::Symbol=:models [FixedMu default]"] &&
            solver_markers["thermo_backend::Symbol=:models [Derivative default]"]
    dod2 = has_line(v2_path, r"^- \[x\] 标注 `src/pnjl/core/Thermodynamics\.jl`") &&
           has_line(v2_path, r"^- \[x\] 完成替代路径接线并去除直连调用点") &&
           has_line(v2_path, r"^- \[x\] 保留 compat API 但移除核心计算承载")
            dod3 = (occursin("pnjl callers: 0", checklist_txt) || occursin("代码库内未发现运行时调用", checklist_txt)) &&
                    (occursin("pnjl solver callers (external): 0", checklist_txt) || occursin("未发现外部运行时调用", checklist_txt))

    ts = Dates.format(now(), "yyyymmdd_HHMMSS")
    out = joinpath(out_dir, "v2_dod_closure_report_$(ts).md")

    open(out, "w") do io
        println(io, "# V2 DoD 收口审计报告")
        println(io)
        println(io, "- generated_at: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "- workspace: ", normalize_path(root))
        println(io)

        println(io, "## DoD 状态")
        println(io)
        println(io, "- ", bool_mark(dod1), " DoD-1 核心求解默认由 models 承载且 compat 不再承载核心实现")
        println(io, "- ", bool_mark(dod2), " DoD-2 Thermodynamics 可替代核心函数迁移收口")
        println(io, "- ", bool_mark(dod3), " DoD-3 scans/workflows 达到实体下线前置条件")
        println(io)

        println(io, "## 证据文件")
        println(io)
        println(io, "- checklist: ", checklist === nothing ? "(missing)" : normalize_path(relpath(checklist, root)))
        println(io, "- prune snapshot: ", prune === nothing ? "(missing)" : normalize_path(relpath(prune, root)))
        println(io, "- allowlist delta: ", allow === nothing ? "(missing)" : normalize_path(relpath(allow, root)))
        println(io, "- include audit: ", include_audit === nothing ? "(missing)" : normalize_path(relpath(include_audit, root)))
        println(io)

        println(io, "## 诊断（DoD-1 未闭环时）")
        println(io)
        println(io, "- ", bool_mark(solver_markers["function solve_with_derivatives("]), " compat core marker: `function solve_with_derivatives(`")
        println(io, "- ", bool_mark(solver_markers["function create_implicit_solver("]), " compat core marker: `function create_implicit_solver(`")
        println(io, "- ", bool_mark(solver_markers["thermo_backend::Symbol=:models [FixedMu default]"]), " models-default marker: `solve(::FixedMu)`")
        println(io, "- ", bool_mark(solver_markers["thermo_backend::Symbol=:models [Derivative default]"]), " models-default marker: implicit derivative wrappers")
        if isempty(scans_solver_refs)
            println(io, "- [x] src/models/scans/*.jl 对 ImplicitSolver 直接引用=0")
        else
            println(io, "- [ ] src/models/scans/*.jl 仍有 ImplicitSolver 直接引用：")
            for hit in scans_solver_refs
                println(io, "  - ", hit)
            end
        end
    end

    println("[v2-dod-closure] output: ", normalize_path(relpath(out, root)))
    println("[v2-dod-closure] DoD-1: ", dod1 ? "ready" : "not-ready")
    println("[v2-dod-closure] DoD-2: ", dod2 ? "ready" : "not-ready")
    println("[v2-dod-closure] DoD-3: ", dod3 ? "ready" : "not-ready")
end

main()
