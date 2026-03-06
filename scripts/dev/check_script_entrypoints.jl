#!/usr/bin/env julia

path_separator() = Sys.iswindows() ? '\\' : '/'

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const TARGET_DIRS = [
    joinpath(ROOT, "benchmark", "pnjl"),
    joinpath(ROOT, "benchmark", "relaxtime"),
    joinpath(ROOT, "scripts", "perf", "pnjl"),
    joinpath(ROOT, "scripts", "perf", "relaxtime"),
    joinpath(ROOT, "scripts", "analysis"),
    joinpath(ROOT, "scripts", "pnjl"),
    joinpath(ROOT, "scripts", "relaxtime"),
]

const SKIP_SUBSTRINGS = [
    string(path_separator(), "_archived", path_separator()),
]

const DIRECT_INCLUDE_RE = r"\binclude\(\"([^\"]+)\"\)"
const JOINPATH_CALL_RE = r"joinpath\(\s*(@__DIR__|[A-Za-z_][A-Za-z0-9_]*)\s*,\s*(.+)\)"
const ROOT_ASSIGN_RE = r"^(?:const\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*=\s*normpath\(joinpath\(@__DIR__,\s*(.+)\)\)"

should_skip(path::String) = any(needle -> occursin(needle, path), SKIP_SUBSTRINGS)

parse_literal_components(expr::AbstractString) = [m.captures[1] for m in eachmatch(r"\"([^\"]+)\"", expr)]

function collect_julia_files(dir::String)
    if !isdir(dir)
        return String[]
    end

    files = String[]
    for (root, _, names) in walkdir(dir)
        should_skip(root) && continue
        for name in names
            endswith(name, ".jl") || continue
            push!(files, joinpath(root, name))
        end
    end
    sort!(files)
    return files
end

function check_ref!(violations::Vector{NamedTuple}, file::String, line_no::Int, rel::AbstractString, kind::String)
    isempty(rel) && return
    isabspath(rel) && return
    target = normpath(joinpath(dirname(file), rel))
    ispath(target) && return
    push!(violations, (file=file, line=line_no, ref=String(rel), resolved=target, kind=kind))
end

function collect_root_aliases(lines::Vector{String}, file::String)
    aliases = Dict{String, String}()
    for line in lines
        stripped = strip(line)
        startswith(stripped, '#') && continue
        m = match(ROOT_ASSIGN_RE, stripped)
        m === nothing && continue
        components = parse_literal_components(m.captures[2])
        isempty(components) && continue
        aliases[m.captures[1]] = normpath(joinpath(dirname(file), components...))
    end
    return aliases
end

function check_joinpath_ref!(violations::Vector{NamedTuple}, file::String, line_no::Int, base::AbstractString, tail_expr::AbstractString, kind::String, aliases::Dict{String, String})
    parts = parse_literal_components(tail_expr)
    isempty(parts) && return

    base_path = if base == "@__DIR__"
        dirname(file)
    else
        get(aliases, base, "")
    end
    isempty(base_path) && return

    target = normpath(joinpath(base_path, parts...))
    ispath(target) && return
    push!(violations, (file=file, line=line_no, ref="joinpath($base, $(join(parts, ", ")))", resolved=target, kind=kind))
end

function scan_file(file::String, violations::Vector{NamedTuple})
    lines = readlines(file)
    aliases = collect_root_aliases(lines, file)

    for (line_no, line) in enumerate(lines)
        stripped = strip(line)
        startswith(stripped, '#') && continue

        m = match(DIRECT_INCLUDE_RE, line)
        if m !== nothing
            check_ref!(violations, file, line_no, m.captures[1], "include")
        end

        if occursin("joinpath(", line)
            join_match = match(JOINPATH_CALL_RE, line)
            if join_match !== nothing
                base = join_match.captures[1]
                tail_expr = join_match.captures[2]
                if occursin("include(", line)
                    check_joinpath_ref!(violations, file, line_no, base, tail_expr, "joinpath-include", aliases)
                elseif occursin("push!(LOAD_PATH", line)
                    check_joinpath_ref!(violations, file, line_no, base, tail_expr, "load-path", aliases)
                elseif occursin("Pkg.activate(", line)
                    check_joinpath_ref!(violations, file, line_no, base, tail_expr, "pkg-activate", aliases)
                end
            end
        end
    end
end

function main()
    files = reduce(vcat, (collect_julia_files(dir) for dir in TARGET_DIRS); init=String[])
    violations = NamedTuple[]

    for file in files
        scan_file(file, violations)
    end

    if !isempty(violations)
        println("[script-entrypoints] FAILED: $(length(violations)) broken relative references")
        for item in violations
            rel_file = relpath(item.file, ROOT)
            println(" - $(rel_file):$(item.line) [$(item.kind)] $(item.ref) -> $(item.resolved)")
        end
        exit(1)
    end

    println("[script-entrypoints] OK: checked $(length(files)) Julia files")
    for dir in TARGET_DIRS
        if isdir(dir)
            println("  - " * relpath(dir, ROOT))
        end
    end
end

main()