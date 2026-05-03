#!/usr/bin/env julia
# Minimal dependency graph generator for Julia project
# Scans `src/` for .jl files, extracts include("..."), using .Module and import .Module
# Produces a Mermaid `graph TD` output into `docs/architecture/dependencies.md`.

using Printf
using Dates

function find_mmdc()
    root = pwd()
    if Sys.which("mmdc") !== nothing
        return Sys.which("mmdc")
    end
    local_bin = joinpath(root, "node_modules", ".bin", Sys.iswindows() ? "mmdc.cmd" : "mmdc")
    return isfile(local_bin) ? local_bin : nothing
end

# regexes
re_include = r"include\([\"']([^\"']+)[\"']\)"
re_using_import = r"(?:using|import)\s+(\.*)([A-Za-z_][A-Za-z0-9_]*)"

function add_edge(a,b)
    if a == b
        return
    end
    adj_a = get!(ACTIVE_ADJ[], a, Set{String}())
    push!(adj_a, b)
end

# adjacency map (string => Set{String}) using node names as repo-relative paths or module names
const ACTIVE_ADJ = Ref(Dict{String, Set{String}}())

# helper: node id sanitize
function id_for(node)
    # replace non-alnum with underscore
    return replace(node, r"[^A-Za-z0-9_]" => "_")
end

function strip_triple_quoted_blocks(text::String)
    io = IOBuffer()
    i = firstindex(text)
    in_triple = false
    while i <= lastindex(text)
        if startswith(SubString(text, i), "\"\"\"")
            print(io, "   ")
            i = nextind(text, i, 3)
            in_triple = !in_triple
            continue
        end

        c = text[i]
        if in_triple
            print(io, c == '\n' ? '\n' : ' ')
        else
            print(io, c)
        end
        i = nextind(text, i)
    end
    return String(take!(io))
end

function strip_line_comments(text::String)
    lines = split(text, '\n'; keepempty=true)
    stripped = String[]
    sizehint!(stripped, length(lines))

    for line in lines
        io = IOBuffer()
        in_string = false
        escaped = false
        for c in line
            if c == '"' && !escaped
                in_string = !in_string
            elseif c == '#' && !in_string
                break
            end
            print(io, c)
            escaped = (c == '\\' && !escaped)
            if c != '\\'
                escaped = false
            end
        end
        push!(stripped, String(take!(io)))
    end

    return join(stripped, '\n')
end

sanitize_source_text(text::String) = strip_line_comments(strip_triple_quoted_blocks(text))

function collect_source_files(srcdir::String)
    files = String[]
    for (d, _ds, fs) in walkdir(srcdir)
        for f in fs
            if endswith(f, ".jl")
                push!(files, joinpath(d, f))
            end
        end
    end
    return files
end

function relpath_for(root::String, path::String)
    return replace(relpath(path, root), "\\" => "/")
end

function collect_dependency_graph(root::String)
    srcdir = joinpath(root, "src")
    files = collect_source_files(srcdir)
    adj = Dict{String, Set{String}}()
    ACTIVE_ADJ[] = adj

    for f in files
        text = sanitize_source_text(read(f, String))
        node = relpath_for(root, f)

        for m in eachmatch(re_include, text)
            inc = m.captures[1]
            inc_path = normpath(joinpath(dirname(f), inc))
            if isfile(inc_path)
                target = relpath_for(root, inc_path)
                add_edge(node, target)
            else
                add_edge(node, inc)
            end
        end

        for m in eachmatch(re_using_import, text)
            leading = m.captures[1]
            modname = m.captures[2]
            if startswith(leading, ".")
                add_edge(node, modname)
            end
        end
    end

    for k in keys(adj)
        get!(adj, k, Set{String}())
        for v in adj[k]
            get!(adj, v, Set{String}())
        end
    end

    return adj
end

# Tarjan's algorithm for SCC detection
function strongly_connected_components(adj::Dict{String, Set{String}})
    index = Dict{String, Int}()
    lowlink = Dict{String, Int}()
    stack = String[]
    onstack = Set{String}()
    idx = Ref(0)
    sccs = Vector{Vector{String}}()

    function strongconnect(v)
        get!(adj, v, Set{String}())
        idx[] += 1
        index[v] = idx[]
        lowlink[v] = idx[]
        push!(stack, v)
        push!(onstack, v)
        for w in get(adj, v, Set{String}())
            if !haskey(index, w)
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elseif w in onstack
                lowlink[v] = min(lowlink[v], index[w])
            end
        end
        if lowlink[v] == index[v]
            comp = String[]
            while true
                w = pop!(stack)
                delete!(onstack, w)
                push!(comp, w)
                if w == v
                    break
                end
            end
            push!(sccs, comp)
        end
    end

    for v in keys(adj)
        if !haskey(index, v)
            strongconnect(v)
        end
    end

    return sccs
end

# build mermaid graph
# group nodes by top-level folder under src if available
function top_level_group(node)
    parts = split(node, '/');
    if length(parts) >= 2 && parts[1] == "src"
        return parts[2]
    else
        return "root"
    end
end

function render_mermaid(adj::Dict{String, Set{String}})
    groups = Dict{String, Vector{String}}()
    for n in keys(adj)
        g = top_level_group(n)
        push!(get!(groups, g, String[]), n)
    end

    idmap = Dict{String,String}()
    for n in keys(adj)
        idmap[n] = id_for(n)
    end

    graphbuf = IOBuffer()
    direction = get(ENV, "MERMAID_DIRECTION", "LR")
    node_spacing = get(ENV, "MERMAID_NODE_SPACING", "40")
    rank_spacing = get(ENV, "MERMAID_RANK_SPACING", "60")
    println(graphbuf, "%%{init: { 'flowchart': { 'nodeSpacing': ", node_spacing, ", 'rankSpacing': ", rank_spacing, ", 'useMaxWidth': false } }}%%")
    println(graphbuf, "flowchart ", direction)

    for (g, nodes) in sort(collect(groups); by=x->x[1])
        println(graphbuf, "  subgraph ", g)
        for n in sort(nodes)
            nid = idmap[n]
            label = replace(n, "src/" => "")
            println(graphbuf, "    ", nid, "[", label, "]")
        end
        println(graphbuf, "  end")
    end

    for a in sort(collect(keys(adj)))
        for b in sort(collect(adj[a]))
            ida = idmap[a]
            idb = haskey(idmap, b) ? idmap[b] : id_for(b)
            println(graphbuf, "  ", ida, " --> ", idb)
        end
    end

    return String(take!(graphbuf))
end

function main()
    root = pwd()
    outdir = joinpath(root, "docs", "architecture")
    mkpath(outdir)
    outfile = joinpath(outdir, "dependencies.md")
    mmdfile = joinpath(outdir, "dependencies.mmd")
    svgfile = joinpath(outdir, "dependencies.svg")
    manualfile = joinpath(outdir, "dependencies.manual.md")

    adj = collect_dependency_graph(root)
    sccs = strongly_connected_components(adj)
    graph_text = render_mermaid(adj)
    manual_content = isfile(manualfile) ? read(manualfile, String) : "## L1/L3 (manual)\n\n请在 docs/architecture/dependencies.manual.md 中补充 L1/L3 内容。\n"

    mermaid = IOBuffer()
    println(mermaid, "# Dependency graph generated: ", Dates.now())
    println(mermaid, "\nRun: julia --project=. scripts/dev/gen_deps.jl\n")
    print(mermaid, manual_content)
    println(mermaid, "\n---\n")
    println(mermaid, "![Dependency graph](dependencies.svg)")
    println(mermaid, "\n---\n")
    println(mermaid, "```mermaid")
    print(mermaid, graph_text)
    println(mermaid, "```\n")

    cycles = [c for c in sccs if length(c) > 1]
    if !isempty(cycles)
        println(mermaid, "### Detected cycles (strongly connected components)")
        for c in cycles
            println(mermaid, "- ", join(c, " -> "))
        end
        println(mermaid, "\n")
    end

    open(mmdfile, "w") do io
        write(io, graph_text)
    end

    open(outfile, "w") do io
        write(io, String(take!(mermaid)))
    end

    println("Wrote dependency graph to: ", outfile)
    println("Wrote Mermaid source to: ", mmdfile)

    svg_width = get(ENV, "MERMAID_WIDTH", "2200")
    svg_height = get(ENV, "MERMAID_HEIGHT", "1400")
    svg_scale = get(ENV, "MERMAID_SCALE", "1.5")
    mmdc_path = find_mmdc()
    if mmdc_path !== nothing
        try
            cmd = Cmd([mmdc_path, "-i", mmdfile, "-o", svgfile, "-w", svg_width, "-H", svg_height, "-s", svg_scale])
            run(cmd)
            println("Wrote SVG to: ", svgfile)
        catch err
            @warn "mmdc failed" err
        end
    else
        println("mmdc not found; skip SVG render. Install with: npm install -g @mermaid-js/mermaid-cli")
        println("or use local dev dependency: npm install -D @mermaid-js/mermaid-cli")
    end

    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
