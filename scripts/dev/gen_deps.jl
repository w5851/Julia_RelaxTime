#!/usr/bin/env julia
# Static include graph plus unresolved relative-import names; never evaluates source.

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
re_using_import = r"(?:using|import)\s+(\.*)([A-Za-z_][A-Za-z0-9_]*)"

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

function base_call_name(ex)
    ex isa Symbol && return ex
    if ex isa Expr && ex.head == :. && ex.args[1] == :Base && ex.args[2] isa QuoteNode
        return ex.args[2].value
    end
    return nothing
end

# Deliberately limited to path construction, not a Julia interpreter.
function static_path(ex, paths::Dict{Symbol,String}, file::String)
    ex isa String && return ex
    ex isa Symbol && return get(paths, ex, nothing)
    ex isa Expr || return nothing
    if ex.head == :macrocall
        ex.args[1] == Symbol("@__DIR__") && return dirname(file)
        ex.args[1] == Symbol("@__FILE__") && return file
    elseif ex.head == :call
        name = base_call_name(ex.args[1])
        name in (:joinpath, :normpath, :dirname, :abspath) || return nothing
        parts = [static_path(arg, paths, file) for arg in ex.args[2:end]]
        any(isnothing, parts) && return nothing
        isempty(parts) && return nothing
        name == :joinpath && return joinpath(parts...)
        name == :normpath && return normpath(joinpath(parts...))
        name == :dirname && length(parts) == 1 && return dirname(only(parts))
        # Relative abspath calls depend on the caller's working directory.
        name == :abspath && isabspath(first(parts)) && return abspath(joinpath(parts...))
    end
    return nothing
end

function forget_bindings!(paths, ex)
    if ex isa Symbol
        delete!(paths, ex)
    elseif ex isa Expr
        for arg in ex.args
            forget_bindings!(paths, arg)
        end
    end
end

function invalidate_assignments!(paths, ex)
    ex isa Expr || return
    ex.head in (:function, :module, :->, :quote, :inert) && return
    if ex.head == :(=) && ex.args[1] isa Symbol
        delete!(paths, ex.args[1])
    end
    for arg in ex.args
        invalidate_assignments!(paths, arg)
    end
end

function collect_includes!(targets, unresolved, ex, paths, file)
    ex isa Expr || return
    ex.head in (:quote, :inert) && return
    if ex.head in (:error, :incomplete)
        push!(unresolved, "parse error: $(first(ex.args))")
        return
    elseif ex.head in (:if, :elseif)
        for arg in ex.args
            collect_includes!(targets, unresolved, arg, copy(paths), file)
        end
        # A conditional assignment cannot establish a unique later path.
        invalidate_assignments!(paths, ex)
        return
    end
    if ex.head == :(=) && ex.args[1] isa Symbol
        value = static_path(ex.args[2], paths, file)
        if value === nothing
            delete!(paths, ex.args[1])
        else
            paths[ex.args[1]] = value
        end
    elseif ex.head == :call && base_call_name(ex.args[1]) == :include
        path = static_path(last(ex.args), paths, file)
        if path === nothing
            push!(unresolved, string(ex))
        else
            push!(targets, normpath(joinpath(dirname(file), path)))
        end
        return
    end
    # Modules have their own globals; function parameters can shadow file bindings.
    if ex.head == :module
        paths = Dict{Symbol,String}()
    elseif ex.head in (:function, :->) ||
           (ex.head == :(=) && ex.args[1] isa Expr)
        paths = copy(paths)
        forget_bindings!(paths, ex.args[1])
    elseif ex.head in (:let, :for, :while)
        paths = copy(paths)
    end
    for arg in ex.args
        collect_includes!(targets, unresolved, arg, paths, file)
    end
end

function collect_dependency_graph(root::String; unresolved::Vector{String}=String[])
    root = abspath(root)
    srcdir = joinpath(root, "src")
    files = collect_source_files(srcdir)
    adj = Dict{String, Set{String}}()

    for f in files
        source = read(f, String)
        text = sanitize_source_text(source)
        node = relpath_for(root, f)
        deps = get!(adj, node, Set{String}())
        targets, dynamic = String[], String[]
        collect_includes!(targets, dynamic, Meta.parseall(source; filename=f), Dict{Symbol,String}(), f)
        for path in targets
            target = relpath_for(root, path)
            target == node || push!(deps, target)
            isfile(path) || push!(unresolved, "$node: missing include target $target")
        end
        append!(unresolved, ["$node: $expr" for expr in dynamic])

        for m in eachmatch(re_using_import, text)
            leading = m.captures[1]
            modname = m.captures[2]
            if startswith(leading, ".")
                push!(deps, modname)
            end
        end
    end

    for k in collect(keys(adj))
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
    if length(parts) >= 3 && parts[1] == "src"
        return parts[2]
    elseif parts[1] == "scripts"
        return "scripts"
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
    unresolved = String[]
    adj = collect_dependency_graph(root; unresolved=unresolved)
    sccs = strongly_connected_components(adj)
    graph_text = render_mermaid(adj)

    mermaid = IOBuffer()
    println(mermaid, "# Dependency graph generated: ", Dates.now())
    println(mermaid, "\nRun: julia --project=. scripts/dev/gen_deps.jl\n")
    println(mermaid, "阅读入口：[职责与调用方向](dependencies.manual.md) · [依赖规则](dependency_rules.md)。\n")
    println(mermaid, "本图静态解析 include 的字面量、路径常量及 joinpath/normpath/dirname/@__DIR__；不执行源码。")
    println(mermaid, "相对 using/import 只显示模块名，尚未解析完整模块身份；条件分支合并展示。无环结果仅适用于已解析的边，不能证明完整运行时依赖无环。\n")
    println(mermaid, "## 未解析或缺失的 include\n")
    if isempty(unresolved)
        println(mermaid, "- 无。\n")
    else
        for item in sort(unique(unresolved))
            println(mermaid, "- `", item, "`")
        end
        println(mermaid)
    end
    println(mermaid, "## 静态依赖图\n")
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
