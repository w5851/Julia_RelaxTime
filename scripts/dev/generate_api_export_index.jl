#!/usr/bin/env julia

"""
生成导出 API 全集索引（纯脚本版）。

用途：
1. 扫描目标 Julia 文件中的 export 列表；
2. 生成一份 Markdown 索引页；
3. 统计这些导出符号是否已在 docs/api 下的人工文档中被提及。

用法示例：

  julia --project=. scripts/dev/generate_api_export_index.jl \
    --module-file src/models/Models.jl \
        --module-file src/models/entrypoints.jl \
        --include-symbol run_phase_pipeline \
        --include-symbol find_cep \
        --output docs/api/models/phase/generated/Exports.md \
    --title "Phase Diagram Export API Index"
"""

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function relpath_display(path::String, root::String)
    return replace(relpath(path, root), '\\' => '/')
end

function parse_args(args::Vector{String})
    options = Dict{String, Any}(
        "module-file" => String[],
        "include-symbol" => String[],
        "output" => "",
        "title" => "Export API Index",
        "doc-root" => "docs/api",
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unknown positional argument: $arg")
        key = arg[3:end]
        haskey(options, key) || error("Unknown option: --$key")
        i == length(args) && error("Missing value for option: --$key")
        if key in ("module-file", "include-symbol")
            push!(options[key], args[i + 1])
        else
            options[key] = args[i + 1]
        end
        i += 2
    end

    isempty(options["module-file"]) && error("at least one --module-file is required")
    isempty(options["output"]) && error("--output is required")
    return options
end

function strip_comment(line::String)
    idx = findfirst(==('#'), line)
    idx === nothing && return line
    return line[1:prevind(line, idx)]
end

function parse_export_symbols(module_files::Vector{String})
    exports = Dict{String, NamedTuple}()
    seen = Set{String}()

    for module_file in module_files
        lines = readlines(module_file)
        in_export = false
        buffer = ""
        start_line = 0

        for (line_no, raw_line) in enumerate(lines)
            line = strip(strip_comment(raw_line))
            isempty(line) && !in_export && continue

            if !in_export && startswith(line, "export ")
                in_export = true
                start_line = line_no
                buffer = strip(line[8:end])
            elseif in_export
                buffer *= " " * line
            else
                continue
            end

            if endswith(strip(line), ",")
                continue
            end

            for item in split(buffer, ',')
                symbol = strip(item)
                isempty(symbol) && continue
                if symbol in seen
                    entry = exports[symbol]
                    files = copy(entry.files)
                    lines_out = copy(entry.lines)
                    module_file in files || push!(files, module_file)
                    start_line in lines_out || push!(lines_out, start_line)
                    exports[symbol] = (symbol=symbol, files=files, lines=lines_out)
                    continue
                end
                exports[symbol] = (symbol=symbol, files=[module_file], lines=[start_line])
                push!(seen, symbol)
            end

            in_export = false
            buffer = ""
            start_line = 0
        end
    end

    values_out = collect(values(exports))
    sort!(values_out; by=entry -> entry.symbol)
    return values_out
end

function collect_markdown_files(dir::String)
    files = String[]
    isdir(dir) || return files
    for (root, _, names) in walkdir(dir)
        for name in names
            endswith(lowercase(name), ".md") || continue
            push!(files, joinpath(root, name))
        end
    end
    sort!(files)
    return files
end

function should_skip_doc(path::String, output_path::String)
    norm = normpath(path)
    out = normpath(output_path)
    norm == out && return true
    parts = split(replace(norm, '\\' => '/'), '/')
    return "generated" in parts
end

function build_doc_mentions(symbols, doc_files::Vector{String}, output_path::String, root::String)
    mentions = Dict{String, Vector{String}}()
    for entry in symbols
        mentions[entry.symbol] = String[]
    end

    for file in doc_files
        should_skip_doc(file, output_path) && continue
        content = read(file, String)
        rel = relpath_display(file, root)
        for entry in symbols
            occursin(entry.symbol, content) || continue
            push!(mentions[entry.symbol], rel)
        end
    end

    return mentions
end

function filter_symbols(symbols, include_symbols::Vector{String})
    isempty(include_symbols) && return symbols
    include_set = Set(include_symbols)
    filtered = [entry for entry in symbols if entry.symbol in include_set]
    isempty(filtered) && error("no exported symbols matched the provided --include-symbol filters")
    return filtered
end

function render_markdown(title::String, module_files::Vector{String}, output_file::String, symbols, mentions, root::String, include_symbols::Vector{String})
    rel_modules = [relpath_display(path, root) for path in module_files]
    rel_output = relpath_display(output_file, root)
    documented = count(entry -> !isempty(mentions[entry.symbol]), symbols)
    undocumented = [entry.symbol for entry in symbols if isempty(mentions[entry.symbol])]

    lines = String[]
    push!(lines, "# " * title)
    push!(lines, "")
    push!(lines, "> Auto-generated by `scripts/dev/generate_api_export_index.jl`. Do not maintain this page by hand.")
    push!(lines, "")
    push!(lines, "- Source module files: " * join(["`$(path)`" for path in rel_modules], ", "))
    push!(lines, "- Generated page: `$(rel_output)`")
    push!(lines, "- Exported symbol count: $(length(symbols))")
    push!(lines, "- Mentioned in non-generated docs/api pages: $(documented)")
    if !isempty(include_symbols)
        push!(lines, "- Topic filter: " * join(["`$(symbol)`" for symbol in include_symbols], ", "))
    end
    push!(lines, "")
    push!(lines, "## Exported Symbols")
    push!(lines, "")
    push!(lines, "| Symbol | Source files | Export lines | Mentioned in docs | Example docs |")
    push!(lines, "| --- | --- | ---: | ---: | --- |")

    for entry in symbols
        docs = mentions[entry.symbol]
        examples = isempty(docs) ? "-" : join(first(docs, min(length(docs), 3)), "<br>")
        source_files = join([relpath_display(path, root) for path in entry.files], "<br>")
        export_lines = join([string(line) for line in entry.lines], ", ")
        push!(lines, "| `$(entry.symbol)` | $(source_files) | $(export_lines) | $(length(docs)) | $(examples) |")
    end

    push!(lines, "")
    push!(lines, "## Undocumented Or Not Yet Mentioned")
    push!(lines, "")
    if isempty(undocumented)
        push!(lines, "All exported symbols were found in at least one non-generated docs/api Markdown page.")
    else
        push!(lines, "The following exported symbols were not found in non-generated docs/api Markdown pages:")
        push!(lines, "")
        for symbol in undocumented
            push!(lines, "- `$(symbol)`")
        end
    end

    push!(lines, "")
    push!(lines, "## Notes")
    push!(lines, "")
    push!(lines, "- This page treats `export` statements as the authoritative public-surface baseline.")
    push!(lines, "- Mention counts are string-match based and serve as a lightweight coverage signal, not a semantic guarantee.")

    return join(lines, "\n") * "\n"
end

function main()
    options = parse_args(ARGS)
    root = project_root()
    module_files = [normpath(joinpath(root, path)) for path in options["module-file"]]
    output_file = normpath(joinpath(root, options["output"]))
    doc_root = normpath(joinpath(root, options["doc-root"]))
    title = options["title"]
    include_symbols = options["include-symbol"]

    for module_file in module_files
        isfile(module_file) || error("module file not found: $module_file")
    end

    symbols = parse_export_symbols(module_files)
    symbols = filter_symbols(symbols, include_symbols)
    isempty(symbols) && error("no exports found in selected module files")

    doc_files = collect_markdown_files(doc_root)
    mentions = build_doc_mentions(symbols, doc_files, output_file, root)
    content = render_markdown(title, module_files, output_file, symbols, mentions, root, include_symbols)

    mkpath(dirname(output_file))
    write(output_file, content)

    println("[generate-api-export-index] OK")
    println(" - module files: " * join([relpath_display(path, root) for path in module_files], ", "))
    println(" - output: " * relpath_display(output_file, root))
    println(" - exports: " * string(length(symbols)))
end

main()