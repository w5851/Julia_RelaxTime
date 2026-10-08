#!/usr/bin/env julia

module ArchiveDocs

using Dates
using JSON3

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const STATUSES = ("completed", "cancelled", "superseded")
const REQUIRED_FIELDS = ("title", "archived", "original", "archived_date")

path_key(path) = Sys.iswindows() ? lowercase(normpath(path)) : normpath(path)
same_path(a, b) = path_key(a) == path_key(b)
exists(path) = ispath(path) || islink(path)

function valid_date(value::AbstractString)
    occursin(r"^\d{4}-\d{2}-\d{2}$", value) || return false
    try
        return Dates.format(Date(value, dateformat"yyyy-mm-dd"), "yyyy-mm-dd") == value
    catch
        return false
    end
end

function document_dir(root, kind)
    base = realpath(root)
    for part in ("docs", "dev", kind)
        base = joinpath(base, part)
        if exists(base)
            isdir(base) && same_path(realpath(base), base) || throw(ArgumentError("document directory must not redirect through a link: $base"))
        end
    end
    return base
end

function resolve_document(root, kind, input)
    directory = document_dir(root, kind)
    path = if isabspath(input)
        normpath(input)
    elseif basename(input) == input
        joinpath(directory, input)
    else
        normpath(joinpath(realpath(root), input))
    end
    same_path(dirname(path), directory) || throw(ArgumentError("file must be directly under docs/dev/$kind: $input"))
    endswith(lowercase(path), ".md") || throw(ArgumentError("expected a Markdown file: $input"))
    isfile(path) && same_path(realpath(path), path) || throw(ArgumentError("file missing or redirected through a link: $input"))
    return path
end

function extract_title(content, filename)
    for line in split(content, '\n')
        heading = match(r"^#\s+(.+)$", strip(line))
        heading === nothing || return String(strip(heading.captures[1]))
    end
    return replace(replace(filename, r"^\d{4}[-_]\d{2}[-_]\d{2}[-_]" => ""), r"\.md$" => "")
end

# This format is a flat scalar mapping, not a general YAML document.
# JSON strings are valid YAML scalars and preserve punctuation and Unicode.
function read_scalar(value)
    value = strip(value)
    isempty(value) && throw(ArgumentError("empty metadata scalar"))
    if startswith(value, '"')
        parsed = JSON3.read(value)
        parsed isa AbstractString || throw(ArgumentError("expected a quoted string"))
        return String(parsed)
    elseif startswith(value, '\'')
        occursin(r"^'(?:[^']|'')*'$", value) || throw(ArgumentError("invalid single-quoted scalar"))
        return replace(chop(value; head=1, tail=1), "''" => "'")
    elseif lowercase(value) == "true"
        return true
    elseif lowercase(value) == "false"
        return false
    elseif lowercase(value) in ("null", "~")
        return nothing
    end
    number = tryparse(Float64, value)
    number === nothing || return number
    (occursin(r"[:#](?:\s|$)|\s#", value) || occursin(r"^[\[\]{}&*!|>@\x60%?]", value)) &&
        throw(ArgumentError("unsupported plain scalar; quote this metadata value"))
    return String(value)
end

function metadata_fields(content)
    lines = split(replace(content, "\r\n" => "\n"), '\n')
    !isempty(lines) && first(lines) == "---" || throw(ArgumentError("missing frontmatter"))
    closing = findnext(==("---"), lines, 2)
    closing === nothing && throw(ArgumentError("unclosed frontmatter"))
    fields = Dict{String,Any}()
    for line in lines[2:closing-1]
        isempty(strip(line)) && continue
        startswith(strip(line), '#') && continue
        entry = match(r"^([A-Za-z_][A-Za-z0-9_-]*):\s*(.+)$", line)
        entry === nothing && throw(ArgumentError("expected flat scalar frontmatter"))
        key = String(entry.captures[1])
        haskey(fields, key) && throw(ArgumentError("duplicate metadata field: $key"))
        fields[key] = read_scalar(entry.captures[2])
    end
    all(key -> haskey(fields, key), REQUIRED_FIELDS) || throw(ArgumentError("missing required metadata fields"))
    fields["archived"] === true || throw(ArgumentError("archived must be true"))
    for key in ("title", "original", "archived_date")
        fields[key] isa AbstractString && !isempty(strip(fields[key])) || throw(ArgumentError("$key must be a nonempty string"))
    end
    valid_date(fields["archived_date"]) || throw(ArgumentError("invalid archived_date"))
    occursin(r"^docs/dev/active/[^/\\]+(?i:\.md)$", fields["original"]) || throw(ArgumentError("invalid original path"))
    if haskey(fields, "task_status")
        fields["task_status"] in STATUSES || throw(ArgumentError("invalid task_status"))
        if fields["task_status"] != "completed"
            reason = get(fields, "archive_reason", "")
            reason isa AbstractString && !isempty(strip(reason)) || throw(ArgumentError("cancelled/superseded tasks require archive_reason"))
        end
    end
    return fields
end

function check_archived_format(path)
    try
        metadata_fields(read(path, String))
        return true
    catch
        return false
    end
end

function prepare_archive(input; root=ROOT, date=Dates.format(today(), "yyyy-mm-dd"), status="completed", reason="")
    valid_date(date) || throw(ArgumentError("--date must be a valid YYYY-MM-DD date"))
    status in STATUSES || throw(ArgumentError("--status must be completed, cancelled or superseded"))
    status == "completed" || !isempty(strip(reason)) || throw(ArgumentError("--reason is required for $status"))
    source = resolve_document(root, "active", input)
    filename = basename(source)
    match_date = match(r"^(\d{4})[-_](\d{2})[-_](\d{2})[-_](.+)$", filename)
    target_name = match_date === nothing ? "$(date)_$filename" : join(match_date.captures[1:3], "-") * "_" * match_date.captures[4]
    target = joinpath(document_dir(root, "archived"), target_name)
    exists(target) && throw(ArgumentError("archive destination already exists: $target"))
    original = read(source)
    title = extract_title(String(copy(original)), filename)
    header = "---\ntitle: $(JSON3.write(title))\narchived: true\noriginal: $(JSON3.write("docs/dev/active/$filename"))\narchived_date: $(JSON3.write(date))\ntask_status: $(JSON3.write(status))\n"
    isempty(reason) || (header *= "archive_reason: $(JSON3.write(reason))\n")
    header *= "---\n\n以下为原始内容（保留，以便审阅与历史参考）：\n\n---\n\n"
    archived = vcat(Vector{UInt8}(codeunits(header)), original)
    metadata_fields(String(copy(archived)))
    return (; root=realpath(root), source, target, original, archived)
end

function execute_archive(plan)
    resolve_document(plan.root, "active", plan.source)
    same_path(dirname(plan.target), document_dir(plan.root, "archived")) ||
        throw(ArgumentError("archive directory changed after preflight"))
    exists(plan.target) && throw(ArgumentError("archive destination already exists: $(plan.target)"))
    read(plan.source) == plan.original || throw(ArgumentError("source changed after archive preflight"))
    mkpath(dirname(plan.target))
    temporary, stream = mktemp(dirname(plan.target))
    try
        write(stream, plan.archived)
        close(stream)
        read(temporary) == plan.archived || error("archive write verification failed")
        # Publish a complete file without replacing a concurrently created target.
        # Unsupported filesystems fail closed, retaining the source.
        hardlink(temporary, plan.target)
        read(plan.target) == plan.archived || error("published archive verification failed; source retained")
        resolve_document(plan.root, "active", plan.source)
        read(plan.source) == plan.original || error("source changed during archive; both files retained")
        rm(plan.source)
    finally
        isopen(stream) && close(stream)
        isfile(temporary) && rm(temporary)
    end
    return plan.target
end

function print_usage(io=stdout)
    println(io, "Usage: julia --project=. scripts/dev/archive_docs.jl [OPTIONS] [FILES...]")
    println(io, "  --dry-run             Validate and show the exact move; write nothing")
    println(io, "  -c, --check [FILES...] Check selected archived files, or all when omitted")
    println(io, "  -d, --date DATE       Valid YYYY-MM-DD date (default: today)")
    println(io, "  --status STATUS       completed (default), cancelled, superseded")
    println(io, "  --reason TEXT         Required for cancelled/superseded tasks")
    println(io, "  -b, --batch           Archive explicitly listed files")
    println(io, "  -i, --interactive     Select active files interactively")
    println(io, "  -h, --help            Show help")
    println(io, "Only direct children of active/ may be archived. Existing targets are never overwritten.")
end

function main(args=collect(String.(ARGS)); root=ROOT, io=stdout, input=stdin)
    try
        isempty(args) && (print_usage(io); return 0)
        date = Dates.format(today(), "yyyy-mm-dd")
        status, reason = "completed", ""
        dry_run, check, interactive = false, false, false
        files = String[]
        index = 1
        while index <= length(args)
            arg = args[index]
            if arg in ("-h", "--help")
                print_usage(io)
                return 0
            elseif arg in ("-d", "--date", "--status", "--reason")
                index == length(args) && throw(ArgumentError("missing value for $arg"))
                index += 1
                value = args[index]
                arg in ("-d", "--date") ? (date = value) : arg == "--status" ? (status = value) : (reason = value)
            elseif arg == "--dry-run"
                dry_run = true
            elseif arg in ("-c", "--check")
                check = true
            elseif arg in ("-i", "--interactive")
                interactive = true
            elseif arg in ("-b", "--batch")
                nothing
            elseif startswith(arg, "-")
                throw(ArgumentError("unknown option: $arg"))
            else
                push!(files, arg)
            end
            index += 1
        end
        check && (dry_run || interactive) && throw(ArgumentError("--check cannot be combined with --dry-run or --interactive"))
        if check
            directory = document_dir(root, "archived")
            targets = isempty(files) ? (isdir(directory) ? sort(filter(f -> endswith(lowercase(f), ".md"), readdir(directory))) : String[]) : files
            invalid = String[]
            for file in targets
                try
                    path = resolve_document(root, "archived", file)
                    metadata_fields(read(path, String))
                catch err
                    push!(invalid, "$file: $(sprint(showerror, err))")
                end
            end
            println(io, "[archive-check] checked=$(length(targets)) invalid=$(length(invalid))")
            foreach(message -> println(io, "  ", message), invalid)
            return isempty(invalid) ? 0 : 1
        end
        if interactive
            isempty(files) || throw(ArgumentError("--interactive cannot be combined with explicit files"))
            directory = document_dir(root, "active")
            choices = isdir(directory) ? sort(filter(f -> endswith(lowercase(f), ".md"), readdir(directory))) : String[]
            isempty(choices) && (println(io, "No active documents"); return 0)
            foreach(pair -> println(io, "$(pair[1]). $(pair[2])"), enumerate(choices))
            println(io, "Select file numbers (comma-separated, or all):")
            selection = strip(readline(input))
            indices = selection == "all" ? collect(eachindex(choices)) : parse.(Int, strip.(split(selection, ',')))
            all(i -> i in eachindex(choices), indices) || throw(ArgumentError("invalid file selection"))
            files = choices[indices]
        end
        isempty(files) && throw(ArgumentError("no files specified"))
        plans = [prepare_archive(file; root, date, status, reason) for file in files]
        length(unique(path_key(plan.target) for plan in plans)) == length(plans) || throw(ArgumentError("duplicate archive destination in batch"))
        for plan in plans
            println(io, "$(dry_run ? "[dry-run]" : "[archive]") $(plan.source) -> $(plan.target)")
            dry_run || execute_archive(plan)
        end
        return 0
    catch err
        println(io, "[archive] FAILED: ", sprint(showerror, err))
        return 1
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(ArchiveDocs.main())
end
