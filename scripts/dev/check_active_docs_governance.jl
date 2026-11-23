#!/usr/bin/env julia

module ActiveDocsGovernance

using Dates
using TOML
using Unicode

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const NAME_RE = r"^\d{4}-\d{2}-\d{2}_.+\.md$"
const REVIEW_AGE_DAYS = 60
const ALLOWLIST_DISPLAY = "config/governance/retired_docs_allowlist.toml"
const ALLOWLIST_RELPATH = joinpath(split(ALLOWLIST_DISPLAY, '/')...)
const ALLOWLIST_SCHEMA_VERSION = 1
const RETIRED_STATUS = "retired_readonly"

_posix_relpath(path::AbstractString) = Unicode.normalize(replace(path, '\\' => '/'), :NFC)

"""Enumerate every regular file below `dir` as sorted POSIX-style relative paths."""
function collect_relative_files(dir::AbstractString)
    files = String[]
    isdir(dir) || return files
    for (dirpath, _, names) in walkdir(dir)
        for name in names
            full = joinpath(dirpath, name)
            isfile(full) || continue
            push!(files, _posix_relpath(relpath(full, dir)))
        end
    end
    return sort!(files)
end

"""Parse the retired docs allowlist into `Dict(dir => Vector(allowed_files))`."""
function load_retired_allowlist(path::AbstractString)
    isfile(path) || throw(ArgumentError("missing retired docs allowlist: $(path)"))
    config = try
        TOML.parsefile(path)
    catch err
        throw(ArgumentError("unreadable retired docs allowlist: $(path) ($(err))"))
    end
    get(config, "schema_version", nothing) == ALLOWLIST_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported schema_version in $(path); expected $(ALLOWLIST_SCHEMA_VERSION)"))
    directories = get(config, "directories", nothing)
    directories isa AbstractDict && !isempty(directories) ||
        throw(ArgumentError("retired docs allowlist $(path) defines no [directories.*] table"))
    allowlist = Dict{String,Vector{String}}()
    for (relative_dir, entry) in directories
        entry isa AbstractDict || throw(ArgumentError("invalid entry for $(relative_dir) in $(path)"))
        status = get(entry, "status", RETIRED_STATUS)
        status == RETIRED_STATUS ||
            throw(ArgumentError("unsupported status=$(status) for $(relative_dir); expected $(RETIRED_STATUS)"))
        files = get(entry, "allowed_files", nothing)
        files isa AbstractVector && all(item -> item isa AbstractString, files) ||
            throw(ArgumentError("allowed_files for $(relative_dir) must be an array of strings"))
        normalized = _posix_relpath.(String.(files))
        length(unique(normalized)) == length(normalized) ||
            throw(ArgumentError("duplicate entries in allowed_files for $(relative_dir)"))
        allowlist[_posix_relpath(relative_dir)] = sort!(normalized)
    end
    return allowlist
end

"""Compare retired directories against the allowlist: new files violate, missing files advise."""
function review_retired_directories(root::AbstractString=ROOT; allowlist_path::AbstractString=joinpath(root, ALLOWLIST_RELPATH))
    violations = String[]
    advisories = String[]
    allowlist = try
        load_retired_allowlist(allowlist_path)
    catch err
        push!(violations, "invalid retired docs allowlist: " * sprint(showerror, err))
        return (; violations, advisories)
    end
    for (relative_dir, allowed) in sort(collect(allowlist); by=first)
        actual = collect_relative_files(joinpath(root, split(relative_dir, '/')...))
        allowed_set = Set(allowed)
        for file in actual
            file in allowed_set && continue
            push!(violations,
                "new file in retired directory: $(relative_dir)/$(file) is not in $(ALLOWLIST_DISPLAY) " *
                "(status=$(RETIRED_STATUS); move it to an active domain per AGENTS.md placement rules — " *
                "do NOT add new docs to the allowlist)")
        end
        actual_set = Set(actual)
        for file in allowed
            file in actual_set && continue
            push!(advisories,
                "retired doc missing from disk: $(relative_dir)/$(file) (listed in $(ALLOWLIST_DISPLAY); " *
                "verify issue130 audit tables and figure_asset_registry provenance before removing the entry)")
        end
    end
    return (; violations, advisories)
end

function review_documents(root::AbstractString=ROOT; current_date::Date=today(), check_retired::Bool=true)
    active_dir = joinpath(root, "docs", "dev", "active")
    violations = String[]
    advisories = String[]
    # Git omits this directory when every active task has been archived.
    if ispath(active_dir)
        for file in sort(readdir(active_dir))
            isfile(joinpath(active_dir, file)) && endswith(file, ".md") || continue
            if !occursin(NAME_RE, file)
                push!(violations, "invalid name format: $(file) (expected YYYY-MM-DD_*.md)")
                continue
            end
            created = try
                Date(first(split(file, '_')), dateformat"yyyy-mm-dd")
            catch
                push!(violations, "invalid date in active document filename: $(file)")
                continue
            end
            if Dates.value(current_date - created) > REVIEW_AGE_DAYS
                push!(advisories, "review active doc (>$(REVIEW_AGE_DAYS)d): $(file); archive by task status, not age")
            end
        end
    end
    if check_retired
        allowlist_path = joinpath(root, ALLOWLIST_RELPATH)
        if isfile(allowlist_path)
            retired = review_retired_directories(root; allowlist_path)
            append!(violations, retired.violations)
            append!(advisories, retired.advisories)
        elseif normpath(root) == normpath(ROOT)
            # Never let a missing allowlist silently disable the gate in the real repository.
            push!(violations, "invalid retired docs allowlist: missing $(ALLOWLIST_DISPLAY)")
        end
    end
    return (; violations, advisories)
end

function main()
    result = review_documents()
    foreach(item -> println("[active-docs-governance] advisory: " * item), result.advisories)
    if !isempty(result.violations)
        println("[active-docs-governance] FAILED")
        foreach(item -> println(" - " * item), result.violations)
        return 1
    end
    println("[active-docs-governance] OK")
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(ActiveDocsGovernance.main())
end
