#!/usr/bin/env julia

"""
Models authority map machine-check guard.

Usage:
  julia --project=. scripts/dev/check_models_authority_map.jl
"""

using TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG_FILE = joinpath(ROOT, "config", "governance", "models_authority_map.toml")

normalize_rel(path::AbstractString) = replace(path, '\\' => '/')

function _require_nonempty_string(value, field::AbstractString)
    value isa AbstractString || error("$(field) must be a non-empty string, got $(typeof(value))")
    text = strip(String(value))
    isempty(text) && error("$(field) must be a non-empty string")
    return text
end

function _validate_required_markers(raw)
    raw isa AbstractVector || error("authority_map.required_markers must be a non-empty string array, got $(typeof(raw))")
    isempty(raw) && error("authority_map.required_markers must be a non-empty string array")

    markers = String[]
    seen = Set{String}()
    for (idx, item) in enumerate(raw)
        marker = _require_nonempty_string(item, "authority_map.required_markers[$idx]")
        marker in seen && error("authority_map.required_markers contains duplicate marker: $(marker)")
        push!(seen, marker)
        push!(markers, marker)
    end
    return markers
end

function main()
    isfile(CONFIG_FILE) || error("missing config file: $(normalize_rel(relpath(CONFIG_FILE, ROOT)))")

    parsed = TOML.parsefile(CONFIG_FILE)
    haskey(parsed, "authority_map") || error("missing [authority_map] section in $(normalize_rel(relpath(CONFIG_FILE, ROOT)))")

    authority_map = parsed["authority_map"]
    haskey(authority_map, "doc") || error("missing authority_map.doc in $(normalize_rel(relpath(CONFIG_FILE, ROOT)))")
    haskey(authority_map, "required_markers") || error("missing authority_map.required_markers in $(normalize_rel(relpath(CONFIG_FILE, ROOT)))")

    doc_rel = _require_nonempty_string(authority_map["doc"], "authority_map.doc")
    markers = _validate_required_markers(authority_map["required_markers"])

    doc_path = joinpath(ROOT, split(doc_rel, '/')...)
    isfile(doc_path) || error("missing authority doc: $(normalize_rel(relpath(doc_path, ROOT)))")

    content = read(doc_path, String)
    violations = String[]
    for marker in markers
        marker_text = String(marker)
        occursin(marker_text, content) || push!(violations, "missing marker in $(normalize_rel(relpath(doc_path, ROOT))): $(marker_text)")
    end

    if !isempty(violations)
        println("[models-authority-map] FAILED")
        for item in violations
            println(" - " * item)
        end
        exit(1)
    end

    println("[models-authority-map] OK")
    println("  doc=$(normalize_rel(relpath(doc_path, ROOT))) markers=$(length(markers))")
end

main()
