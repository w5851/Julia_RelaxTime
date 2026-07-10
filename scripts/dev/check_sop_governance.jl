#!/usr/bin/env julia

module SopGovernance

using Dates
using TOML

export validate_registry, main

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_REGISTRY_REL = joinpath("config", "governance", "docs_authority_map.toml")
const SOP_ROOT_REL = replace(joinpath("docs", "guides", "sop"), '\\' => '/') * "/"

normalize_rel(path::AbstractString) = replace(normpath(String(path)), '\\' => '/')

function _nonempty_string(value, field::AbstractString, violations::Vector{String})
    if !(value isa AbstractString) || isempty(strip(String(value)))
        push!(violations, "$(field) must be a non-empty string")
        return nothing
    end
    return strip(String(value))
end

function _string_vector(table::AbstractDict, field::AbstractString, label::AbstractString, violations::Vector{String})
    if !haskey(table, field)
        push!(violations, "$(label) missing field: $(field)")
        return String[]
    end
    raw = table[field]
    if !(raw isa AbstractVector)
        push!(violations, "$(label).$(field) must be a string array")
        return String[]
    end

    values = String[]
    seen = Set{String}()
    for (index, item) in enumerate(raw)
        value = _nonempty_string(item, "$(label).$(field)[$(index)]", violations)
        value === nothing && continue
        if value in seen
            push!(violations, "$(label).$(field) contains duplicate value: $(value)")
            continue
        end
        push!(seen, value)
        push!(values, value)
    end
    return values
end

function _existing_repo_path!(violations::Vector{String}, root::AbstractString, rel::AbstractString, label::AbstractString)
    normalized = normalize_rel(rel)
    if isabspath(normalized) || startswith(normalized, "../")
        push!(violations, "$(label) must be a repository-relative path: $(rel)")
        return nothing
    end
    path = joinpath(root, split(normalized, '/')...)
    ispath(path) || push!(violations, "$(label) does not exist: $(normalized)")
    return path
end

function _parse_date(value, label::AbstractString, violations::Vector{String})
    text = _nonempty_string(value, label, violations)
    text === nothing && return nothing
    try
        return Date(text, dateformat"yyyy-mm-dd")
    catch
        push!(violations, "$(label) must use YYYY-MM-DD: $(text)")
        return nothing
    end
end

function _load_registry(root::AbstractString, registry_rel::AbstractString, violations::Vector{String})
    registry_path = joinpath(root, split(normalize_rel(registry_rel), '/')...)
    if !isfile(registry_path)
        push!(violations, "missing SOP authority registry: $(normalize_rel(registry_rel))")
        return nothing
    end
    try
        return TOML.parsefile(registry_path)
    catch err
        push!(violations, "failed to parse $(normalize_rel(registry_rel)): $(sprint(showerror, err))")
        return nothing
    end
end

function validate_registry(
    root::AbstractString=PROJECT_ROOT;
    registry_rel::AbstractString=DEFAULT_REGISTRY_REL,
    current_date::Date=today(),
)
    violations = String[]
    parsed = _load_registry(root, registry_rel, violations)
    parsed === nothing && return violations

    get(parsed, "schema_version", nothing) == "v1" ||
        push!(violations, "schema_version must be v1")

    allowed_statuses = Set(_string_vector(parsed, "allowed_statuses", "registry", violations))
    required_sections = _string_vector(parsed, "required_sections", "registry", violations)
    forbidden_patterns = _string_vector(parsed, "forbidden_patterns", "registry", violations)

    index_rel = _nonempty_string(
        get(parsed, "stable_entrypoint_index", nothing),
        "registry.stable_entrypoint_index",
        violations,
    )
    index_content = ""
    if index_rel !== nothing
        index_path = _existing_repo_path!(violations, root, index_rel, "stable entrypoint index")
        if index_path !== nothing && isfile(index_path)
            index_content = read(index_path, String)
        end
    end

    raw_sops = get(parsed, "sop", nothing)
    if !(raw_sops isa AbstractVector) || isempty(raw_sops)
        push!(violations, "registry.sop must be a non-empty array of tables")
        return violations
    end

    seen_ids = Set{String}()
    authority_owner = Dict{String, String}()

    for (index, raw) in enumerate(raw_sops)
        label = "sop[$(index)]"
        if !(raw isa AbstractDict)
            push!(violations, "$(label) must be a TOML table")
            continue
        end

        id = _nonempty_string(get(raw, "id", nothing), "$(label).id", violations)
        path_rel = _nonempty_string(get(raw, "path", nothing), "$(label).path", violations)
        status = _nonempty_string(get(raw, "status", nothing), "$(label).status", violations)
        item_label = id === nothing ? label : "sop[$(id)]"

        if id !== nothing
            if id in seen_ids
                push!(violations, "duplicate SOP id: $(id)")
            else
                push!(seen_ids, id)
            end
        end

        if status !== nothing && !(status in allowed_statuses)
            push!(violations, "$(item_label).status is not allowed: $(status)")
        end

        authorities = _string_vector(raw, "authoritative_for", item_label, violations)
        isempty(authorities) && push!(violations, "$(item_label).authoritative_for must not be empty")
        if status != "deprecated"
            for authority in authorities
                if haskey(authority_owner, authority)
                    push!(violations, "authoritative_for '$(authority)' is claimed by both $(authority_owner[authority]) and $(item_label)")
                else
                    authority_owner[authority] = item_label
                end
            end
        end

        entrypoints = _string_vector(raw, "entrypoints", item_label, violations)
        stable_entrypoints = _string_vector(raw, "stable_entrypoints", item_label, violations)
        configs = _string_vector(raw, "configs", item_label, violations)
        commands = _string_vector(raw, "verification_commands", item_label, violations)
        status == "active" && isempty(commands) &&
            push!(violations, "$(item_label).verification_commands must not be empty for active SOP")

        for rel in vcat(entrypoints, configs)
            _existing_repo_path!(violations, root, rel, "$(item_label) referenced path")
        end
        for rel in stable_entrypoints
            _existing_repo_path!(violations, root, rel, "$(item_label) stable entrypoint")
            isempty(index_content) || occursin(normalize_rel(rel), normalize_rel(index_content)) ||
                push!(violations, "$(item_label) stable entrypoint is absent from $(index_rel): $(normalize_rel(rel))")
        end

        sop_path = nothing
        if path_rel !== nothing
            normalized_path = normalize_rel(path_rel)
            startswith(normalized_path, SOP_ROOT_REL) ||
                push!(violations, "$(item_label).path must stay under docs/guides/sop/: $(normalized_path)")
            sop_path = _existing_repo_path!(violations, root, normalized_path, "$(item_label) document")
        end

        review_days = get(raw, "review_cycle_days", nothing)
        if !(review_days isa Integer) || review_days <= 0
            push!(violations, "$(item_label).review_cycle_days must be a positive integer")
            review_days = nothing
        end
        last_verified = _parse_date(get(raw, "last_verified", nothing), "$(item_label).last_verified", violations)
        if status == "active" && review_days !== nothing && last_verified !== nothing
            if last_verified > current_date
                push!(violations, "$(item_label).last_verified is in the future: $(last_verified)")
            elseif Dates.value(current_date - last_verified) > review_days
                push!(violations, "$(item_label) review is overdue: last_verified=$(last_verified), review_cycle_days=$(review_days)")
            end
        end

        if status == "active" && sop_path !== nothing && isfile(sop_path)
            content = read(sop_path, String)
            for section in required_sections
                occursin(section, content) ||
                    push!(violations, "$(item_label) missing required section: $(section)")
            end
            for pattern in forbidden_patterns
                occursin(pattern, content) &&
                    push!(violations, "$(item_label) contains forbidden pattern: $(pattern)")
            end
        end
    end

    return violations
end

function main(args::Vector{String}=collect(String.(ARGS)))
    isempty(args) || error("check_sop_governance.jl does not accept arguments")
    violations = validate_registry()
    if !isempty(violations)
        println("[sop-governance] FAILED: $(length(violations)) violation(s)")
        for item in violations
            println(" - " * item)
        end
        return 1
    end

    println("[sop-governance] OK")
    println("  registry=$(normalize_rel(DEFAULT_REGISTRY_REL))")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

end # module SopGovernance
