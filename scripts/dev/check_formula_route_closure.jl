#!/usr/bin/env julia

"""
    check_formula_route_closure.jl

Machine-readable governance check for formula-route closure packets.  A route
packet must identify its starting model, approximations, external sources,
formula-to-code tests, production authorization state, and unresolved items.
This checker validates documentation provenance only; it does not certify a
numerical result or a physical production route.
"""
module FormulaRouteClosure

using TOML

export validate_registry, main

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_REGISTRY_REL = joinpath("config", "governance", "formula_route_closure.toml")
const FORMULA_ROOT = "docs/reference/formula/"
const DOI_PATTERN = r"^10\.\d{4,9}/\S+$"i
const ARXIV_PATTERN = r"^(?:arxiv:)?(?:\d{4}\.\d{4,5}|[a-z][a-z.\-]+/\d{7})(?:v\d+)?$"i

normalize_rel(path::AbstractString) = replace(normpath(String(path)), '\\' => '/')

function _valid_source_identifier(source::AbstractString)
    value = strip(String(source))
    doi_value = replace(value, r"^https?://(?:dx\.)?doi\.org/"i => "")
    occursin(DOI_PATTERN, doi_value) && return true
    arxiv_value = replace(value, r"^https?://arxiv\.org/(?:abs|html)/"i => "")
    return occursin(ARXIV_PATTERN, arxiv_value)
end

function _nonempty_string(value, field::AbstractString, violations::Vector{String})
    if !(value isa AbstractString) || isempty(strip(String(value)))
        push!(violations, "$(field) must be a non-empty string")
        return nothing
    end
    return strip(String(value))
end

function _string_vector(table::AbstractDict, field::AbstractString, label::AbstractString, violations::Vector{String})
    raw = get(table, field, nothing)
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

function _existing_repo_file!(violations::Vector{String}, root::AbstractString, rel::AbstractString, label::AbstractString)
    normalized = normalize_rel(rel)
    if isabspath(normalized) || normalized == ".." || startswith(normalized, "../")
        push!(violations, "$(label) must be a repository-relative path: $(rel)")
        return nothing
    end
    path = joinpath(root, split(normalized, '/')...)
    if isfile(path)
        root_real = normpath(realpath(root))
        path_real = normpath(realpath(path))
        root_cmp = Sys.iswindows() ? lowercase(root_real) : root_real
        path_cmp = Sys.iswindows() ? lowercase(path_real) : path_real
        prefix = endswith(root_cmp, string(Base.Filesystem.path_separator)) ?
            root_cmp : root_cmp * string(Base.Filesystem.path_separator)
        if path_cmp != root_cmp && !startswith(path_cmp, prefix)
            push!(violations, "$(label) resolves outside repository root: $(normalized)")
            return nothing
        end
    else
        push!(violations, "$(label) does not exist: $(normalized)")
    end
    return path
end

function _load_registry(root::AbstractString, registry_rel::AbstractString, violations::Vector{String})
    normalized = normalize_rel(registry_rel)
    if isabspath(normalized) || normalized == ".." || startswith(normalized, "../")
        push!(violations, "formula-route registry must be a repository-relative path: $(registry_rel)")
        return nothing
    end
    path = joinpath(root, split(normalized, '/')...)
    if !isfile(path)
        push!(violations, "missing formula-route registry: $(normalized)")
        return nothing
    end
    root_real = normpath(realpath(root))
    path_real = normpath(realpath(path))
    root_cmp = Sys.iswindows() ? lowercase(root_real) : root_real
    path_cmp = Sys.iswindows() ? lowercase(path_real) : path_real
    prefix = endswith(root_cmp, string(Base.Filesystem.path_separator)) ?
        root_cmp : root_cmp * string(Base.Filesystem.path_separator)
    if path_cmp != root_cmp && !startswith(path_cmp, prefix)
        push!(violations, "formula-route registry resolves outside repository root: $(normalized)")
        return nothing
    end
    try
        return TOML.parsefile(path)
    catch err
        push!(violations, "failed to parse $(normalized): $(sprint(showerror, err))")
        return nothing
    end
end

function validate_registry(
    root::AbstractString=PROJECT_ROOT;
    registry_rel::AbstractString=DEFAULT_REGISTRY_REL,
)
    violations = String[]
    parsed = _load_registry(root, registry_rel, violations)
    parsed === nothing && return violations

    get(parsed, "schema_version", nothing) == "v1" ||
        push!(violations, "schema_version must be v1")
    allowed_statuses = Set(_string_vector(parsed, "allowed_statuses", "registry", violations))
    required_fields = _string_vector(parsed, "required_fields", "registry", violations)

    raw_routes = get(parsed, "route", nothing)
    if !(raw_routes isa AbstractVector) || isempty(raw_routes)
        push!(violations, "registry.route must be a non-empty array of tables")
        return violations
    end

    seen_ids = Set{String}()
    for (index, raw) in enumerate(raw_routes)
        label = "route[$(index)]"
        if !(raw isa AbstractDict)
            push!(violations, "$(label) must be a TOML table")
            continue
        end

        missing = filter(field -> !haskey(raw, field), required_fields)
        for field in missing
            push!(violations, "$(label) missing required field: $(field)")
        end

        id = _nonempty_string(get(raw, "id", nothing), "$(label).id", violations)
        id === nothing && continue
        route_label = "route[$(id)]"
        if id in seen_ids
            push!(violations, "duplicate formula-route id: $(id)")
        else
            push!(seen_ids, id)
        end

        status = _nonempty_string(get(raw, "status", nothing), "$(route_label).status", violations)
        if status !== nothing && !(status in allowed_statuses)
            push!(violations, "$(route_label).status is not allowed: $(status)")
        end

        document = _nonempty_string(get(raw, "document", nothing), "$(route_label).document", violations)
        document_path = nothing
        if document !== nothing
            normalized = normalize_rel(document)
            startswith(normalized, FORMULA_ROOT) ||
                push!(violations, "$(route_label).document must stay under $(FORMULA_ROOT): $(normalized)")
            document_path = _existing_repo_file!(violations, root, normalized, "$(route_label).document")
        end

        for field in ("scope", "model_start")
            _nonempty_string(get(raw, field, nothing), "$(route_label).$(field)", violations)
        end

        approximations = _string_vector(raw, "approximations", route_label, violations)
        isempty(approximations) && push!(violations, "$(route_label).approximations must not be empty")
        sources = _string_vector(raw, "external_sources", route_label, violations)
        isempty(sources) && push!(violations, "$(route_label).external_sources must not be empty")
        for source in sources
            _valid_source_identifier(source) || push!(violations,
                "$(route_label).external_sources has invalid DOI/arXiv identifier: $(source)")
        end
        tests = _string_vector(raw, "formula_tests", route_label, violations)
        isempty(tests) && push!(violations, "$(route_label).formula_tests must not be empty")
        unresolved = _string_vector(raw, "unresolved_items", route_label, violations)
        if status == "production_authorized"
            !isempty(unresolved) && push!(violations,
                "$(route_label).unresolved_items must be empty for production_authorized")
        else
            isempty(unresolved) && push!(violations, "$(route_label).unresolved_items must not be empty")
        end

        markers = _string_vector(raw, "required_document_markers", route_label, violations)
        isempty(markers) && push!(violations, "$(route_label).required_document_markers must not be empty")

        if id == "charged_rpa_bu_quark_only"
            for field in ("density_algorithms", "comparison_scheme", "bose_domain_policy")
                haskey(raw, field) || push!(violations,
                    "$(route_label) missing required field: $(field)")
            end
            density_algorithms = _string_vector(raw, "density_algorithms", route_label, violations)
            isempty(density_algorithms) && push!(violations, "$(route_label).density_algorithms must not be empty")
            required_density_algorithms = Set((
                "stable_particle_limit",
                "reduced_strict_bw",
                "q_pole_strict_bw",
                "phase_shift_bu",
            ))
            for algorithm in sort!(collect(setdiff(required_density_algorithms, Set(density_algorithms))))
                push!(violations, "$(route_label).density_algorithms missing required algorithm: $(algorithm)")
            end
            comparison_scheme = _nonempty_string(
                get(raw, "comparison_scheme", nothing),
                "$(route_label).comparison_scheme",
                violations,
            )
            if comparison_scheme !== nothing && comparison_scheme != "phase_shift_gbu_reference"
                push!(violations,
                    "$(route_label).comparison_scheme must be phase_shift_gbu_reference")
            end
            bose_policy = _nonempty_string(
                get(raw, "bose_domain_policy", nothing),
                "$(route_label).bose_domain_policy",
                violations,
            )
            if bose_policy !== nothing && bose_policy != "normal_phase_gate_x_min_cut_diagnostic"
                push!(violations,
                    "$(route_label).bose_domain_policy has unsupported value: $(bose_policy)")
            end
        end

        authorized = get(raw, "production_authorized", nothing)
        authorized isa Bool || push!(violations, "$(route_label).production_authorized must be boolean")
        if authorized === true && status != "production_authorized"
            push!(violations, "$(route_label) cannot authorize production while status=$(status)")
        end
        if status == "production_authorized" && authorized !== true
            push!(violations, "$(route_label) status=production_authorized requires production_authorized=true")
        end

        if document_path !== nothing && isfile(document_path)
            content = read(document_path, String)
            for marker in markers
                occursin(marker, content) ||
                    push!(violations, "$(route_label).document missing required marker: $(marker)")
            end
            for source in sources
                occursin(source, content) ||
                    push!(violations, "$(route_label).document does not cite registered source: $(source)")
            end
        end

        for test_path in tests
            _existing_repo_file!(violations, root, test_path, "$(route_label).formula_tests entry")
        end
    end

    return violations
end

function _registry_from_args(args::Vector{String})
    isempty(args) && return DEFAULT_REGISTRY_REL
    length(args) == 1 && startswith(args[1], "--registry=") ||
        error("usage: check_formula_route_closure.jl [--registry=<repository-relative-path>]")
    value = strip(split(args[1], "="; limit=2)[2])
    isempty(value) && error("--registry requires a non-empty path")
    return value
end

function main(args::Vector{String}=collect(String.(ARGS)))
    registry_rel = _registry_from_args(args)
    violations = validate_registry(PROJECT_ROOT; registry_rel=registry_rel)
    if !isempty(violations)
        println("[formula-route-closure] FAILED: $(length(violations)) violation(s)")
        for item in violations
            println(" - " * item)
        end
        return 1
    end
    println("[formula-route-closure] OK")
    println("  registry=$(normalize_rel(registry_rel))")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

end # module FormulaRouteClosure
