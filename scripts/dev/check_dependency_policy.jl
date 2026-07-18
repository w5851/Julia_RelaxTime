#!/usr/bin/env julia

module DependencyPolicyGovernance

using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_POLICY = joinpath(REPO_ROOT, "config", "ci", "dependency_policy.toml")

normalize_rel(path::AbstractString) = replace(normpath(path), '\\' => '/')

function _project_dependencies(project_path::AbstractString)
    isfile(project_path) || return nothing
    project = TOML.parsefile(project_path)
    return Dict{String,Any}(
        "deps" => get(project, "deps", Dict{String,Any}()),
        "weakdeps" => get(project, "weakdeps", Dict{String,Any}()),
        "extras" => get(project, "extras", Dict{String,Any}()),
        "compat" => get(project, "compat", Dict{String,Any}()),
    )
end

function _julia_files(root::AbstractString, scan_roots)
    files = String[]
    for relative_root in scan_roots
        directory = joinpath(root, relative_root)
        isdir(directory) || continue
        for (dir, _, names) in walkdir(directory)
            for name in names
                endswith(name, ".jl") || continue
                push!(files, joinpath(dir, name))
            end
        end
    end
    return sort(files)
end

function _runtime_reference_violations(root::AbstractString, files, dependency::String, symbols)
    violations = String[]
    import_pattern = Regex("^\\s*(?:using|import)\\s+" * dependency * "(?:\\s*:|\\s*\\z)")
    qualified_pattern = Regex("\\b" * dependency * "\\s*\\.")
    symbol_patterns = [Regex("\\b" * String(symbol) * "\\s*\\(") for symbol in symbols]

    for file in files
        relative = normalize_rel(relpath(file, root))
        for (line_number, line) in enumerate(eachline(file))
            stripped = strip(line)
            startswith(stripped, '#') && continue
            if occursin(import_pattern, line) || occursin(qualified_pattern, line) ||
               any(pattern -> occursin(pattern, line), symbol_patterns)
                push!(violations, "$(relative):$(line_number) references forbidden root dependency $(dependency)")
            end
        end
    end
    return violations
end

function validate_repository(root::AbstractString; policy_path::AbstractString=joinpath(root, "config", "ci", "dependency_policy.toml"))
    violations = String[]
    isfile(policy_path) || return ["missing dependency policy: $(normalize_rel(relpath(policy_path, root)))"]

    policy = TOML.parsefile(policy_path)
    get(policy, "schema_version", 0) == 1 || push!(violations, "dependency policy schema_version must be 1")

    root_project_rel = String(get(policy, "root_project", "Project.toml"))
    root_project_path = joinpath(root, root_project_rel)
    root_project = _project_dependencies(root_project_path)
    if root_project === nothing
        push!(violations, "missing root project: $(normalize_rel(root_project_rel))")
        return violations
    end

    forbidden = String.(get(policy, "forbidden_root_direct_dependencies", String[]))
    isolated = get(policy, "isolated_dependencies", Dict{String,Any}())
    scan_roots = String.(get(policy, "runtime_scan_roots", String[]))
    runtime_files = _julia_files(root, scan_roots)

    for dependency in forbidden
        for section in ("deps", "weakdeps", "extras", "compat")
            haskey(root_project[section], dependency) && push!(violations,
                "$(dependency) must not appear in root $(root_project_rel) [$(section)]")
        end

        dependency_policy = get(isolated, dependency, Dict{String,Any}())
        symbols = String.(get(dependency_policy, "forbidden_symbols", String[]))
        append!(violations, _runtime_reference_violations(root, runtime_files, dependency, symbols))

        isolated_project_rel = get(dependency_policy, "project", nothing)
        if isolated_project_rel === nothing
            push!(violations, "$(dependency) is missing an isolated project declaration")
            continue
        end
        isolated_project = _project_dependencies(joinpath(root, String(isolated_project_rel)))
        if isolated_project === nothing
            push!(violations, "missing isolated project for $(dependency): $(isolated_project_rel)")
        elseif !haskey(isolated_project["deps"], dependency)
            push!(violations, "$(dependency) must be declared in isolated project $(isolated_project_rel)")
        end
    end

    return violations
end

function main()
    violations = validate_repository(REPO_ROOT; policy_path=DEFAULT_POLICY)
    if !isempty(violations)
        println("[dependency-policy] FAILED: $(length(violations)) violation(s)")
        foreach(item -> println(" - " * item), violations)
        return 1
    end
    println("[dependency-policy] OK")
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(DependencyPolicyGovernance.main())
end
