#!/usr/bin/env julia

"""
Validate repository-owned Codex skills and optionally audit duplicate registrations.

Usage:
  julia --project=. scripts/dev/check_skill_governance.jl
  julia --project=. scripts/dev/check_skill_governance.jl --extra-root C:\\Users\\name\\.agents\\skills
  julia --project=. scripts/dev/check_skill_governance.jl --extra-root PATH --strict-duplicates
"""

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SKILLS_ROOT = joinpath(ROOT, ".agents", "skills")
const ALLOWED_FRONTMATTER_KEYS = Set(["name", "description"])
const TRIGGER_HEADING_RE = r"(?im)^##\s+(?:when\s+to(?:\s+use|\s+apply)?|apply\s+when|何时使用|适用场景|触发与启动|[0-9]+\)\s*适用场景).*$"
const MARKDOWN_LINK_RE = r"\[[^\]]+\]\(([^)]+)\)"

normalize_text(text::AbstractString) = replace(text, "\r\n" => "\n", "\r" => "\n")

function parse_args(args)
    extra_roots = String[]
    strict_duplicates = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--extra-root"
            i == length(args) && error("--extra-root requires a path")
            i += 1
            push!(extra_roots, abspath(args[i]))
        elseif startswith(arg, "--extra-root=")
            push!(extra_roots, abspath(split(arg, '='; limit=2)[2]))
        elseif arg == "--strict-duplicates"
            strict_duplicates = true
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/check_skill_governance.jl [--extra-root PATH] [--strict-duplicates]")
            return nothing
        else
            error("unknown option: $(arg)")
        end
        i += 1
    end
    return (; extra_roots, strict_duplicates)
end

function skill_files(root::String)
    files = Dict{String,String}()
    isdir(root) || return files
    for entry in sort(readdir(root))
        skill_md = joinpath(root, entry, "SKILL.md")
        isfile(skill_md) || continue
        files[entry] = skill_md
    end
    return files
end

function extract_frontmatter(content::String)
    lines = split(normalize_text(content), '\n'; keepempty=true)
    isempty(lines) && error("empty SKILL.md")
    lines[1] == "---" || error("frontmatter must start on line 1")
    closing = findnext(==("---"), lines, 2)
    closing === nothing && error("frontmatter closing delimiter is missing")
    closing == 2 && error("frontmatter is empty")
    return lines[2:closing-1], join(lines[closing+1:end], "\n")
end

function top_level_frontmatter(frontmatter_lines)
    fields = Dict{String,String}()
    keys = String[]
    for line in frontmatter_lines
        (!isempty(line) && isspace(first(line))) && continue
        m = match(r"^([A-Za-z0-9_-]+):\s*(.*)$", line)
        m === nothing && continue
        key = m.captures[1]
        value = strip(m.captures[2])
        push!(keys, key)
        fields[key] = value
    end
    return keys, fields
end

function unquote_yaml_scalar(value::String)
    length(value) >= 2 || return value
    if (startswith(value, '"') && endswith(value, '"')) ||
       (startswith(value, '\'') && endswith(value, '\''))
        return value[2:end-1]
    end
    return value
end

function local_link_target(raw_target::AbstractString)
    target = strip(raw_target)
    if startswith(target, "<") && endswith(target, ">")
        target = target[2:end-1]
    end
    target = first(split(target, '#'; limit=2))
    isempty(target) && return nothing
    lowercase_target = lowercase(target)
    startswith(lowercase_target, "http://") && return nothing
    startswith(lowercase_target, "https://") && return nothing
    startswith(lowercase_target, "mailto:") && return nothing
    isabspath(target) && return nothing
    return target
end

function validate_skill(folder_name::String, skill_md::String)
    violations = String[]
    warnings = String[]
    content = read(skill_md, String)

    frontmatter_lines = String[]
    body = ""
    try
        frontmatter_lines, body = extract_frontmatter(content)
    catch err
        push!(violations, "$(folder_name): $(sprint(showerror, err))")
        return violations, warnings
    end

    keys, fields = top_level_frontmatter(frontmatter_lines)
    key_set = Set(keys)
    key_set == ALLOWED_FRONTMATTER_KEYS || push!(violations,
        "$(folder_name): frontmatter keys must be exactly name, description; found $(join(sort(collect(key_set)), ", "))")
    length(keys) == length(key_set) || push!(violations, "$(folder_name): duplicate frontmatter key")

    declared_name = unquote_yaml_scalar(get(fields, "name", ""))
    declared_name == folder_name || push!(violations,
        "$(folder_name): frontmatter name is $(repr(declared_name))")

    description = unquote_yaml_scalar(get(fields, "description", ""))
    isempty(strip(description)) && push!(violations, "$(folder_name): description is empty")
    occursin("TODO", description) && push!(violations, "$(folder_name): description still contains TODO")

    occursin(TRIGGER_HEADING_RE, body) && push!(violations,
        "$(folder_name): move positive trigger headings into the frontmatter description")

    for m in eachmatch(MARKDOWN_LINK_RE, body)
        target = local_link_target(m.captures[1])
        target === nothing && continue
        resolved = normpath(joinpath(dirname(skill_md), target))
        ispath(resolved) || push!(violations,
            "$(folder_name): missing linked resource $(target)")
    end

    openai_yaml = joinpath(dirname(skill_md), "agents", "openai.yaml")
    if !isfile(openai_yaml)
        push!(violations, "$(folder_name): missing agents/openai.yaml")
    else
        yaml = read(openai_yaml, String)
        occursin("default_prompt:", yaml) || push!(violations,
            "$(folder_name): agents/openai.yaml is missing default_prompt")
        occursin("\$" * folder_name, yaml) || push!(violations,
            "$(folder_name): default_prompt must mention \$$(folder_name)")
    end

    if folder_name == "xjtu-institution-literature-access" && isfile(openai_yaml)
        yaml = read(openai_yaml, String)
        occursin("allow_implicit_invocation: false", yaml) || push!(violations,
            "$(folder_name): institution access must disable implicit invocation")
    end

    if folder_name == "mcp-inspire" && isfile(openai_yaml)
        yaml = read(openai_yaml, String)
        occursin("value: \"mcp-inspire\"", yaml) || push!(violations,
            "$(folder_name): declare the MCP tool dependency in agents/openai.yaml")
    end

    line_count = count(==('\n'), normalize_text(content)) + 1
    line_count > 200 && push!(warnings, "$(folder_name): SKILL.md has $(line_count) lines; consider progressive disclosure")
    return violations, warnings
end

function audit_external_duplicates(repo_skills, extra_roots, strict_duplicates)
    violations = String[]
    warnings = String[]
    for root in extra_roots
        if !isdir(root)
            push!(warnings, "extra skill root does not exist: $(root)")
            continue
        end
        for (name, external_file) in skill_files(root)
            haskey(repo_skills, name) || continue
            repo_file = repo_skills[name]
            same = normalize_text(read(repo_file, String)) == normalize_text(read(external_file, String))
            detail = same ? "identical duplicate registration" : "conflicting duplicate registration"
            message = "$(name): $(detail) in $(root)"
            if strict_duplicates
                push!(violations, message)
            else
                push!(warnings, message)
            end
        end
    end
    return violations, warnings
end

function main()
    parsed = parse_args(ARGS)
    parsed === nothing && return

    repo_skills = skill_files(SKILLS_ROOT)
    isempty(repo_skills) && error("no skills found under $(SKILLS_ROOT)")

    violations = String[]
    warnings = String[]
    for name in sort(collect(keys(repo_skills)))
        skill_violations, skill_warnings = validate_skill(name, repo_skills[name])
        append!(violations, skill_violations)
        append!(warnings, skill_warnings)
    end

    duplicate_violations, duplicate_warnings = audit_external_duplicates(
        repo_skills,
        parsed.extra_roots,
        parsed.strict_duplicates,
    )
    append!(violations, duplicate_violations)
    append!(warnings, duplicate_warnings)

    if !isempty(violations)
        println("[skill-governance] FAILED: $(length(violations)) violation(s)")
        for item in violations
            println(" - " * item)
        end
        if !isempty(warnings)
            println("[skill-governance] WARNINGS: $(length(warnings))")
            for item in warnings
                println(" - " * item)
            end
        end
        exit(1)
    end

    println("[skill-governance] OK")
    println("  skills=$(length(repo_skills))")
    if !isempty(warnings)
        println("[skill-governance] WARNINGS: $(length(warnings))")
        for item in warnings
            println(" - " * item)
        end
    end
end

main()
