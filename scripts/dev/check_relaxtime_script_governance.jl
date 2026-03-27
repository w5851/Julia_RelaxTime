module RelaxTimeScriptGovernance

using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MANIFEST_PATH = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "workflow", "classification_manifest.toml")

function _read(path::String)
    isfile(path) || error("governance manifest not found: $path")
    return TOML.parsefile(path)
end

function _has_direct_relaxtime_include(path::String)
    if !isfile(path)
        return false
    end
    for line in eachline(path)
        s = strip(line)
        occursin("include", s) || continue
        if occursin("\"src\"", s) && occursin("\"relaxtime\"", s)
            return true
        end
    end
    return false
end

function check_governance()::Bool
    manifest = _read(MANIFEST_PATH)
    scripts = get(manifest, "scripts", Dict{String,Any}())
    allowlist = Set(String.(get(get(manifest, "allowlist", Dict{String,Any}()), "public_authoritative_direct_relaxtime_includes", Any[])))

    violations = String[]
    for (rel, role_any) in scripts
        role = String(role_any)
        role == "public-authoritative" || continue
        rel_s = String(rel)
        rel_s in allowlist && continue
        abs_path = joinpath(PROJECT_ROOT, rel_s)
        if _has_direct_relaxtime_include(abs_path)
            push!(violations, rel_s)
        end
    end

    if !isempty(violations)
        println("[governance][FAIL] public-authoritative scripts must not include src/relaxtime directly")
        for v in violations
            println("  - ", v)
        end
        return false
    end

    println("[governance][PASS] relaxtime script classification and include policy are valid")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    ok = check_governance()
    ok || error("relaxtime script governance check failed")
end

end # module
