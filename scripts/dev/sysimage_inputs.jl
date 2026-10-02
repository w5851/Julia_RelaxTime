#!/usr/bin/env julia

module SysimageInputs

using SHA

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REQUIRED_FILES = (
    "Project.toml",
    "scripts/dev/build_sysimage.jl",
    "scripts/dev/precompile_workload.jl",
    "scripts/dev/sysimage_inputs.jl",
)

"""List source, dependency and configuration files used by the sysimage build."""
function input_files(root::AbstractString)
    files = String[]
    for rel in REQUIRED_FILES
        isfile(joinpath(root, rel)) || throw(ArgumentError("sysimage build input missing: $rel"))
        push!(files, rel)
    end
    isfile(joinpath(root, "Manifest.toml")) && push!(files, "Manifest.toml")
    for (directory, julia_only) in (("src", true), ("config/models", false), ("config/physics", false))
        full = joinpath(root, directory)
        isdir(full) || continue
        for (parent, _, names) in walkdir(full), name in names
            julia_only && !endswith(name, ".jl") && continue
            push!(files, replace(relpath(joinpath(parent, name), root), '\\' => '/'))
        end
    end
    return sort!(files)
end

"""Fingerprint file paths and bytes, including uncommitted additions/deletions."""
function fingerprint(root::AbstractString=ROOT)
    context = SHA.SHA2_256_CTX()
    SHA.update!(context, codeunits("sysimage-inputs-v1\n"))
    for rel in input_files(root)
        bytes = read(joinpath(root, rel))
        # Length prefixes keep both names and file boundaries unambiguous.
        SHA.update!(context, codeunits("$(ncodeunits(rel)):$rel:$(length(bytes)):"))
        SHA.update!(context, bytes)
    end
    return bytes2hex(SHA.digest!(context))
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ARGS == ["--fingerprint"] || error("Usage: julia scripts/dev/sysimage_inputs.jl --fingerprint")
    println(SysimageInputs.fingerprint())
end
