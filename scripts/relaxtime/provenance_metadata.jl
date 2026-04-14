module ProvenanceMetadata

using Dates
using SHA

export RunContext
export new_run_context
export write_run_sidecars
export csv_stats

Base.@kwdef struct RunContext
    run_id::String
    timestamp_utc::String
    script::String
    argv::Vector{String}
end

@inline function _json_escape(s::AbstractString)
    out = IOBuffer()
    for c in s
        if c == '"'
            print(out, "\\\"")
        elseif c == '\\'
            print(out, "\\\\")
        elseif c == '\n'
            print(out, "\\n")
        elseif c == '\r'
            print(out, "\\r")
        elseif c == '\t'
            print(out, "\\t")
        else
            print(out, c)
        end
    end
    return String(take!(out))
end

function _to_json(x)
    if x === nothing
        return "null"
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa Integer || x isa AbstractFloat
        return string(x)
    elseif x isa AbstractString
        return "\"$(_json_escape(x))\""
    elseif x isa AbstractVector
        return "[" * join((_to_json(v) for v in x), ",") * "]"
    elseif x isa AbstractDict
        pairs_sorted = sort(collect(pairs(x)); by=kv -> String(kv.first))
        parts = String[]
        for (k, v) in pairs_sorted
            push!(parts, "\"$(_json_escape(String(k)))\":" * _to_json(v))
        end
        return "{" * join(parts, ",") * "}"
    else
        return "\"$(_json_escape(string(x)))\""
    end
end

function _write_json(path::String, x)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, _to_json(x))
    end
end

function _current_git_commit(project_root::String)
    try
        return readchomp(`git -C $(project_root) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _current_git_branch(project_root::String)
    try
        return readchomp(`git -C $(project_root) rev-parse --abbrev-ref HEAD`)
    catch
        return "unknown"
    end
end

function _git_dirty(project_root::String)
    try
        return !isempty(readchomp(`git -C $(project_root) status --porcelain`))
    catch
        return false
    end
end

function _sha256_file(path::String)
    open(path, "r") do io
        return bytes2hex(sha256(read(io)))
    end
end

function _count_csv_rows(path::String)
    rows = 0
    seen_header = false
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !seen_header
                seen_header = true
                continue
            end
            rows += 1
        end
    end
    return rows
end

function _artifact_entry(path::String, project_root::String)
    abs_path = abspath(path)
    return Dict{String,Any}(
        "path" => relpath(abs_path, project_root),
        "sha256" => _sha256_file(abs_path),
        "rows" => endswith(lowercase(path), ".csv") ? _count_csv_rows(abs_path) : nothing,
    )
end

function new_run_context(script::String, argv::Vector{String}=copy(ARGS))::RunContext
    run_id = string(Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS"), "_", bytes2hex(rand(UInt8, 4)))
    timestamp_utc = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ")
    return RunContext(run_id=run_id, timestamp_utc=timestamp_utc, script=script, argv=copy(argv))
end

function csv_stats(path::String; converged_col::String="converged")
    header = String[]
    idx_conv = 0
    total = 0
    success = 0
    errors = 0
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if isempty(header)
                header = [strip(x) for x in split(s, ',')]
                idx_conv = something(findfirst(==(converged_col), header), 0)
                continue
            end
            total += 1
            if idx_conv > 0
                parts = split(s, ',')
                if idx_conv <= length(parts)
                    v = lowercase(strip(parts[idx_conv]))
                    if v == "true"
                        success += 1
                    else
                        errors += 1
                    end
                end
            end
        end
    end
    return (points_total=total, success_count=success, error_count=errors)
end

function write_run_sidecars(outdir::String;
    ctx::RunContext,
    effective_config::Dict{String,Any},
    artifacts::Vector{String},
    summary::Dict{String,Any}=Dict{String,Any}(),
)
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    effective = Dict{String,Any}(
        "schema_version" => "v1",
        "run_id" => ctx.run_id,
        "script" => ctx.script,
        "options" => effective_config,
    )
    effective_json = _to_json(effective)
    config_hash = bytes2hex(sha256(effective_json))

    project_toml = joinpath(project_root, "Project.toml")
    manifest_toml = joinpath(project_root, "Manifest.toml")

    manifest = Dict{String,Any}(
        "schema_version" => "v1",
        "run_id" => ctx.run_id,
        "timestamp_utc" => ctx.timestamp_utc,
        "script" => ctx.script,
        "argv" => ctx.argv,
        "cwd" => pwd(),
        "project_path" => project_root,
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
        "git_commit" => _current_git_commit(project_root),
        "git_branch" => _current_git_branch(project_root),
        "git_dirty" => _git_dirty(project_root),
        "project_toml_hash" => isfile(project_toml) ? _sha256_file(project_toml) : "",
        "manifest_toml_hash" => isfile(manifest_toml) ? _sha256_file(manifest_toml) : "",
        "config_hash" => config_hash,
        "artifacts" => [_artifact_entry(p, project_root) for p in artifacts],
        "summary" => summary,
    )

    _write_json(joinpath(outdir, "effective_config.json"), effective)
    _write_json(joinpath(outdir, "run_manifest.json"), manifest)
end

end # module ProvenanceMetadata
