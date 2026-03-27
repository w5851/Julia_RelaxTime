module WorkflowConfig

export normalize_merge_validate

@inline function _to_dict_string_any(x)::Dict{String,Any}
    out = Dict{String,Any}()
    for (k, v) in pairs(x)
        ks = String(k)
        if v isa AbstractDict
            out[ks] = _to_dict_string_any(v)
        else
            out[ks] = v
        end
    end
    return out
end

@inline function _deep_merge(left::Dict{String,Any}, right::Dict{String,Any})::Dict{String,Any}
    out = Dict{String,Any}(left)
    for (k, rv) in right
        if haskey(out, k) && out[k] isa Dict{String,Any} && rv isa Dict{String,Any}
            out[k] = _deep_merge(out[k], rv)
        else
            out[k] = rv
        end
    end
    return out
end

@inline function _normalize_aliases(cfg::Dict{String,Any}, aliases::Dict{String,Any})::Dict{String,Any}
    out = Dict{String,Any}(cfg)
    proc_alias_map = get(get(aliases, "process_aliases", Dict{String,Any}()), "mappings", Dict{String,Any}())

    if haskey(out, "scan") && out["scan"] isa Dict{String,Any}
        scan = Dict{String,Any}(out["scan"])
        if haskey(scan, "cross_section") && scan["cross_section"] isa Dict{String,Any}
            xs = Dict{String,Any}(scan["cross_section"])
            if haskey(xs, "processes") && xs["processes"] isa AbstractVector
                xs["processes"] = [String(get(proc_alias_map, String(p), String(p))) for p in xs["processes"]]
            end
            scan["cross_section"] = xs
        end
        out["scan"] = scan
    end

    return out
end

@inline function _validate_source_shape(cfg::Dict{String,Any})
    return nothing
end

@inline function _validate_effective(cfg::Dict{String,Any})
    haskey(cfg, "scan") || throw(ArgumentError("effective config missing required key: scan"))
    return nothing
end

function normalize_merge_validate(default_cfg, toml_cfg, cli_cfg, aliases)
    default_cfg_d = _to_dict_string_any(default_cfg)
    toml_cfg_d = _to_dict_string_any(toml_cfg)
    cli_cfg_d = _to_dict_string_any(cli_cfg)
    aliases_d = _to_dict_string_any(aliases)

    trace = Symbol[]

    # raw(TOML+CLI) -> alias normalization -> source validation -> precedence merge -> effective validation
    push!(trace, :alias_normalization)
    toml_norm = _normalize_aliases(toml_cfg_d, aliases_d)
    cli_norm = _normalize_aliases(cli_cfg_d, aliases_d)

    push!(trace, :source_validation)
    _validate_source_shape(default_cfg_d)
    _validate_source_shape(toml_norm)
    _validate_source_shape(cli_norm)

    push!(trace, :precedence_merge)
    effective = _deep_merge(default_cfg_d, toml_norm)
    effective = _deep_merge(effective, cli_norm)

    push!(trace, :effective_validation)
    _validate_effective(effective)

    return (effective=effective, trace=trace)
end

end # module
