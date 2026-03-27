module WorkflowConfigAudit

export flatten_effective_keys, build_consumption_report

function flatten_effective_keys(cfg; prefix::String="")::Set{String}
    keys_out = Set{String}()
    for (k, v) in pairs(cfg)
        path = isempty(prefix) ? String(k) : string(prefix, ".", k)
        if v isa AbstractDict
            union!(keys_out, flatten_effective_keys(v; prefix=path))
        else
            push!(keys_out, path)
        end
    end
    return keys_out
end

function build_consumption_report(
    effective_cfg,
    consumed_keys::Set{String};
    overridden::Set{String}=Set{String}(),
    fallback_used::Bool=false,
    strict::Bool=false,
)
    all_keys = flatten_effective_keys(effective_cfg)
    unused = sort(collect(setdiff(all_keys, consumed_keys)))
    consumed = sort(collect(intersect(all_keys, consumed_keys)))
    overridden_sorted = sort(collect(overridden))

    if strict && !isempty(unused)
        throw(ArgumentError("strict mode: unconsumed keys detected: $(join(unused, ", "))"))
    end

    return Dict{String,Any}(
        "consumed_keys" => consumed,
        "unused_keys" => unused,
        "overridden_keys" => overridden_sorted,
        "fallback_used" => fallback_used,
    )
end

end # module
