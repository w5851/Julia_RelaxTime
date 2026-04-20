using SHA

@inline function adapter_repo_root()
    return normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
end

function adapter_git_commit()
    try
        return readchomp(`git -C $(adapter_repo_root()) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function adapter_config_hash(kind::Symbol, kwargs::NamedTuple)
    payload = String[]
    push!(payload, String(kind))
    for key in sort(collect(keys(kwargs)); by=String)
        push!(payload, String(key) * "=" * sprint(show, getproperty(kwargs, key)))
    end
    return bytes2hex(SHA.sha2_256(join(payload, "|")))
end

function resolve_adapter_output_dir(
    kwargs::NamedTuple,
    key_priority::Tuple{Vararg{Symbol}};
    path_keys::Tuple{Vararg{Symbol}}=(:output_path,),
)
    for key in key_priority
        hasproperty(kwargs, key) || continue
        value = String(getproperty(kwargs, key))
        if key in path_keys
            return dirname(value)
        end
        return value
    end
    return mktempdir()
end

function normalize_adapter_kwargs(kwargs::NamedTuple, alias_map::Dict{String, Symbol})
    normalized_pairs = Pair{Symbol, Any}[]
    canonical_present = Set{Symbol}(keys(kwargs))
    emitted_keys = Set{Symbol}()

    for key in keys(kwargs)
        value = getproperty(kwargs, key)
        canonical_key = get(alias_map, String(key), key)

        if canonical_key != key && canonical_key in canonical_present
            continue
        end
        canonical_key in emitted_keys && continue

        push!(normalized_pairs, canonical_key => value)
        push!(emitted_keys, canonical_key)
    end

    return (; normalized_pairs...)
end

function classify_pipeline_error(err)::Symbol
    if err isa ArgumentError
        return :input_validation_error
    end
    if err isa Base.IOError
        return :artifact_io_error
    end

    error_message = sprint(showerror, err)
    if occursin(r"converge"i, error_message)
        return :numerical_convergence_error
    end

    return :unexpected_error
end
