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
