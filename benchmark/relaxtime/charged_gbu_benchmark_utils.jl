module ChargedGBUBenchmarkUtils

using Statistics

export parse_points, validate_settings, summarize_timed, timing_record, measure_call

function parse_points(spec::AbstractString)
    entries = split(strip(spec), ',')
    isempty(entries) && throw(ArgumentError("benchmark point list must not be empty"))
    points = Tuple{Float64,Float64}[]
    for entry in entries
        fields = split(strip(entry), ':')
        length(fields) == 2 || throw(ArgumentError(
            "point must be T_MeV:muB_MeV, got $(entry)"))
        T_MeV = parse(Float64, strip(fields[1]))
        muB_MeV = parse(Float64, strip(fields[2]))
        isfinite(T_MeV) && T_MeV > 0.0 || throw(ArgumentError(
            "T_MeV must be finite and positive, got $(T_MeV)"))
        isfinite(muB_MeV) && muB_MeV >= 0.0 || throw(ArgumentError(
            "muB_MeV must be finite and nonnegative, got $(muB_MeV)"))
        push!(points, (T_MeV, muB_MeV))
    end
    return points
end

function validate_settings(settings)
    settings.mesh >= 16 && settings.cut_nodes >= 8 && settings.tail_nodes >= 16 &&
        settings.omega_nodes >= 16 && settings.q_nodes >= 4 &&
        isfinite(settings.qmax) && settings.qmax > 0.0 ||
        throw(ArgumentError("invalid infinite-kernel screening settings"))
    settings.repeats >= 1 && settings.repeats <= 100 ||
        throw(ArgumentError("repeats must be between 1 and 100"))
    settings.density_repeats >= 1 && settings.density_repeats <= 100 ||
        throw(ArgumentError("density_repeats must be between 1 and 100"))
    return settings
end

function _numeric_values(samples, field)
    values = Float64[]
    for sample in samples
        value = getproperty(sample, field)
        value === nothing && continue
        isfinite(Float64(value)) || continue
        push!(values, Float64(value))
    end
    return values
end

function _summary(values)
    isempty(values) && return (count=0, mean=nothing, median=nothing,
        min=nothing, max=nothing)
    return (count=length(values), mean=mean(values), median=median(values),
        min=minimum(values), max=maximum(values))
end

function summarize_timed(samples)
    successes = [sample for sample in samples if sample.success]
    failures = [sample.error for sample in samples if !sample.success]
    return (count=length(samples), success_count=length(successes),
        failure_count=length(failures), elapsed_s=_summary(_numeric_values(successes, :elapsed_s)),
        timed_s=_summary(_numeric_values(successes, :timed_s)),
        bytes=_summary(_numeric_values(successes, :bytes)),
        gc_s=_summary(_numeric_values(successes, :gc_s)),
        compile_s=_summary(_numeric_values(successes, :compile_s)),
        recompile_s=_summary(_numeric_values(successes, :recompile_s)),
        compilation_observed_count=count(s -> s.compile_s !== nothing && s.compile_s > 0, successes),
        errors=failures)
end

function timing_record(sample; include_value=false)
    record = (success=sample.success, elapsed_s=sample.elapsed_s, timed_s=sample.timed_s,
        bytes=sample.bytes, gc_s=sample.gc_s, compile_s=sample.compile_s,
        recompile_s=sample.recompile_s, error=sample.error)
    return include_value ? merge(record, (value=sample.value,)) : record
end

"""Measure one invocation; cancellation is never converted into a failed sample.

Compiler fields describe work observed inside @timed. Caller inference may take
place earlier, so these fields are not a complete process-level JIT accounting.
"""
function measure_call(f)
    started = Base.time_ns()
    try
        timed = @timed f()
        optional(key) = hasproperty(timed, key) ? Float64(getproperty(timed, key)) : nothing
        return (success=true, value=timed.value,
            elapsed_s=(Base.time_ns() - started) / 1.0e9,
            timed_s=Float64(timed.time), bytes=optional(:bytes), gc_s=optional(:gctime),
            compile_s=optional(:compile_time), recompile_s=optional(:recompile_time), error=nothing)
    catch err
        err isa InterruptException && rethrow()
        return (success=false, value=nothing, elapsed_s=(Base.time_ns() - started) / 1.0e9,
            timed_s=nothing, bytes=nothing, gc_s=nothing, compile_s=nothing,
            recompile_s=nothing, error=sprint(showerror, err))
    end
end

end
