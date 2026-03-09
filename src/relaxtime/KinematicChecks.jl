module KinematicChecks

export KINEMATIC_ERROR_MODES, kinematic_error_mode, check_kinematic_threshold
export emit_regularization_notice, emit_negative_radicand_notice, emit_runtime_failure_notice

const KINEMATIC_ERROR_MODES = (:strict, :warn, :silent)
const EPS_THRESHOLD = 1e-12

@inline function kinematic_error_mode()::Symbol
    mode = Symbol(lowercase(strip(get(ENV, "RELAXTIME_KINEMATIC_ERROR_MODE", "warn"))))
    mode in KINEMATIC_ERROR_MODES || error(
        "Unknown RELAXTIME_KINEMATIC_ERROR_MODE=$(mode). Allowed: $(KINEMATIC_ERROR_MODES)"
    )
    return mode
end

function check_kinematic_threshold(
    s::Float64,
    m1::Float64,
    m2::Float64;
    warn_close::Bool=true,
    mode::Symbol=kinematic_error_mode(),
    context::AbstractString="Kinematic threshold violation",
)::Bool
    s_threshold = (m1 + m2)^2

    if s < s_threshold
        if mode === :strict
            throw(DomainError(s, "$(context): s = $(s) < threshold = $(s_threshold)"))
        elseif mode === :warn
            @warn context s=s threshold=s_threshold deficit=(s_threshold - s)
        end
        return false
    end

    s_plus = s - s_threshold
    if warn_close && s_plus < EPS_THRESHOLD && mode === :warn
        @warn "s is very close to threshold" context=context s=s threshold=s_threshold s_12_plus=s_plus maxlog=10
    end

    return true
end

function emit_regularization_notice(
    context::AbstractString,
    value::Float64,
    threshold::Float64;
    mode::Symbol=kinematic_error_mode(),
)
    message = "$(context): |value| = $(abs(value)) < $(threshold), applying regularization"
    if mode === :strict
        throw(DomainError(value, message))
    elseif mode === :warn
        @warn message value=value threshold=threshold maxlog=10
    end
    return nothing
end

function emit_negative_radicand_notice(
    context::AbstractString,
    delta::Float64,
    tolerance::Float64;
    mode::Symbol=kinematic_error_mode(),
)
    if delta > -tolerance
        message = "$(context): delta = $(delta) in (-$(tolerance), 0), setting momentum to 0"
    else
        message = "$(context): delta = $(delta) < -$(tolerance), setting momentum to 0"
    end

    if mode === :strict
        throw(DomainError(delta, message))
    elseif mode === :warn
        @warn message delta=delta tolerance=tolerance maxlog=10
    end
    return nothing
end

function emit_runtime_failure_notice(
    context::AbstractString,
    err;
    mode::Symbol=kinematic_error_mode(),
)
    if mode === :strict
        rethrow(err)
    elseif mode === :warn
        @warn context exception=err maxlog=10
    end
    return nothing
end

end # module KinematicChecks