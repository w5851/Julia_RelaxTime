using NLsolve

"""
    RootProblemSpec

通用根求解问题描述。
"""
struct RootProblemSpec{F,C}
    residual!::F
    ctx::C
    x_dim::Int
    branch_kind::Symbol
end

"""
    RootPolicy

通用根求解策略配置。
"""
struct RootPolicy
    primary_method::Symbol
    fallback_method::Symbol
    use_fallback::Bool
    use_multiseed::Bool
    residual_norm_max::Float64
    require_converged::Bool
    diagnostics_level::Symbol
end

function RootPolicy(;
    primary_method::Symbol=:newton,
    fallback_method::Symbol=:trust_region,
    use_fallback::Bool=true,
    use_multiseed::Bool=false,
    residual_norm_max::Real=1e-6,
    require_converged::Bool=true,
    diagnostics_level::Symbol=:basic,
)
    return RootPolicy(
        primary_method,
        fallback_method,
        use_fallback,
        use_multiseed,
        Float64(residual_norm_max),
        require_converged,
        diagnostics_level,
    )
end

mutable struct ContinuationState
    seed_by_branch::Dict{Symbol,Vector{Float64}}
    last_point::Union{Nothing,NamedTuple}
end

ContinuationState() = ContinuationState(Dict{Symbol,Vector{Float64}}(), nothing)

struct RootAttempt
    method::Symbol
    seed_source::Symbol
    converged::Bool
    residual_norm::Float64
    quality_tag::Symbol
end

struct RootDiagnostics
    selected_attempt::Int
    attempts::Vector{RootAttempt}
    branch_key::Symbol
    notes::Vector{String}
end

struct RootSolveResult
    x::Vector{Float64}
    converged::Bool
    residual_norm::Float64
    quality_tag::Symbol
    diagnostics::RootDiagnostics
end

@inline function _normalize_branch_key(branch_key, x, ctx)::Symbol
    if branch_key === nothing
        return :default
    elseif branch_key isa Symbol
        return branch_key
    elseif branch_key isa Function
        key = branch_key(x, ctx)
        key isa Symbol || throw(ArgumentError("branch_key(x, ctx) must return Symbol, got $(typeof(key))"))
        return key
    else
        throw(ArgumentError("branch_key must be nothing, Symbol, or Function"))
    end
end

@inline function _domain_quality_ok(domain_quality, x, ctx)
    domain_quality === nothing && return true
    out = domain_quality(x, ctx)
    if out isa NamedTuple
        return Bool(get(out, :ok, false))
    elseif out isa Bool
        return out
    end
    throw(ArgumentError("domain_quality(x, ctx) must return Bool or NamedTuple with :ok"))
end

@inline function _quality_tag(converged::Bool, residual_norm::Float64, domain_ok::Bool, policy::RootPolicy; is_fallback::Bool=false)
    residual_ok = isfinite(residual_norm) && residual_norm <= policy.residual_norm_max
    converged_ok = !policy.require_converged || converged
    good = converged_ok && residual_ok && domain_ok
    if good
        return is_fallback ? :fallback : :good
    end
    if converged || residual_ok
        return :degraded
    end
    return :bad
end

function _solve_once(spec::RootProblemSpec, seed::Vector{Float64}, method::Symbol)
    residual_fn! = (F, x) -> spec.residual!(F, x, spec.ctx)
    res = NLsolve.nlsolve(
        residual_fn!,
        seed;
        autodiff=:forward,
        method=method,
        xtol=1e-9,
        ftol=1e-9,
    )
    return res
end

function _candidate_from_result(res)
    return Vector{Float64}(res.zero), Bool(res.f_converged), Float64(res.residual_norm)
end

function _multiseed_candidates(seed::Vector{Float64})
    out = Vector{Vector{Float64}}()
    push!(out, copy(seed))
    if !isempty(seed)
        step = 1e-2
        push!(out, [seed[i] + (i == 1 ? step : 0.0) for i in eachindex(seed)])
        push!(out, [seed[i] - (i == 1 ? step : 0.0) for i in eachindex(seed)])
    end
    return out
end

function solve_root_with_policy(
    spec::RootProblemSpec,
    x0::AbstractVector{<:Real};
    policy::RootPolicy=RootPolicy(),
    continuation_state::Union{Nothing,ContinuationState}=nothing,
    domain_quality=nothing,
    branch_key=nothing,
)
    length(x0) == spec.x_dim || throw(ArgumentError("x0 length $(length(x0)) does not match spec.x_dim $(spec.x_dim)"))

    seed = Float64.(x0)
    branch = _normalize_branch_key(branch_key, seed, spec.ctx)
    seed_source = :default
    if continuation_state !== nothing && haskey(continuation_state.seed_by_branch, branch)
        cseed = continuation_state.seed_by_branch[branch]
        if length(cseed) == spec.x_dim
            seed = copy(cseed)
            seed_source = :continuation
        end
    end

    attempts = RootAttempt[]
    notes = String[]
    selected_idx = 0
    selected_x = copy(seed)
    selected_conv = false
    selected_resn = Inf
    selected_tag = :bad

    function register_attempt(method::Symbol, source::Symbol, x::Vector{Float64}, conv::Bool, resn::Float64; is_fallback::Bool=false)
        domain_ok = _domain_quality_ok(domain_quality, x, spec.ctx)
        qtag = _quality_tag(conv, resn, domain_ok, policy; is_fallback=is_fallback)
        push!(attempts, RootAttempt(method, source, conv, resn, qtag))
        return qtag
    end

    res_primary = _solve_once(spec, seed, policy.primary_method)
    x_primary, conv_primary, resn_primary = _candidate_from_result(res_primary)
    tag_primary = register_attempt(policy.primary_method, seed_source, x_primary, conv_primary, resn_primary)
    selected_idx = 1
    selected_x = x_primary
    selected_conv = conv_primary
    selected_resn = resn_primary
    selected_tag = tag_primary

    if policy.use_fallback && tag_primary in (:degraded, :bad)
        res_fb = _solve_once(spec, seed, policy.fallback_method)
        x_fb, conv_fb, resn_fb = _candidate_from_result(res_fb)
        tag_fb = register_attempt(policy.fallback_method, seed_source, x_fb, conv_fb, resn_fb; is_fallback=true)
        if tag_fb in (:good, :fallback) || (tag_primary in (:degraded, :bad) && resn_fb < selected_resn)
            selected_idx = length(attempts)
            selected_x = x_fb
            selected_conv = conv_fb
            selected_resn = resn_fb
            selected_tag = tag_fb
        end
    end

    if policy.use_multiseed && selected_tag in (:degraded, :bad)
        for s in _multiseed_candidates(seed)
            res_ms = _solve_once(spec, s, policy.primary_method)
            x_ms, conv_ms, resn_ms = _candidate_from_result(res_ms)
            tag_ms = register_attempt(policy.primary_method, :multiseed, x_ms, conv_ms, resn_ms)
            if tag_ms in (:good, :fallback) || resn_ms < selected_resn
                selected_idx = length(attempts)
                selected_x = x_ms
                selected_conv = conv_ms
                selected_resn = resn_ms
                selected_tag = tag_ms
            end
        end
        push!(notes, "multiseed attempted")
    end

    if continuation_state !== nothing
        continuation_state.seed_by_branch[branch] = copy(selected_x)
    end

    diag = RootDiagnostics(selected_idx, attempts, branch, notes)
    return RootSolveResult(selected_x, selected_conv, selected_resn, selected_tag, diag)
end

function solve_root_continuation(
    scan_points,
    spec_factory::Function;
    policy::RootPolicy=RootPolicy(),
    tracker::Union{Nothing,ContinuationState}=ContinuationState(),
    x0::AbstractVector{<:Real},
    domain_quality=nothing,
    branch_key=nothing,
)
    results = RootSolveResult[]
    current_seed = Float64.(x0)

    for point in scan_points
        spec = spec_factory(point)
        res = solve_root_with_policy(
            spec,
            current_seed;
            policy=policy,
            continuation_state=tracker,
            domain_quality=domain_quality,
            branch_key=branch_key,
        )
        push!(results, res)
        current_seed = copy(res.x)
        if tracker !== nothing
            tracker.last_point = (point=point,)
        end
    end

    return results
end
