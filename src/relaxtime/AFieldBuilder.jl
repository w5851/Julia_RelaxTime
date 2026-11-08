module AFieldBuilder

using ..GaussLegendre: gauleg, DEFAULT_COSΘ_NODES
using ..OneLoopIntegrals: A
using ..OneLoopIntegralsCorrection: A_aniso, A_aniso_reduced

export build_A_triplet, ensure_quark_params_has_A, A_BUILD_METHODS

const A_BUILD_METHODS = (:quadrature, :reduced)

const _AUTO_A_WARNED_ANISO = Ref(false)
const _AUTO_A_WARNED_ISO = Ref(false)

@inline _auto_a_warn_enabled() = get(ENV, "RELAXTIME_AUTO_A_WARN", "1") != "0"

@inline function _should_emit_auto_a_warn(use_aniso::Bool)
    flag = use_aniso ? _AUTO_A_WARNED_ANISO : _AUTO_A_WARNED_ISO
    flag[] && return false
    flag[] = true
    return true
end

@inline _xi(thermo_params::NamedTuple)::Float64 = hasproperty(thermo_params, :ξ) ? Float64(thermo_params.ξ) : 0.0

@inline function _normalize_a_build_method(method)::Symbol
    normalized = method isa Symbol ? method : Symbol(method)
    normalized in A_BUILD_METHODS ||
        throw(ArgumentError("unknown a_method=$(method); expected one of $(A_BUILD_METHODS)"))
    return normalized
end

function _build_A_value(
    m::Float64,
    μ::Float64,
    thermo_params::NamedTuple,
    nodes_p::AbstractVector{<:Real},
    weights_p::AbstractVector{<:Real},
    nodes_cos::AbstractVector{<:Real},
    weights_cos::AbstractVector{<:Real};
    a_method::Symbol,
    reduced_p_max::Union{Nothing,Real},
    reduced_radial_nodes::Int,
    reduced_scale::Float64,
)
    if a_method === :quadrature
        return A_aniso(m, μ, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
            _xi(thermo_params), collect(Float64, nodes_p), collect(Float64, weights_p),
            collect(Float64, nodes_cos), collect(Float64, weights_cos))
    end

    reduced_p_max === nothing && throw(ArgumentError(
        "a_method=:reduced requires an explicit reduced_p_max for a supplied momentum grid"))
    return A_aniso_reduced(m, μ, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
        _xi(thermo_params); p_max=Float64(reduced_p_max),
        radial_nodes=reduced_radial_nodes, semi_infinite_scale=reduced_scale)
end

@inline function _validate_inputs(quark_params::NamedTuple, thermo_params::NamedTuple)
    hasproperty(quark_params, :m) || error("quark_params is missing :m")
    hasproperty(quark_params, :μ) || error("quark_params is missing :μ")
    hasproperty(thermo_params, :T) || error("thermo_params is missing :T")
    hasproperty(thermo_params, :Φ) || error("thermo_params is missing :Φ")
    hasproperty(thermo_params, :Φbar) || error("thermo_params is missing :Φbar")
    return nothing
end

function build_A_triplet(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    nodes_p::AbstractVector,
    weights_p::AbstractVector;
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    use_aniso::Bool=true,
    a_method::Symbol=:quadrature,
    reduced_p_max::Union{Nothing,Real}=nothing,
    reduced_radial_nodes::Int=length(nodes_p),
    reduced_scale::Float64=5.0,
)::NamedTuple
    _validate_inputs(quark_params, thermo_params)

    length(nodes_p) == length(weights_p) || error("build_A_triplet: nodes_p and weights_p length mismatch")
    a_method = _normalize_a_build_method(a_method)
    reduced_radial_nodes > 0 || throw(ArgumentError("reduced_radial_nodes must be positive"))

    ξ = _xi(thermo_params)
    if use_aniso && abs(ξ) > 0.0
        nodes_cos, weights_cos = gauleg(-1.0, 1.0, cos_nodes)
        A_u = _build_A_value(quark_params.m.u, quark_params.μ.u, thermo_params,
            nodes_p, weights_p, nodes_cos, weights_cos;
            a_method, reduced_p_max, reduced_radial_nodes, reduced_scale)
        A_d = _build_A_value(quark_params.m.d, quark_params.μ.d, thermo_params,
            nodes_p, weights_p, nodes_cos, weights_cos;
            a_method, reduced_p_max, reduced_radial_nodes, reduced_scale)
        A_s = _build_A_value(quark_params.m.s, quark_params.μ.s, thermo_params,
            nodes_p, weights_p, nodes_cos, weights_cos;
            a_method, reduced_p_max, reduced_radial_nodes, reduced_scale)
        return (u=A_u, d=A_d, s=A_s)
    end

    p_nodes_f64 = collect(Float64, nodes_p)
    p_weights_f64 = collect(Float64, weights_p)
    A_u = A(quark_params.m.u, quark_params.μ.u, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, p_nodes_f64, p_weights_f64)
    A_d = A(quark_params.m.d, quark_params.μ.d, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, p_nodes_f64, p_weights_f64)
    A_s = A(quark_params.m.s, quark_params.μ.s, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, p_nodes_f64, p_weights_f64)
    return (u=A_u, d=A_d, s=A_s)
end

function build_A_triplet(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    p_nodes::Int=16,
    p_max::Float64=20.0,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    use_aniso::Bool=true,
    a_method::Symbol=:quadrature,
    reduced_p_max::Union{Nothing,Real}=nothing,
    reduced_radial_nodes::Int=p_nodes,
    reduced_scale::Float64=5.0,
)::NamedTuple
    a_method = _normalize_a_build_method(a_method)
    if a_method === :reduced && isinf(p_max)
        _validate_inputs(quark_params, thermo_params)
        ξ = _xi(thermo_params)
        if use_aniso && abs(ξ) > 0.0
            limit = reduced_p_max === nothing ? Inf : Float64(reduced_p_max)
            return (
                u=A_aniso_reduced(Float64(quark_params.m.u), Float64(quark_params.μ.u),
                    Float64(thermo_params.T), Float64(thermo_params.Φ), Float64(thermo_params.Φbar), ξ;
                    p_max=limit, radial_nodes=reduced_radial_nodes, semi_infinite_scale=reduced_scale),
                d=A_aniso_reduced(Float64(quark_params.m.d), Float64(quark_params.μ.d),
                    Float64(thermo_params.T), Float64(thermo_params.Φ), Float64(thermo_params.Φbar), ξ;
                    p_max=limit, radial_nodes=reduced_radial_nodes, semi_infinite_scale=reduced_scale),
                s=A_aniso_reduced(Float64(quark_params.m.s), Float64(quark_params.μ.s),
                    Float64(thermo_params.T), Float64(thermo_params.Φ), Float64(thermo_params.Φbar), ξ;
                    p_max=limit, radial_nodes=reduced_radial_nodes, semi_infinite_scale=reduced_scale),
            )
        end
        p_max = 20.0
    end
    isfinite(p_max) && p_max > 0.0 || throw(ArgumentError("p_max must be finite and positive for this A builder"))
    nodes_p, weights_p = gauleg(0.0, p_max, p_nodes)
    return build_A_triplet(
        quark_params,
        thermo_params,
        nodes_p,
        weights_p;
        cos_nodes=cos_nodes,
        use_aniso=use_aniso,
        a_method=a_method,
        reduced_p_max=a_method === :reduced && reduced_p_max === nothing ? p_max : reduced_p_max,
        reduced_radial_nodes=reduced_radial_nodes,
        reduced_scale=reduced_scale,
    )
end

function ensure_quark_params_has_A(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    p_nodes::Int=16,
    p_max::Float64=20.0,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    use_aniso::Bool=true,
    warn_on_auto::Bool=true,
    a_method::Symbol=:quadrature,
    reduced_p_max::Union{Nothing,Real}=nothing,
    reduced_radial_nodes::Int=p_nodes,
    reduced_scale::Float64=5.0,
)::NamedTuple
    if hasproperty(quark_params, :A)
        return quark_params
    end

    ξ = _xi(thermo_params)
    if warn_on_auto && abs(ξ) > 0.0 && _auto_a_warn_enabled() && _should_emit_auto_a_warn(use_aniso)
        if use_aniso
            @warn "quark_params missing :A under ξ≠0; auto-building anisotropic A via A_aniso" ξ=ξ p_nodes=p_nodes p_max=p_max cos_nodes=cos_nodes
        else
            @warn "quark_params missing :A under ξ≠0; use_aniso=false so isotropic A will be used" ξ=ξ p_nodes=p_nodes p_max=p_max
        end
    end

    A_vals = build_A_triplet(
        quark_params,
        thermo_params;
        p_nodes=p_nodes,
        p_max=p_max,
        cos_nodes=cos_nodes,
        use_aniso=use_aniso,
        a_method=a_method,
        reduced_p_max=reduced_p_max,
        reduced_radial_nodes=reduced_radial_nodes,
        reduced_scale=reduced_scale,
    )
    return merge(quark_params, (A=A_vals,))
end

end # module AFieldBuilder
