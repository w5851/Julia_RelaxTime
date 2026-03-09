module WorkflowParamAdapters

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "..", "types", "ParameterTypes.jl"))
IncludeOnce.include_once!(Main, :ParameterTypes, _PARAMETER_TYPES_PATH)
using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

export normalize_quark_params, normalize_thermo_params, as_legacy_inputs

@inline function _ensure_real_field(scope::AbstractString, name::Symbol, value)
    value isa Real || throw(ArgumentError("$(scope).$(name) must be Real, got $(typeof(value))"))
    isfinite(Float64(value)) || throw(ArgumentError("$(scope).$(name) must be finite, got $(value)"))
    return nothing
end

@inline function _get_required_field(nt::NamedTuple, scope::AbstractString, name::Symbol)
    hasproperty(nt, name) || throw(ArgumentError("$(scope) must contain field :$(name)"))
    return getproperty(nt, name)
end

@inline function _normalize_flavor_triplet(value, scope::AbstractString, name::Symbol)
    value isa NamedTuple || throw(ArgumentError("$(scope).$(name) must be NamedTuple with keys (:u,:d,:s), got $(typeof(value))"))
    for f in (:u, :d, :s)
        hasproperty(value, f) || throw(ArgumentError("$(scope).$(name) must contain field :$(f)"))
        _ensure_real_field("$(scope).$(name)", f, getproperty(value, f))
    end
    return (
        u=Float64(getproperty(value, :u)),
        d=Float64(getproperty(value, :d)),
        s=Float64(getproperty(value, :s)),
    )
end

@inline function _normalize_quark_namedtuple(q::NamedTuple)
    m = _normalize_flavor_triplet(_get_required_field(q, "quark_params", :m), "quark_params", :m)
    mu = _normalize_flavor_triplet(_get_required_field(q, "quark_params", :μ), "quark_params", :μ)
    return (m=m, μ=mu)
end

@inline function _normalize_thermo_namedtuple(t::NamedTuple)
    T = _get_required_field(t, "thermo_params", :T)
    Phi = _get_required_field(t, "thermo_params", :Φ)
    Phibar = _get_required_field(t, "thermo_params", :Φbar)
    xi = _get_required_field(t, "thermo_params", :ξ)

    _ensure_real_field("thermo_params", :T, T)
    _ensure_real_field("thermo_params", :Φ, Phi)
    _ensure_real_field("thermo_params", :Φbar, Phibar)
    _ensure_real_field("thermo_params", :ξ, xi)

    return (
        T=Float64(T),
        Φ=Float64(Phi),
        Φbar=Float64(Phibar),
        ξ=Float64(xi),
    )
end

@inline function normalize_quark_params(q::QuarkParams)
    return _normalize_quark_namedtuple(as_namedtuple(q))
end

@inline function normalize_quark_params(q)
    throw(ArgumentError("quark_params must be QuarkParams, got $(typeof(q))"))
end

@inline function normalize_thermo_params(t::ThermoParams)
    return _normalize_thermo_namedtuple(as_namedtuple(t))
end

@inline function normalize_thermo_params(t)
    throw(ArgumentError("thermo_params must be ThermoParams, got $(typeof(t))"))
end

@inline as_legacy_inputs(q, t) = (quark_params=normalize_quark_params(q), thermo_params=normalize_thermo_params(t))

end # module