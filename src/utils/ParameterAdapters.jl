module ParameterAdapters

const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "types", "ParameterTypes.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PARAMETER_TYPES_PATH)
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

export normalize_quark_input, normalize_thermo_input
export normalize_symbol_mapping_input, lookup_symbol_value

@inline normalize_quark_input(q::QuarkParams) = as_namedtuple(q)
@inline normalize_quark_input(q::NamedTuple) = q

@inline normalize_thermo_input(t::ThermoParams) = as_namedtuple(t)
@inline normalize_thermo_input(t::NamedTuple) = t

@inline normalize_symbol_mapping_input(values::NamedTuple, context::AbstractString) = values

function normalize_symbol_mapping_input(values::AbstractDict, context::AbstractString)
    return (; (Symbol(k) => v for (k, v) in pairs(values))...)
end

function normalize_symbol_mapping_input(values, context::AbstractString)
    throw(ArgumentError("$context must be a NamedTuple or Dict-like mapping"))
end

@inline function lookup_symbol_value(values::NamedTuple, key::Symbol, context::AbstractString)
    hasproperty(values, key) || throw(ArgumentError("$context is missing $(key)"))
    return getproperty(values, key)
end

end # module ParameterAdapters