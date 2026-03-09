module ParameterAdapters

const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "types", "ParameterTypes.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PARAMETER_TYPES_PATH)
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

export normalize_quark_input, normalize_thermo_input

@inline normalize_quark_input(q::QuarkParams) = as_namedtuple(q)
@inline normalize_quark_input(q::NamedTuple) = q

@inline normalize_thermo_input(t::ThermoParams) = as_namedtuple(t)
@inline normalize_thermo_input(t::NamedTuple) = t

end # module ParameterAdapters