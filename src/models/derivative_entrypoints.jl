"""derivative_entrypoints.jl

Models 统一导数入口（Option B 过渡层）：
- 将历史 `PNJL.ThermoDerivatives` 能力通过 Models 暴露，避免调用方直接 include `src/pnjl/PNJL.jl`。
"""

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

@inline function _pnjl_thermo_derivatives_module()
    module_path = normpath(joinpath(@__DIR__, "derivatives", "ThermoDerivatives.jl"))
    td = IncludeOnce.include_once!(Main, :ModelsThermoDerivatives, module_path)
    return Base.invokelatest(td.derivatives_module_ref)
end

@inline function _invoke_td_worldsafe(symbol::Symbol, args...; kwargs...)
    mod = _pnjl_thermo_derivatives_module()
    fn = getproperty(mod, symbol)
    return fn(args...; kwargs...)
end

mass_derivatives(args...; kwargs...) = _invoke_td_worldsafe(:mass_derivatives, args...; kwargs...)
thermo_derivatives(args...; kwargs...) = _invoke_td_worldsafe(:thermo_derivatives, args...; kwargs...)
bulk_derivative_coeffs(args...; kwargs...) = _invoke_td_worldsafe(:bulk_derivative_coeffs, args...; kwargs...)
bulk_viscosity_coefficients(args...; kwargs...) = _invoke_td_worldsafe(:bulk_viscosity_coefficients, args...; kwargs...)
compute_B_bracket(args...; kwargs...) = _invoke_td_worldsafe(:compute_B_bracket, args...; kwargs...)
dP_dT(args...; kwargs...) = _invoke_td_worldsafe(:dP_dT, args...; kwargs...)
dP_dmu(args...; kwargs...) = _invoke_td_worldsafe(:dP_dmu, args...; kwargs...)