"""derivative_entrypoints.jl

Models 统一导数入口（Option B 过渡层）：
- 将历史 `PNJL.ThermoDerivatives` 能力通过 Models 暴露，避免调用方直接 include `src/pnjl/PNJL.jl`。
"""

@inline function _pnjl_module_for_derivatives()
    include_once_path = normpath(joinpath(@__DIR__, "..", "utils", "IncludeOnce.jl"))
    if !isdefined(Main, :IncludeOnce)
        Base.include(Main, include_once_path)
    end
    include_once = Main.IncludeOnce
    pnjl_entry_path = normpath(joinpath(@__DIR__, "..", "pnjl", "PNJL.jl"))
    return include_once.include_once!(Main, :PNJL, pnjl_entry_path)
end

@inline function _pnjl_thermo_derivatives_module()
    pnjl = _pnjl_module_for_derivatives()
    isdefined(pnjl, :ThermoDerivatives) || error("PNJL.ThermoDerivatives is not available")
    return getproperty(pnjl, :ThermoDerivatives)
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