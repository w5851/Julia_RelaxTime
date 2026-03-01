if !isdefined(Main, :ModelsThermoDerivatives)
    @eval module ModelsThermoDerivatives

    const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
    if !isdefined(Main, :IncludeOnce)
        Base.include(Main, _INCLUDE_ONCE_PATH)
    end
    const IncludeOnce = Main.IncludeOnce

    const _PNJL_ENTRY_PATH = normpath(joinpath(@__DIR__, "..", "pnjl", "PNJL.jl"))

    function derivatives_module_ref()
        pnjl = IncludeOnce.include_once!(Main, :PNJL, _PNJL_ENTRY_PATH)
        isdefined(pnjl, :ThermoDerivatives) || error("PNJL.ThermoDerivatives is not available")
        return getproperty(pnjl, :ThermoDerivatives)
    end

    end
end
