if !isdefined(Main, :ModelsScanEntrypoints)
    @eval module ModelsScanEntrypoints

    const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
    if !isdefined(Main, :IncludeOnce)
        Base.include(Main, _INCLUDE_ONCE_PATH)
    end
    const IncludeOnce = Main.IncludeOnce

    const _PNJL_ENTRY_PATH = normpath(joinpath(@__DIR__, "..", "pnjl", "PNJL.jl"))

    function pnjl_module_ref()
        return IncludeOnce.include_once!(Main, :PNJL, _PNJL_ENTRY_PATH)
    end

    end
end
