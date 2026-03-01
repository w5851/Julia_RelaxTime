"""factory.jl

模型工厂：根据 Symbol / Type 构造模型实例。

此处先满足 MVP：支持 `:NJL`。
"""

export create_model

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _COMPAT_PNJL_MODEL_PATH = normpath(joinpath(@__DIR__, "legacy", "LegacyPNJLModel.jl"))
const _LEGACY_NJL_MODEL_PATH = normpath(joinpath(@__DIR__, "legacy", "LegacyNJLModel.jl"))

@inline function _compat_pnjl_ctor()
    if !isdefined(@__MODULE__, :LegacyPNJLModel)
        Base.include(@__MODULE__, _COMPAT_PNJL_MODEL_PATH)
    end
    return getfield(@__MODULE__, :LegacyPNJLModel)
end

@inline function _legacy_njl_ctor()
    if !isdefined(@__MODULE__, :LegacyNJLModel)
        Base.include(@__MODULE__, _LEGACY_NJL_MODEL_PATH)
    end
    return getfield(@__MODULE__, :LegacyNJLModel)
end

@inline function _construct_worldsafe(ctor; kwargs...)
    return Base.invokelatest(ctor; kwargs...)
end

"""create_model(kind; kwargs...) -> model

示例：
- `create_model(:NJL; profile="default")`
"""
function create_model(kind::Symbol; kwargs...)
    if kind === :NJL
        return NJLModel(; kwargs...)
    elseif kind === :NJL2
        return NJL2Model(; kwargs...)
    elseif kind === :PNJL
        return PNJLModel(; kwargs...)
    elseif kind === :PNJLMagnetic
        return PNJLMagneticModel(; kwargs...)
    elseif kind === :RPNJL
        return RPNJLModel(; kwargs...)
    elseif kind === :LegacyNJL
        ctor = _legacy_njl_ctor()
        return _construct_worldsafe(ctor; kwargs...)
    elseif kind === :LegacyPNJL
        ctor = _compat_pnjl_ctor()
        return _construct_worldsafe(ctor)
    end
    error("Unknown model kind: ", kind)
end

"""create_model(::Type{T}; kwargs...) -> T

允许 `create_model(NJLModel; profile=...)` 这类调用。
"""
function create_model(::Type{T}; kwargs...) where {T<:AbstractQCDModel}
    return T(; kwargs...)
end
