"""factory.jl

模型工厂：根据 Symbol / Type 构造模型实例。

此处先满足 MVP：支持 `:NJL`。
"""

export create_model

const _LEGACY_PNJL_MODEL_PATH = joinpath(@__DIR__, "legacy", "LegacyPNJLModel.jl")
const _LEGACY_NJL_MODEL_PATH = joinpath(@__DIR__, "legacy", "LegacyNJLModel.jl")

@inline function _ensure_legacy_pnjl_loaded()
    if !isdefined(@__MODULE__, :LegacyPNJLModel)
        include(_LEGACY_PNJL_MODEL_PATH)
    end
    return nothing
end

@inline function _ensure_legacy_njl_loaded()
    if !isdefined(@__MODULE__, :LegacyNJLModel)
        include(_LEGACY_NJL_MODEL_PATH)
    end
    return nothing
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
    elseif kind === :RPNJL
        return RPNJLModel(; kwargs...)
    elseif kind === :LegacyNJL
        _ensure_legacy_njl_loaded()
        return Base.invokelatest(LegacyNJLModel; kwargs...)
    elseif kind === :LegacyPNJL
        _ensure_legacy_pnjl_loaded()
        return Base.invokelatest(LegacyPNJLModel)
    end
    error("Unknown model kind: ", kind)
end

"""create_model(::Type{T}; kwargs...) -> T

允许 `create_model(NJLModel; profile=...)` 这类调用。
"""
function create_model(::Type{T}; kwargs...) where {T<:AbstractQCDModel}
    return T(; kwargs...)
end
