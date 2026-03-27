"""factory.jl

模型工厂：根据 Symbol / Type 构造模型实例。

此处先满足 MVP：支持 `:NJL`。
"""

using Base.Threads: ReentrantLock, lock, unlock

export create_model
export get_cached_model, clear_model_cache!

const _MODEL_CACHE = Dict{Symbol, AbstractQCDModel}()
const _MODEL_CACHE_LOCK = ReentrantLock()

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
    elseif kind === :Rotation
        return RotationModel(; kwargs...)
    elseif kind === :GasLiquid
        return GasLiquidModel(; kwargs...)
    end
    error("Unknown model kind: ", kind)
end

"""create_model(::Type{T}; kwargs...) -> T

允许 `create_model(NJLModel; profile=...)` 这类调用。
"""
function create_model(::Type{T}; kwargs...) where {T<:AbstractQCDModel}
    return T(; kwargs...)
end

@inline function get_cached_model(kind::Symbol)
    lock(_MODEL_CACHE_LOCK)
    try
        return get!(_MODEL_CACHE, kind) do
            create_model(kind)
        end
    finally
        unlock(_MODEL_CACHE_LOCK)
    end
end

function clear_model_cache!(; kinds=nothing)
    lock(_MODEL_CACHE_LOCK)
    try
        if isnothing(kinds)
            empty!(_MODEL_CACHE)
            return nothing
        end
        for kind in kinds
            delete!(_MODEL_CACHE, kind)
        end
    finally
        unlock(_MODEL_CACHE_LOCK)
    end
    return nothing
end
