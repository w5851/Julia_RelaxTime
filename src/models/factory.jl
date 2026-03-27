"""factory.jl

模型工厂：根据 Symbol / Type 构造模型实例。

此处先满足 MVP：支持 `:NJL`。
"""

using Base.Threads: ReentrantLock, lock, unlock

export create_model
export get_cached_model, clear_model_cache!
export register_model!, unregister_model!, registered_model_kinds

const _MODEL_CACHE = Dict{Symbol, AbstractQCDModel}()
const _MODEL_CACHE_LOCK = ReentrantLock()
const _MODEL_REGISTRY = Dict{Symbol, Function}()
const _MODEL_REGISTRY_LOCK = ReentrantLock()

@inline function _set_model_builder!(kind::Symbol, builder::Function)
    lock(_MODEL_REGISTRY_LOCK)
    try
        _MODEL_REGISTRY[kind] = builder
    finally
        unlock(_MODEL_REGISTRY_LOCK)
    end
    return nothing
end

function _register_builtin_models!()
    _set_model_builder!(:NJL, (; kwargs...) -> NJLModel(; kwargs...))
    _set_model_builder!(:NJL2, (; kwargs...) -> NJL2Model(; kwargs...))
    _set_model_builder!(:PNJL, (; kwargs...) -> PNJLModel(; kwargs...))
    _set_model_builder!(:PNJLMagnetic, (; kwargs...) -> PNJLMagneticModel(; kwargs...))
    _set_model_builder!(:RPNJL, (; kwargs...) -> RPNJLModel(; kwargs...))
    _set_model_builder!(:Rotation, (; kwargs...) -> RotationModel(; kwargs...))
    _set_model_builder!(:GasLiquid, (; kwargs...) -> GasLiquidModel(; kwargs...))
    return nothing
end

"""create_model(kind; kwargs...) -> model

示例：
- `create_model(:NJL; profile="default")`
"""
function create_model(kind::Symbol; kwargs...)
    lock(_MODEL_REGISTRY_LOCK)
    local builder
    try
        builder = get(_MODEL_REGISTRY, kind, nothing)
    finally
        unlock(_MODEL_REGISTRY_LOCK)
    end

    if isnothing(builder)
        error("Unknown model kind: ", kind)
    end

    return builder(; kwargs...)
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

function register_model!(kind::Symbol, builder::Function)
    _set_model_builder!(kind, builder)
    clear_model_cache!(kinds=(kind,))
    return nothing
end

function unregister_model!(kind::Symbol)
    lock(_MODEL_REGISTRY_LOCK)
    try
        delete!(_MODEL_REGISTRY, kind)
    finally
        unlock(_MODEL_REGISTRY_LOCK)
    end
    clear_model_cache!(kinds=(kind,))
    return nothing
end

function registered_model_kinds()
    lock(_MODEL_REGISTRY_LOCK)
    try
        return sort!(collect(keys(_MODEL_REGISTRY)))
    finally
        unlock(_MODEL_REGISTRY_LOCK)
    end
end

_register_builtin_models!()
