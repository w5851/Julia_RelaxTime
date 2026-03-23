"""factory.jl

模型工厂：根据 Symbol / Type 构造模型实例。

此处先满足 MVP：支持 `:NJL`。
"""

export create_model

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
