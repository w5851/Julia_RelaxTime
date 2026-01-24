"""
    ParameterTypes

项目级可复用的“参数结构体”模块。

设计目标：
- 给多个顶层模块（如 `PNJL`/`RelaxationTime`/`TransportCoefficients` 等）提供共享的参数类型；
- 这些类型挂在 `Main.ParameterTypes` 下，其他模块可用 `using ..ParameterTypes` 直接引用，避免在子模块里重复定义导致类型不一致。

当前包含：
- `QuarkParams`：封装质量 `m=(u,d,s)` 与化学势 `μ=(u,d,s)`
- `ThermoParams`：封装 `T, Φ, Φbar, ξ`（`ξ` 默认 0）
- `as_namedtuple`：将结构体转换为旧接口所需的 `NamedTuple`
"""

module ParameterTypes

export QuarkParams, ThermoParams, as_namedtuple

struct QuarkParams
    m::NamedTuple
    μ::NamedTuple
end

function QuarkParams(quark_params::NamedTuple)
    hasproperty(quark_params, :m) || error("QuarkParams: input is missing field :m")
    hasproperty(quark_params, :μ) || error("QuarkParams: input is missing field :μ")
    return QuarkParams(quark_params.m, quark_params.μ)
end

@inline as_namedtuple(q::QuarkParams) = (m=q.m, μ=q.μ)

struct ThermoParams
    T::Float64
    Φ::Float64
    Φbar::Float64
    ξ::Float64
end

function ThermoParams(thermo_params::NamedTuple)
    hasproperty(thermo_params, :T) || error("ThermoParams: input is missing field :T")
    hasproperty(thermo_params, :Φ) || error("ThermoParams: input is missing field :Φ")
    hasproperty(thermo_params, :Φbar) || error("ThermoParams: input is missing field :Φbar")
    ξ = get(thermo_params, :ξ, 0.0)
    return ThermoParams(Float64(thermo_params.T), Float64(thermo_params.Φ), Float64(thermo_params.Φbar), Float64(ξ))
end

@inline as_namedtuple(t::ThermoParams) = (T=t.T, Φ=t.Φ, Φbar=t.Φbar, ξ=t.ξ)

end # module
