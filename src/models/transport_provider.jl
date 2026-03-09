"""transport_provider

为输运层提供一个“models-side 可注入 provider”（分布函数/色散关系）。

设计目标（阶段 4）：
- Transport 层只依赖 provider 的函数集合，而不直接依赖 PNJL 旧实现细节。
- 目前最小切片：仍复用现有 PNJL 分布函数实现；后续可替换为 models 侧实现。

返回值：NamedTuple，字段与 `TransportCoefficients` 的 provider 约定一致：
- energy_from_p(p, m)
- quark_distribution(E, μ, T, Φ, Φbar)
- antiquark_distribution(E, μ, T, Φ, Φbar)
- quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
- antiquark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)

可选扩展（阶段 4 逐步启用）：
- mass_for_species(species, quark_params, thermo_params)
- mu_for_species(species, quark_params, thermo_params)
"""

@inline energy_from_p(p::Float64, m::Float64) = sqrt(p * p + m * m)

@inline energy_from_p_aniso(p::Float64, m::Float64, ξ::Float64, cosθ::Float64) = sqrt(p * p + m * m + ξ * (p * cosθ)^2)

"""TransportProvider

对象型 provider：字段名与 TransportCoefficients 约定一致，允许携带额外上下文（ctx）。

约定：Transport 层只依赖字段函数（`.energy_from_p` 等），不直接依赖 ctx。
"""
struct TransportProvider{EF,EFA,QD,AQD,QDA,AQDA,PEA,CTX}
    energy_from_p::EF
    energy_from_p_aniso::EFA
    quark_distribution::QD
    antiquark_distribution::AQD
    quark_distribution_aniso::QDA
    antiquark_distribution_aniso::AQDA
    prefer_energy_aniso::PEA
    mass_for_species
    mu_for_species
    ctx::CTX
end

@inline function _mass_from_quark_params(species::Symbol, quark_params)::Float64
    if species in (:u, :d)
        return species === :u ? Float64(quark_params.m.u) : Float64(quark_params.m.d)
    elseif species === :s
        return Float64(quark_params.m.s)
    elseif species in (:ubar, :dbar)
        return species === :ubar ? Float64(quark_params.m.u) : Float64(quark_params.m.d)
    elseif species === :sbar
        return Float64(quark_params.m.s)
    end
    throw(ArgumentError("unknown species=$species"))
end

@inline function _mu_from_quark_params(species::Symbol, quark_params)::Float64
    if species in (:u, :ubar)
        return Float64(quark_params.μ.u)
    elseif species in (:d, :dbar)
        return Float64(quark_params.μ.d)
    elseif species in (:s, :sbar)
        return Float64(quark_params.μ.s)
    end
    throw(ArgumentError("unknown species=$species"))
end

@inline function _mass_from_ctx_or_params(ctx, species::Symbol, quark_params, thermo_params)::Float64
    masses = get(ctx, :masses, nothing)
    if masses !== nothing
        if species in (:u, :ubar)
            return Float64(masses.u)
        elseif species in (:d, :dbar)
            return Float64(masses.d)
        elseif species in (:s, :sbar)
            return Float64(masses.s)
        end
    end
    return _mass_from_quark_params(species, quark_params)
end

@inline function _mu_from_ctx_or_params(ctx, species::Symbol, quark_params, thermo_params)::Float64
    mu_cache = get(ctx, :mu, nothing)
    if mu_cache !== nothing
        if species in (:u, :ubar)
            return Float64(mu_cache.u)
        elseif species in (:d, :dbar)
            return Float64(mu_cache.d)
        elseif species in (:s, :sbar)
            return Float64(mu_cache.s)
        end
    end
    return _mu_from_quark_params(species, quark_params)
end

@inline function transport_provider(model::AbstractQCDModel)
    if model isa PNJLModel
        ctx0 = (backend=:models, impl=:models_pnjl_distributions, model=model)
        mass_for_species = (species, quark_params, thermo_params) -> _mass_from_ctx_or_params(ctx0, species, quark_params, thermo_params)
        mu_for_species = (species, quark_params, thermo_params) -> _mu_from_ctx_or_params(ctx0, species, quark_params, thermo_params)

        return TransportProvider(
            energy_from_p,
            energy_from_p_aniso,
            pnjl_quark_distribution,
            pnjl_antiquark_distribution,
            pnjl_quark_distribution_aniso,
            pnjl_antiquark_distribution_aniso,
            true,
            mass_for_species,
            mu_for_species,
            ctx0,
        )
    end
    throw(ArgumentError("transport_provider(model) is not implemented for model type $(typeof(model))"))
end

"""prepare_transport_provider(provider, equilibrium; quark_params, thermo_params, masses=nothing)

用平衡态输出（或显式 masses）为 provider 填充/缓存 ctx，使 provider.ctx 真正可用于色散/有效质量等依赖。

当前最小落地：缓存 (u,d,s) 有效质量与 (u,d,s) 化学势，并通过 provider.mass_for_species / provider.mu_for_species 生效。
"""
function prepare_transport_provider(
    provider::TransportProvider,
    equilibrium;
    quark_params,
    thermo_params,
    masses=nothing,
)
    masses_val = masses
    if masses_val === nothing
        if equilibrium !== nothing && hasproperty(equilibrium, :masses)
            masses_val = getproperty(equilibrium, :masses)
        end
    end

    masses_nt = if masses_val === nothing
        (u=Float64(quark_params.m.u), d=Float64(quark_params.m.d), s=Float64(quark_params.m.s))
    elseif masses_val isa NamedTuple
        (u=Float64(masses_val.u), d=Float64(masses_val.d), s=Float64(masses_val.s))
    else
        (u=Float64(masses_val[1]), d=Float64(masses_val[2]), s=Float64(masses_val[3]))
    end
    mu_nt = (u=Float64(quark_params.μ.u), d=Float64(quark_params.μ.d), s=Float64(quark_params.μ.s))

    ctx1 = merge(provider.ctx, (masses=masses_nt, mu=mu_nt))
    mass_for_species = (species, qp, tp) -> _mass_from_ctx_or_params(ctx1, species, qp, tp)
    mu_for_species = (species, qp, tp) -> _mu_from_ctx_or_params(ctx1, species, qp, tp)

    return TransportProvider(
        provider.energy_from_p,
        provider.energy_from_p_aniso,
        provider.quark_distribution,
        provider.antiquark_distribution,
        provider.quark_distribution_aniso,
        provider.antiquark_distribution_aniso,
        provider.prefer_energy_aniso,
        mass_for_species,
        mu_for_species,
        ctx1,
    )
end

function prepare_transport_provider(
    provider::NamedTuple,
    equilibrium;
    quark_params,
    thermo_params,
    masses=nothing,
)
    masses_val = masses
    if masses_val === nothing
        if equilibrium !== nothing && hasproperty(equilibrium, :masses)
            masses_val = getproperty(equilibrium, :masses)
        end
    end

    masses_nt = if masses_val === nothing
        (u=Float64(quark_params.m.u), d=Float64(quark_params.m.d), s=Float64(quark_params.m.s))
    elseif masses_val isa NamedTuple
        (u=Float64(masses_val.u), d=Float64(masses_val.d), s=Float64(masses_val.s))
    else
        (u=Float64(masses_val[1]), d=Float64(masses_val[2]), s=Float64(masses_val[3]))
    end
    mu_nt = (u=Float64(quark_params.μ.u), d=Float64(quark_params.μ.d), s=Float64(quark_params.μ.s))
    ctx1 = (masses=masses_nt, mu=mu_nt)

    out = merge(provider, (ctx=ctx1,))
    if !hasproperty(out, :energy_from_p_aniso)
        out = merge(out, (energy_from_p_aniso=energy_from_p_aniso,))
    end
    if !hasproperty(out, :prefer_energy_aniso)
        out = merge(out, (prefer_energy_aniso=true,))
    end
    if !hasproperty(out, :mass_for_species)
        mass_for_species = (species, qp, tp) -> _mass_from_ctx_or_params(ctx1, species, qp, tp)
        out = merge(out, (mass_for_species=mass_for_species,))
    end
    if !hasproperty(out, :mu_for_species)
        mu_for_species = (species, qp, tp) -> _mu_from_ctx_or_params(ctx1, species, qp, tp)
        out = merge(out, (mu_for_species=mu_for_species,))
    end

    return out
end
