module TransportWorkflow

"""
将"各向异性 PNJL 能隙方程求解"与"RTA 输运系数计算"串联起来，并可选在内部计算弛豫时间 τ。

工作流：
1) 调用 `PNJL.solve(FixedMu(), T_fm, μ_fm)` 在给定 `(T, μ, ξ)` 下求平衡（能隙/Polyakov 环）。
2) 从平衡解构造 `quark_params` 与 `thermo_params`。
3) （可选）计算数密度与平均散射率，从而得到 `tau`（见 `RelaxationTime.relaxation_times`）。
4) （可选）通过 `ThermoDerivatives.bulk_viscosity_coefficients` 生成体粘滞需要的等熵系数。
5) 调用 `TransportCoefficients.transport_coefficients` 返回 (η, ζ, σ)。

单位约定：
- `T_fm`, `mu_fm`, 质量/动量均为 fm⁻¹
- `tau` 为 fm

注意：仓库目前采用 include 组织方式，本模块会 include 相关依赖文件。
- `prefer_energy_aniso`: ξ≠0 时是否优先走“能量直通”分布路径。
    - `true`：复用已计算的 `E_aniso`，直接调用 `quark_distribution(E,...) / antiquark_distribution(E,...)`（更省一次 `sqrt`）。
    - `false`：优先使用 provider 的 `*_distribution_aniso(p,m,...)`（适合你实现了非平凡 aniso 分布接口的情况）。
    - 可通过 `transport_kwargs.prefer_energy_aniso` 覆写；若两者都提供，以显式 keyword 为准。
"""

# --- 依赖：通过 RelaxTime.jl 入口加载全部 relaxtime 子模块 ---
const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end
const RelaxationTime = Main.RelaxationTime
const AFieldBuilder = Main.AFieldBuilder

# Shared config loader (TOML)
const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

const _PREFER_ENERGY_ANISO_CACHE_LOCK = ReentrantLock()
const _PREFER_ENERGY_ANISO_CACHE = Dict{String, Bool}()
const _A_BUILDER_CONFIG_CACHE_LOCK = ReentrantLock()
const _A_BUILDER_CONFIG_CACHE = Dict{String, NamedTuple}()

"""reset_transport_workflow_config_cache!()

仅供测试/调试：清空 workflow 内部对 toml 默认值的轻量缓存。

用途：
- 在同一 Julia session 内切换 `PHYSICS_PARAM_PROFILE` 后，强制让 workflow 重新读取 `config/physics/<profile>.toml`。
"""
function reset_transport_workflow_config_cache!()
    lock(_PREFER_ENERGY_ANISO_CACHE_LOCK)
    try
        empty!(_PREFER_ENERGY_ANISO_CACHE)
    finally
        unlock(_PREFER_ENERGY_ANISO_CACHE_LOCK)
    end

    lock(_A_BUILDER_CONFIG_CACHE_LOCK)
    try
        empty!(_A_BUILDER_CONFIG_CACHE)
    finally
        unlock(_A_BUILDER_CONFIG_CACHE_LOCK)
    end
    return nothing
end

# Unified equilibrium facade (solve_gap + state_vector + masses): reuse Main.* to avoid module duplication.
const _EQUILIBRIUM_FACADE_PATH = normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "EquilibriumFacade.jl"))
if !isdefined(Main, :EquilibriumFacade)
    Base.include(Main, _EQUILIBRIUM_FACADE_PATH)
end
const EquilibriumFacade = Main.EquilibriumFacade

# Shared parameter structs (QuarkParams/ThermoParams)
const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "..", "types", "ParameterTypes.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PARAMETER_TYPES_PATH)
end
using Main.ParameterTypes: QuarkParams, ThermoParams

const _WORKFLOW_PARAM_ADAPTERS_PATH = normpath(joinpath(@__DIR__, "..", "..", "models", "workflows", "WorkflowParamAdapters.jl"))
if !isdefined(Main, :WorkflowParamAdapters)
    Base.include(Main, _WORKFLOW_PARAM_ADAPTERS_PATH)
end
const WorkflowParamAdapters = Main.WorkflowParamAdapters
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params, as_legacy_inputs

# TransportCoefficients available via RelaxTime backward-compat promotion
const TransportCoefficients = Main.TransportCoefficients

using StaticArrays

const HADRON_SEED_5 = Main.Models.HADRON_SEED_5
const DEFAULT_MOMENTUM_COUNT = Main.Models.PNJLCore.DEFAULT_MOMENTUM_COUNT
const DEFAULT_THETA_COUNT = Main.Models.PNJLCore.DEFAULT_THETA_COUNT
const bulk_viscosity_coefficients = Main.Models.bulk_viscosity_coefficients
const DEFAULT_MOMENTUM_NODES = Main.Models.PNJLIntegrals.THERMAL_DEFAULT_NODES
const DEFAULT_MOMENTUM_WEIGHTS = Main.Models.PNJLIntegrals.THERMAL_DEFAULT_WEIGHTS
using Main.RelaxationTime: relaxation_times
using .TransportCoefficients: transport_coefficients, TransportIntegrationConfig
using .TransportCoefficients: rho_mass_from_densities

const PNJL = Main.Models

const _MODEL_CACHE = Dict{Symbol, Any}()

@inline function _get_model(model_kind::Symbol)
    return get!(_MODEL_CACHE, model_kind) do
        Main.Models.create_model(model_kind)
    end
end

export solve_gap_and_transport, build_equilibrium_params
export TransportIntegrationConfig
export solve_transport_from_equilibrium

const TRANSPORT_INTEGRATION_KEYS = (
    :p_nodes, :p_max,
    :p_grid, :p_w,
    :cos_nodes, :cos_grid, :cos_w,
)

const TRANSPORT_PROVIDER_KEYS = (
    :prefer_energy_aniso,
)

const A_BUILDER_KEYS = (
    :p_nodes,
    :p_max,
    :cos_nodes,
    :use_aniso,
)

@inline function _drop_transport_integration_keys(kwargs::NamedTuple)::NamedTuple
    return (; (k => v for (k, v) in pairs(kwargs) if !(k in TRANSPORT_INTEGRATION_KEYS))...)
end

@inline function _extract_transport_integration_kwargs(kwargs::NamedTuple)::NamedTuple
    return (; (k => v for (k, v) in pairs(kwargs) if (k in TRANSPORT_INTEGRATION_KEYS))...)
end

@inline function _provider_prefer_energy_aniso_from_kwargs(kwargs::NamedTuple)
    return hasproperty(kwargs, :prefer_energy_aniso) ? kwargs.prefer_energy_aniso : nothing
end

@inline function _default_prefer_energy_aniso_from_toml()::Bool
    physics_profile = get(ENV, "PHYSICS_PARAM_PROFILE", "default")

    lock(_PREFER_ENERGY_ANISO_CACHE_LOCK)
    try
        if haskey(_PREFER_ENERGY_ANISO_CACHE, physics_profile)
            return _PREFER_ENERGY_ANISO_CACHE[physics_profile]
        end
    finally
        unlock(_PREFER_ENERGY_ANISO_CACHE_LOCK)
    end

    physics_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "physics"))

    default_physics = Dict{String, Any}(
        "physical" => Dict(
            "hbarc" => 197.327,
            "alpha_em" => 1.0 / 137.035999084,
        ),
        "transport_workflow" => Dict(
            "prefer_energy_aniso" => true,
        ),
    )

    data = load_config(physics_dir, default_physics; profile=physics_profile)
    cfg = data.config
    tw = get(cfg, "transport_workflow", Dict{String, Any}())
    val = Bool(get(tw, "prefer_energy_aniso", true))

    lock(_PREFER_ENERGY_ANISO_CACHE_LOCK)
    try
        _PREFER_ENERGY_ANISO_CACHE[physics_profile] = val
    finally
        unlock(_PREFER_ENERGY_ANISO_CACHE_LOCK)
    end

    return val
end

@inline function _default_a_builder_config_from_toml()::NamedTuple
    physics_profile = get(ENV, "PHYSICS_PARAM_PROFILE", "default")

    lock(_A_BUILDER_CONFIG_CACHE_LOCK)
    try
        if haskey(_A_BUILDER_CONFIG_CACHE, physics_profile)
            return _A_BUILDER_CONFIG_CACHE[physics_profile]
        end
    finally
        unlock(_A_BUILDER_CONFIG_CACHE_LOCK)
    end

    physics_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "physics"))

    default_physics = Dict{String, Any}(
        "physical" => Dict(
            "hbarc" => 197.327,
            "alpha_em" => 1.0 / 137.035999084,
        ),
        "transport_workflow" => Dict(
            "prefer_energy_aniso" => true,
            "a_builder" => Dict(
                "p_nodes" => 16,
                "p_max" => 20.0,
                "cos_nodes" => 4,
                "use_aniso" => true,
            ),
        ),
    )

    data = load_config(physics_dir, default_physics; profile=physics_profile)
    cfg = data.config
    tw = get(cfg, "transport_workflow", Dict{String, Any}())
    a_builder = get(tw, "a_builder", Dict{String, Any}())

    val = (
        p_nodes=Int(get(a_builder, "p_nodes", 16)),
        p_max=Float64(get(a_builder, "p_max", 20.0)),
        cos_nodes=Int(get(a_builder, "cos_nodes", 4)),
        use_aniso=Bool(get(a_builder, "use_aniso", true)),
    )

    lock(_A_BUILDER_CONFIG_CACHE_LOCK)
    try
        _A_BUILDER_CONFIG_CACHE[physics_profile] = val
    finally
        unlock(_A_BUILDER_CONFIG_CACHE_LOCK)
    end

    return val
end

@inline function _merge_a_builder_config(base::NamedTuple, override::NamedTuple)::NamedTuple
    unknown = Symbol[]
    for k in keys(override)
        k in A_BUILDER_KEYS || push!(unknown, k)
    end
    isempty(unknown) || error("Unknown a_builder_config key(s): $(unknown). Allowed keys: $(A_BUILDER_KEYS)")

    return (
        p_nodes=get(override, :p_nodes, base.p_nodes),
        p_max=get(override, :p_max, base.p_max),
        cos_nodes=get(override, :cos_nodes, base.cos_nodes),
        use_aniso=get(override, :use_aniso, base.use_aniso),
    )
end

@inline function _effective_a_builder_config(a_builder_config::Union{Nothing,NamedTuple})::NamedTuple
    base = _default_a_builder_config_from_toml()
    return a_builder_config === nothing ? base : _merge_a_builder_config(base, a_builder_config)
end

@inline function _drop_transport_provider_keys(kwargs::NamedTuple)::NamedTuple
    return (; (k => v for (k, v) in pairs(kwargs) if !(k in TRANSPORT_PROVIDER_KEYS))...)
end

@inline function _apply_prefer_energy_aniso(provider, prefer_energy_aniso)
    prefer_energy_aniso === nothing && return provider

    if isdefined(Main, :Models) && isdefined(Main.Models, :TransportProvider) && provider isa Main.Models.TransportProvider
        return Main.Models.TransportProvider(
            provider.energy_from_p,
            provider.energy_from_p_aniso,
            provider.quark_distribution,
            provider.antiquark_distribution,
            provider.quark_distribution_aniso,
            provider.antiquark_distribution_aniso,
            prefer_energy_aniso,
            provider.mass_for_species,
            provider.mu_for_species,
            provider.ctx,
        )
    end

    # NamedTuple provider (or other mergeable tuples)
    if provider isa NamedTuple
        return merge(provider, (prefer_energy_aniso=prefer_energy_aniso,))
    end

    return provider
end

@inline function _default_transport_provider_for_backend()
    if isdefined(Main, :Models) && isdefined(Main.Models, :transport_provider)
        m = _get_model(:PNJL)
        try
            return Main.Models.transport_provider(m)
        catch
            return nothing
        end
    end
    return nothing
end

"""将平衡求解结果转换成 (quark_params, thermo_params)。"""
function build_equilibrium_params(base, T_fm::Real, mu_fm::Real; xi::Real=0.0)
    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    return (
        quark_params=QuarkParams((
            m=(u=NaN, d=NaN, s=NaN),
            μ=(u=Float64(mu_fm), d=Float64(mu_fm), s=Float64(mu_fm)),
        )),
        thermo_params=ThermoParams((T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi))),
    )
end

@inline function _densities_from_equilibrium(x_state, mu_vec, T_fm, thermal_nodes, xi; p_num::Int, t_num::Int)
    nd = Main.Models.number_densities(
        _get_model(:PNJL),
        x_state,
        T_fm,
        mu_vec,
        ;
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
    return (
        u=Float64(nd.quark[1]),
        d=Float64(nd.quark[2]),
        s=Float64(nd.quark[3]),
        ubar=Float64(nd.antiquark[1]),
        dbar=Float64(nd.antiquark[2]),
        sbar=Float64(nd.antiquark[3]),
    )
end

@inline function _thermo_background_from_equilibrium(base, T_fm::Real; xi::Real, p_num::Int, t_num::Int)
    if hasproperty(base, :pressure) && hasproperty(base, :energy) && hasproperty(base, :entropy)
        return (pressure=Float64(base.pressure), entropy=Float64(base.entropy), energy=Float64(base.energy))
    end

    pressure, _, entropy, energy = Main.Models.model_thermo(
        _get_model(:PNJL),
        base.x_state,
        base.mu_vec,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        xi=xi,
    )
    return (pressure=Float64(pressure), entropy=Float64(entropy), energy=Float64(energy))
end

@inline _solve_equilibrium(args...; kwargs...) = EquilibriumFacade.solve_equilibrium_backend(args...; kwargs...)

@inline function _transport_inputs_from_equilibrium(
    base,
    T_fm::Real,
    mu_fm::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
)
    masses = base.masses

    params0 = build_equilibrium_params(base, T_fm, mu_fm; xi=xi)
    quark_params_basic = QuarkParams((
        m=(u=Float64(masses[1]), d=Float64(masses[2]), s=Float64(masses[3])),
        μ=params0.quark_params.μ,
    ))
    thermo_params = params0.thermo_params

    densities = _densities_from_equilibrium(base.x_state, base.mu_vec, T_fm, nothing, Float64(xi);
        p_num=p_num,
        t_num=t_num,
    )

    return (masses=masses, quark_params=quark_params_basic, thermo_params=thermo_params, densities=densities)
end

@inline function _A_from_equilibrium(T_fm::Real, quark_params, thermo_params;
                                     a_builder_config::Union{Nothing,NamedTuple}=nothing)
    qp = normalize_quark_params(quark_params)
    tp = normalize_thermo_params(thermo_params)

    cfg = _effective_a_builder_config(a_builder_config)
    return AFieldBuilder.build_A_triplet(
        qp,
        tp;
        p_nodes=cfg.p_nodes,
        p_max=cfg.p_max,
        cos_nodes=cfg.cos_nodes,
        use_aniso=cfg.use_aniso,
    )
end

"""一次性完成：平衡求解 →（可选）τ 计算 →（可选）ζ 导数 → 输运系数。

返回值（NamedTuple）包含：
- `equilibrium`: 求解输出（pressure/energy/rho/.../x_state 等）
- `quark_params`, `thermo_params`
- `masses`: 三味有效质量
- `densities`: 用于 τ 计算的六种粒子/反粒子数密度
- `tau`, `tau_inv`: 弛豫时间及其倒数（若内部计算或调用方提供）
- `rates`: 平均散射率（若内部计算 τ 则返回，便于复用/诊断）
- `bulk_coeffs`: 若 `compute_bulk=true`，则给出体粘滞需要的导数组合
- `transport`: (eta, zeta, sigma)

关键词参数：
- `xi`: 各向异性参数 ξ
- `compute_tau`: 是否在内部计算 τ（需要 `K_coeffs`）
- `K_coeffs`: τ 所需的有效耦合系数（传给截面/平均散射率）
- `tau`: 若提供则直接用于输运计算，并仍会被原样回传（适合你自己在外部算 τ 或做对照）
- `compute_bulk`: 是否计算 ζ 需要的导数（会额外触发多次自动微分+求解，较慢）
- `p_num`, `t_num`: 热积分节点数（传给能隙求解/密度计算）
- `solver_kwargs`: 透传到 `solve`（例如 `iterations` 等）
- `tau_kwargs`: 透传到 `RelaxationTime.relaxation_times`（例如 `p_nodes/angle_nodes/phi_nodes/n_sigma_points/cs_caches/existing_rates` 等）
- `a_builder_config`: A 构造配置（`p_nodes/p_max/cos_nodes/use_aniso`），用于 `compute_tau=true` 时构建 `quark_params.A`。
- `transport_config`: 输运系数积分配置（推荐），例如 `TransportIntegrationConfig(p_nodes=64, p_max=15.0, cos_nodes=32)`。
    - 若同时在 `transport_kwargs` 里提供了 `p_nodes/p_max/...`，会被自动提取并用于构造 config（便于平滑迁移）。
- `c_p`, `rho_mass`: 可选热背景量；若提供则 `transport.prandtl_number` 可直接给出有限值。
- `transport_kwargs`: 透传到 `transport_coefficients` 的其它参数（建议只放 `degeneracy/charges` 等非积分配置项）。
"""
function solve_gap_and_transport(
    T_fm::Real,
    mu_fm::Real;
    xi::Real=0.0,
    solver_backend::Symbol=:legacy,
    equilibrium::Union{Nothing,Any}=nothing,
    compute_tau::Bool=false,
    K_coeffs::Union{Nothing,NamedTuple}=nothing,
    tau::Union{Nothing,NamedTuple}=nothing,
    compute_bulk::Bool=true,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    seed_state=HADRON_SEED_5,
    solver_kwargs::NamedTuple=(;),
    models_solver=nothing,
    models_residual_norm_max::Real=1e-4,
    tau_kwargs::NamedTuple=(;),
    a_builder_config::Union{Nothing,NamedTuple}=nothing,
    transport_config::Union{Nothing,TransportIntegrationConfig}=nothing,
    c_p::Union{Nothing,Real}=nothing,
    rho_mass::Union{Nothing,Real}=nothing,
    transport_kwargs::NamedTuple=(;),
    provider=nothing,
    prefer_energy_aniso=nothing
)
    base = equilibrium === nothing ? _solve_equilibrium(
        T_fm,
        mu_fm;
        xi=xi,
        solver_backend=solver_backend,
        p_num=p_num,
        t_num=t_num,
        seed_state=seed_state,
        solver_kwargs=solver_kwargs,
        models_solver=models_solver,
        models_residual_norm_max=models_residual_norm_max,
    ) : equilibrium

    return solve_transport_from_equilibrium(
        base,
        T_fm,
        mu_fm;
        xi=xi,
        compute_tau=compute_tau,
        K_coeffs=K_coeffs,
        tau=tau,
        compute_bulk=compute_bulk,
        p_num=p_num,
        t_num=t_num,
        tau_kwargs=tau_kwargs,
        a_builder_config=a_builder_config,
        transport_config=transport_config,
        c_p=c_p,
        rho_mass=rho_mass,
        transport_kwargs=transport_kwargs,
        provider=provider,
        prefer_energy_aniso=prefer_energy_aniso,
    )
end

"""solve_transport_from_equilibrium(equilibrium, T_fm, mu_fm; ...) -> NamedTuple

纯后处理层：在已给定平衡态解的前提下，完成密度/τ/ζ 导数与输运系数计算。
用于阶段 4 的 workflow 解耦：允许外部独立求平衡态，然后复用同一输运计算链。
"""
function solve_transport_from_equilibrium(
    base,
    T_fm::Real,
    mu_fm::Real;
    xi::Real=0.0,
    compute_tau::Bool=false,
    K_coeffs::Union{Nothing,NamedTuple}=nothing,
    tau::Union{Nothing,NamedTuple}=nothing,
    compute_bulk::Bool=true,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    tau_kwargs::NamedTuple=(;),
    a_builder_config::Union{Nothing,NamedTuple}=nothing,
    transport_config::Union{Nothing,TransportIntegrationConfig}=nothing,
    c_p::Union{Nothing,Real}=nothing,
    rho_mass::Union{Nothing,Real}=nothing,
    transport_kwargs::NamedTuple=(;),
    provider=nothing,
    prefer_energy_aniso=nothing
)
    inputs = _transport_inputs_from_equilibrium(base, T_fm, mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    masses = inputs.masses
    quark_params_basic = inputs.quark_params
    thermo_params = inputs.thermo_params
    densities = inputs.densities
    thermo_background = _thermo_background_from_equilibrium(base, T_fm; xi=xi, p_num=p_num, t_num=t_num)
    rho_mass_background = rho_mass_from_densities(quark_params_basic.m, densities)

    tau_inv = nothing
    rates = nothing

    if compute_tau
        K_coeffs === nothing && error("compute_tau=true requires K_coeffs")

        # 为截面/传播子准备 A 字段（TotalPropagator 会用到）
        A_vals = _A_from_equilibrium(T_fm, quark_params_basic, thermo_params; a_builder_config=a_builder_config)
        quark_params_full = (m=quark_params_basic.m, μ=quark_params_basic.μ, A=A_vals)

        tau_res = relaxation_times(
            quark_params_full,
            thermo_params,
            K_coeffs;
            densities=densities,
            tau_kwargs...,
        )
        tau = tau_res.tau
        tau_inv = tau_res.tau_inv
        rates = tau_res.rates
    else
        # 若不内部计算 τ，则要求调用方提供 τ
        tau === nothing && error("either provide tau=... or set compute_tau=true")
    end

    bulk_coeffs = nothing
    if compute_bulk
        bulk_coeffs = try
            bulk_viscosity_coefficients(
                T_fm,
                mu_fm;
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        catch bc_err
            @warn "bulk_viscosity_coefficients failed — ζ will be NaN for this point" T_fm=T_fm mu_fm=mu_fm xi=xi err=bc_err
            nothing
        end
    end

    c_p_background = if c_p !== nothing
        Float64(c_p)
    elseif bulk_coeffs !== nothing && hasproperty(bulk_coeffs, :c_p)
        Float64(bulk_coeffs.c_p)
    else
        NaN
    end
    rho_mass_effective = rho_mass === nothing ? rho_mass_background : Float64(rho_mass)
    thermo_background = merge(thermo_background, (rho_mass=rho_mass_effective, c_p=c_p_background))

    # Backward/ergonomic compatibility:
    prefer_from_kwargs = _provider_prefer_energy_aniso_from_kwargs(transport_kwargs)
    prefer_effective = if prefer_energy_aniso !== nothing
        prefer_energy_aniso
    elseif prefer_from_kwargs !== nothing
        prefer_from_kwargs
    else
        _default_prefer_energy_aniso_from_toml()
    end

    effective_provider = provider === nothing ? _default_transport_provider_for_backend() : provider
    legacy_inputs = as_legacy_inputs(quark_params_basic, thermo_params)
    if effective_provider === nothing
        # No backend default provider: only materialize a provider if the desired
        # behavior differs from TransportCoefficients default.
        if prefer_effective != true
            effective_provider = merge(TransportCoefficients.default_transport_provider(), (prefer_energy_aniso=prefer_effective,))
        end
    else
        effective_provider = _apply_prefer_energy_aniso(effective_provider, prefer_effective)
    end
    if effective_provider !== nothing && isdefined(Main, :Models) && isdefined(Main.Models, :TransportProvider)
        if effective_provider isa Main.Models.TransportProvider && isdefined(Main.Models, :prepare_transport_provider)
            effective_provider = Main.Models.prepare_transport_provider(
                effective_provider,
                base;
                quark_params=legacy_inputs.quark_params,
                thermo_params=legacy_inputs.thermo_params,
                masses=masses,
            )
        end
    end
    integration_kwargs = _extract_transport_integration_kwargs(transport_kwargs)
    transport_kwargs_clean = _drop_transport_integration_keys(transport_kwargs)
    transport_kwargs_clean = _drop_transport_provider_keys(transport_kwargs_clean)
    if effective_provider !== nothing
        transport_kwargs_clean = (; (k => v for (k, v) in pairs(transport_kwargs_clean) if k != :provider)...)
    end
    effective_transport_config = if transport_config !== nothing
        transport_config
    elseif length(keys(integration_kwargs)) > 0
        TransportIntegrationConfig(; integration_kwargs...)
    else
        nothing
    end

    pass_provider = effective_provider === nothing ? (; ) : (; provider=effective_provider)

    tr = if effective_transport_config === nothing
        transport_coefficients(
            legacy_inputs.quark_params,
            legacy_inputs.thermo_params;
            tau=tau,
            bulk_coeffs=bulk_coeffs,
            densities=densities,
            pressure=thermo_background.pressure,
            entropy=thermo_background.entropy,
            energy=thermo_background.energy,
            c_p=thermo_background.c_p,
            rho_mass=thermo_background.rho_mass,
            pass_provider...,
            transport_kwargs_clean...,
        )
    else
        transport_coefficients(
            legacy_inputs.quark_params,
            legacy_inputs.thermo_params;
            tau=tau,
            bulk_coeffs=bulk_coeffs,
            densities=densities,
            pressure=thermo_background.pressure,
            entropy=thermo_background.entropy,
            energy=thermo_background.energy,
            c_p=thermo_background.c_p,
            rho_mass=thermo_background.rho_mass,
            config=effective_transport_config,
            pass_provider...,
            transport_kwargs_clean...,
        )
    end

    return (
        equilibrium=base,
        quark_params=quark_params_basic,
        thermo_params=thermo_params,
        thermo_background=thermo_background,
        masses=masses,
        densities=densities,
        tau=tau,
        tau_inv=tau_inv,
        rates=rates,
        bulk_coeffs=bulk_coeffs,
        transport=tr,
    )
end

end # module
