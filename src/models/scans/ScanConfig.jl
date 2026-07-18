"""
    ScanConfig

扫描配置对象：为 `run_tmu_scan` / `run_trho_scan` 提供结构化参数入口，
同时保持现有 kwargs 接口不变。
"""
module ScanConfig

export TmuScanConfig, TrhoScanConfig
export FreezeoutScanConfig
export scan_kwargs

Base.@kwdef struct TmuScanConfig
    T_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    mu_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    xi_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    output_path::Union{Nothing, AbstractString} = nothing
    overwrite::Union{Nothing, Bool} = nothing
    resume::Union{Nothing, Bool} = nothing
    use_phase_aware::Union{Nothing, Bool} = nothing
    bootstrap_multiseed::Union{Nothing, Bool} = nothing
    thermo_backend::Union{Nothing, Symbol} = nothing
    solver_backend::Union{Nothing, Symbol} = nothing
    auto_pnjl_backend::Union{Nothing, Symbol} = nothing
    semantic_mode::Union{Nothing, Symbol} = nothing
    selector::Union{Nothing, Function} = nothing
    p_num::Union{Nothing, Int} = nothing
    t_num::Union{Nothing, Int} = nothing
    thermo_quadrature_policy::Union{Nothing, Symbol} = nothing
    thermo_quadrature_rtol::Union{Nothing, Float64} = nothing
    thermo_quadrature_atol::Union{Nothing, Float64} = nothing
    thermo_quadrature_maxevals::Union{Nothing, Int} = nothing
    progress_cb::Union{Nothing, Function} = nothing
    nlsolve_kwargs::NamedTuple = (;)
end

Base.@kwdef struct TrhoScanConfig
    T_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    rho_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    xi_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    output_path::Union{Nothing, AbstractString} = nothing
    overwrite::Union{Nothing, Bool} = nothing
    resume::Union{Nothing, Bool} = nothing
    reverse_rho::Union{Nothing, Bool} = nothing
    seed_policy::Union{Nothing, Symbol} = nothing
    hybrid_weighted_fallback::Union{Nothing, Bool} = nothing
    hybrid_weighted_max_seed_candidates::Union{Nothing, Int} = nothing
    constraint_mode::Union{Nothing, Symbol} = nothing
    asym_ud_ratio_target::Union{Nothing, Float64} = nothing
    asym_s_target::Union{Nothing, Float64} = nothing
    thermo_backend::Union{Nothing, Symbol} = nothing
    solver_backend::Union{Nothing, Symbol} = nothing
    auto_pnjl_backend::Union{Nothing, Symbol} = nothing
    semantic_mode::Union{Nothing, Symbol} = nothing
    selector::Union{Nothing, Function} = nothing
    model_kind::Union{Nothing, Symbol} = nothing
    p_num::Union{Nothing, Int} = nothing
    t_num::Union{Nothing, Int} = nothing
    thermo_quadrature_policy::Union{Nothing, Symbol} = nothing
    thermo_quadrature_rtol::Union{Nothing, Float64} = nothing
    thermo_quadrature_atol::Union{Nothing, Float64} = nothing
    thermo_quadrature_maxevals::Union{Nothing, Int} = nothing
    progress_cb::Union{Nothing, Function} = nothing
    nlsolve_kwargs::NamedTuple = (;)
end

Base.@kwdef struct FreezeoutScanConfig
    sqrt_s_NN_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    xi_values::Union{Nothing, AbstractVector{<:Real}} = nothing
    output_path::Union{Nothing, AbstractString} = nothing
    overwrite::Union{Nothing, Bool} = nothing
    resume::Union{Nothing, Bool} = nothing
    bootstrap_multiseed::Union{Nothing, Bool} = nothing
    solver_backend::Union{Nothing, Symbol} = nothing
    auto_pnjl_backend::Union{Nothing, Symbol} = nothing
    semantic_mode::Union{Nothing, Symbol} = nothing
    selector::Union{Nothing, Function} = nothing
    model_kind::Union{Nothing, Symbol} = nothing
    p_num::Union{Nothing, Int} = nothing
    t_num::Union{Nothing, Int} = nothing
    profile_name::Union{Nothing, AbstractString} = nothing
    path_profile_name::Union{Nothing, AbstractString} = nothing
    traversal::Union{Nothing, Symbol} = nothing
    progress_cb::Union{Nothing, Function} = nothing
    nlsolve_kwargs::NamedTuple = (;)
end

@inline function _drop_nothing(values)
    return (; (k => v for (k, v) in pairs(values) if v !== nothing)...)
end

function scan_kwargs(cfg::TmuScanConfig)::NamedTuple
    return _drop_nothing((
        T_values=cfg.T_values,
        mu_values=cfg.mu_values,
        xi_values=cfg.xi_values,
        output_path=cfg.output_path,
        overwrite=cfg.overwrite,
        resume=cfg.resume,
        use_phase_aware=cfg.use_phase_aware,
        bootstrap_multiseed=cfg.bootstrap_multiseed,
        thermo_backend=cfg.thermo_backend,
        solver_backend=cfg.solver_backend,
        auto_pnjl_backend=cfg.auto_pnjl_backend,
        semantic_mode=cfg.semantic_mode,
        selector=cfg.selector,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        thermo_quadrature_policy=cfg.thermo_quadrature_policy,
        thermo_quadrature_rtol=cfg.thermo_quadrature_rtol,
        thermo_quadrature_atol=cfg.thermo_quadrature_atol,
        thermo_quadrature_maxevals=cfg.thermo_quadrature_maxevals,
        progress_cb=cfg.progress_cb,
    ))
end

function scan_kwargs(cfg::TrhoScanConfig)::NamedTuple
    return _drop_nothing((
        T_values=cfg.T_values,
        rho_values=cfg.rho_values,
        xi_values=cfg.xi_values,
        output_path=cfg.output_path,
        overwrite=cfg.overwrite,
        resume=cfg.resume,
        reverse_rho=cfg.reverse_rho,
        seed_policy=cfg.seed_policy,
        hybrid_weighted_fallback=cfg.hybrid_weighted_fallback,
        hybrid_weighted_max_seed_candidates=cfg.hybrid_weighted_max_seed_candidates,
        constraint_mode=cfg.constraint_mode,
        asym_ud_ratio_target=cfg.asym_ud_ratio_target,
        asym_s_target=cfg.asym_s_target,
        thermo_backend=cfg.thermo_backend,
        solver_backend=cfg.solver_backend,
        auto_pnjl_backend=cfg.auto_pnjl_backend,
        semantic_mode=cfg.semantic_mode,
        selector=cfg.selector,
        model_kind=cfg.model_kind,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        thermo_quadrature_policy=cfg.thermo_quadrature_policy,
        thermo_quadrature_rtol=cfg.thermo_quadrature_rtol,
        thermo_quadrature_atol=cfg.thermo_quadrature_atol,
        thermo_quadrature_maxevals=cfg.thermo_quadrature_maxevals,
        progress_cb=cfg.progress_cb,
    ))
end

function scan_kwargs(cfg::FreezeoutScanConfig)::NamedTuple
    return _drop_nothing((
        sqrt_s_NN_values=cfg.sqrt_s_NN_values,
        xi_values=cfg.xi_values,
        output_path=cfg.output_path,
        overwrite=cfg.overwrite,
        resume=cfg.resume,
        bootstrap_multiseed=cfg.bootstrap_multiseed,
        solver_backend=cfg.solver_backend,
        auto_pnjl_backend=cfg.auto_pnjl_backend,
        semantic_mode=cfg.semantic_mode,
        selector=cfg.selector,
        model_kind=cfg.model_kind,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        profile_name=cfg.profile_name,
        path_profile_name=cfg.path_profile_name,
        traversal=cfg.traversal,
        progress_cb=cfg.progress_cb,
    ))
end

end # module
