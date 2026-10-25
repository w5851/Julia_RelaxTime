# Included in AverageScatteringRate after the quadrature helpers are defined.

struct ThermalSamplingConfig
    rtol::Float64
    atol::Float64
    weighted_window::Float64
    seed_order::Int
    max_evaluations::Int
    max_nodes::Int
    max_rounds::Int
    max_events::Int
    relative_velocity::Symbol
    function ThermalSamplingConfig(;rtol=1e-3,atol=1e-12,weighted_window=0.6,
        seed_order=48,max_evaluations=40_000,max_nodes=12_000,max_rounds=40,
        max_events=30_000_000,relative_velocity=:cm)
        all(isfinite,(rtol,atol,weighted_window)) && rtol>=0 && atol>=0 &&
            rtol+atol>0 && weighted_window>=0 ||
            throw(ArgumentError("thermal sampling requires finite nonnegative tolerances and window"))
        seed_order>=8 && max_evaluations>0 && max_nodes>=2 && max_rounds>0 && max_events>0 ||
            throw(ArgumentError("invalid thermal sampling node/evaluation/event budget"))
        velocity_mode=_normalize_relative_velocity_policy(relative_velocity)
        new(rtol,atol,weighted_window,seed_order,max_evaluations,max_nodes,max_rounds,max_events,velocity_mode)
    end
end

const DEFAULT_THERMAL_SAMPLING_CONFIG = ThermalSamplingConfig()

"""A qualified table and its average on exactly the recorded thermal quadrature.

Its identity covers both external distributions and propagator inputs, all
effective quadrature nodes/weights, the physical cutoff and error budgets.
The large discrete measure is released after construction. This cache is
immutable, but callers must also treat its table arrays as read-only.
"""
struct ThermalWeightedCache{I,C,S}
    process::Symbol
    table::ThermalWeightedSampling.WeightedSigmaTable
    rate::Float64
    estimated_error::Float64
    identity::I
    context::C
    statistics::S
end

function _thermal_sampling_identity(config::ThermalSamplingConfig)
    return NamedTuple{fieldnames(ThermalSamplingConfig)}(
        Tuple(getfield(config,k) for k in fieldnames(ThermalSamplingConfig)))
end

function _prepare_thermal_cache_problem(process,quark_params,thermo_params,K_coeffs;
    p_nodes=DEFAULT_P_NODES,angle_nodes=DEFAULT_ANGLE_NODES,phi_nodes=DEFAULT_PHI_NODES,
    p_grid=nothing,p_w=nothing,cos_grid=nothing,cos_w=nothing,phi_grid=nothing,phi_w=nothing,
    scale=DEFAULT_SEMI_INF_SCALE,density_p_grid=nothing,density_p_w=nothing,
    density_p_nodes=DEFAULT_SEMI_INF_NODES,density_scale=DEFAULT_SEMI_INF_SCALE,
    sigma_cutoff=Λ_inv_fm,apply_s_domain_cut=true,n_sigma_points=DEFAULT_N_SIGMA_POINTS,
    propagator_xi_policy=:match_thermo,propagator_quark_params=nothing,
    sampling_config::ThermalSamplingConfig=DEFAULT_THERMAL_SAMPLING_CONFIG)
    qp=normalize_quark_input(quark_params);tp=normalize_thermo_input(thermo_params)
    policy=_normalize_propagator_xi_policy(propagator_xi_policy)
    pt=_resolve_propagator_thermo_params(tp,policy)
    pq=_resolve_propagator_quark_params(qp,pt,policy,propagator_quark_params)
    pq=ensure_quark_params_has_A(pq,pt;cos_nodes=DEFAULT_A_COS_NODES,warn_on_auto=false)
    all(n->n>0,(p_nodes,angle_nodes,phi_nodes,density_p_nodes,n_sigma_points)) ||
        throw(ArgumentError("thermal quadrature node counts must be positive"))
    all(v->isfinite(v) && v>0,(scale,density_scale,tp.T)) && isfinite(tp.ξ) && tp.ξ>-1 ||
        throw(ArgumentError("thermal cache requires positive finite scales/T and xi > -1"))
    apply_s_domain_cut || throw(ArgumentError("thermal-weighted sampling requires the physical s-domain cut"))
    cutoff=sigma_cutoff===nothing ? Λ_inv_fm : Float64(sigma_cutoff)
    isfinite(cutoff) && cutoff>0 || throw(ArgumentError("sigma_cutoff must be positive and finite"))
    for (gn,g,wn,w) in (("p_grid",p_grid,"p_w",p_w),("cos_grid",cos_grid,"cos_w",cos_w),
        ("phi_grid",phi_grid,"phi_w",phi_w),("density_p_grid",density_p_grid,"density_p_w",density_p_w))
        validate_grid_weight_pair("thermal-weighted cache",gn,g,wn,w)
    end
    p,pw=p_grid===nothing ? _cached_semi_infinite_momentum_grid(p_nodes,Float64(scale)) : (p_grid,p_w)
    c,cw=cos_grid===nothing ? _cached_interval_grid(-1.0,1.0,angle_nodes) : (cos_grid,cos_w)
    az,aw=phi_grid===nothing ? _cached_interval_grid(0.0,TWO_PI,phi_nodes) : (phi_grid,phi_w)
    dp,dw=density_p_grid===nothing ?
        _cached_semi_infinite_momentum_grid(density_p_nodes,Float64(density_scale)) : (density_p_grid,density_p_w)
    all(v->all(isfinite,v),(p,pw,c,cw,az,aw,dp,dw)) &&
        all(v->all(>(0),v),(pw,cw,aw,dw)) && all(>=(0),p) && all(>=(0),dp) && all(x->-1<=x<=1,c) ||
        throw(ArgumentError("thermal quadrature requires finite nodes, positive weights and physical cosines"))
    incoming_i,incoming_j,outgoing_i,outgoing_j=parse_particles_from_process(process)
    mi,mj,mc,md=(get_mass(s,qp) for s in (incoming_i,incoming_j,outgoing_i,outgoing_j))
    _,lo,hi=_resolve_sigma_support_bounds(mi,mj,mc,md,cutoff,p_grid,p_w)
    mui,muj=get_mu(incoming_i,qp),get_mu(incoming_j,qp)
    ei,ej=energy_from_p.(p,mi),energy_from_p.(p,mj)
    fi=[distribution_with_anisotropy(incoming_i,v,mi,mui,tp.T,tp.Φ,tp.Φbar,tp.ξ,x) for v in p,x in c]
    fj=[distribution_with_anisotropy(incoming_j,v,mj,muj,tp.T,tp.Φ,tp.Φbar,tp.ξ,x) for v in p,x in c]
    st=sqrt.(max.(1.0 .-c.*c,0.0));ca=cos.(az)
    cr=[c[a]*c[b]+st[a]*st[b]*ca[k] for k in eachindex(az),a in eachindex(c),b in eachindex(c)]
    ni=number_density(incoming_i,mi,mui,tp.T,tp.Φ,tp.Φbar,tp.ξ;
        p_grid=dp,p_w=dw,angle_nodes=angle_nodes)
    nj=number_density(incoming_j,mj,muj,tp.T,tp.Φ,tp.Φbar,tp.ξ;
        p_grid=dp,p_w=dw,angle_nodes=angle_nodes)
    isfinite(ni) && isfinite(nj) && ni>0 && nj>0 ||
        throw(ArgumentError("thermal-weighted cache requires positive finite incoming densities"))
    plan=(mi=mi,mj=mj,p=p,pw=pw,ei=ei,ej=ej,fi=fi,fj=fj,cw=cw,phiw=aw,
        cos_relative=cr,prefactor=DQ^2/(32pi^5*ni*nj))
    identity=(version=1,process=process,outer_quark=_quark_fingerprint(_base_quark_params_without_A(qp)),
        outer_thermo=_thermo_fingerprint(tp),propagator_quark=_quark_fingerprint(pq),
        propagator_thermo=_thermo_fingerprint(pt),couplings=_fingerprint_value(K_coeffs),
        propagator_xi_policy=policy,n_sigma_points=n_sigma_points,cutoff=cutoff,bounds=(lo,hi),
        momentum=(_s_grid_fingerprint(p),_s_grid_fingerprint(pw)),
        polar=(_s_grid_fingerprint(c),_s_grid_fingerprint(cw)),
        azimuth=(_s_grid_fingerprint(az),_s_grid_fingerprint(aw)),
        density=(_s_grid_fingerprint(dp),_s_grid_fingerprint(dw),angle_nodes),
        sampling=_thermal_sampling_identity(sampling_config))
    context=(quark_params=pq,thermo_params=pt,K_coeffs=K_coeffs,n_sigma_points=n_sigma_points,bounds=(lo,hi))
    return (plan=plan,identity=identity,context=context,sampling=sampling_config)
end

function _build_thermal_weighted_cache(process,problem)
    ts=ThermalWeightedSampling;cfg=problem.sampling;ctx=problem.context
    measure=ts.thermal_measure(problem.plan,ctx.bounds;max_events=cfg.max_events,
        relative_velocity=cfg.relative_velocity)
    state=(quark_params=ctx.quark_params,thermo_params=ctx.thermo_params,K_coeffs=ctx.K_coeffs)
    s0=first(ctx.bounds)
    features=ts.scattering_feature_seeds(process,state,measure,s0)
    qp=ctx.quark_params
    thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
    direct=s->total_cross_section(process,s,ctx.quark_params,ctx.thermo_params,ctx.K_coeffs;
        n_points=ctx.n_sigma_points)
    built=ts.build_thermal_weighted_table(direct,measure,s0;rtol=cfg.rtol,atol=cfg.atol,
        weighted_window=cfg.weighted_window,seed_order=cfg.seed_order,thresholds=thresholds,
        extra_seeds=features.seeds,max_evaluations=cfg.max_evaluations,max_nodes=cfg.max_nodes,
        max_rounds=cfg.max_rounds)
    table=ts.validated_table(built) # An exhausted budget never yields a production cache.
    stats=(evaluations=built.evaluations,nodes=length(table.nodes),events=measure.events,
        unique_events=length(measure.nodes),guard_phase=built.guard_phase,
        feature_count=length(features.features),reason=built.reason)
    return ThermalWeightedCache(process,table,built.rate,built.estimated_error,
        problem.identity,problem.context,stats)
end

"""Build a thermal-weighted cache; failure to reach the requested budget throws.

The error estimate controls interpolation at the specified outer quadrature.
It does not certify convergence of that quadrature, t integration or B0.
"""
function build_thermal_weighted_cache(process::Symbol,quark_params,thermo_params,K_coeffs;kwargs...)
    problem=_prepare_thermal_cache_problem(process,quark_params,thermo_params,K_coeffs;kwargs...)
    return _build_thermal_weighted_cache(process,problem)
end

function _thermal_weighted_average(process,qp,tp,kc,cache;kwargs...)
    problem=_prepare_thermal_cache_problem(process,qp,tp,kc;kwargs...)
    if cache===nothing
        cache=_build_thermal_weighted_cache(process,problem)
    elseif !(cache isa ThermalWeightedCache)
        throw(ArgumentError("thermal-weighted sampling requires a ThermalWeightedCache"))
    elseif !isequal(cache.identity,problem.identity)
        throw(ArgumentError("thermal cache context or quadrature/error budget changed; rebuild the cache"))
    end
    return cache.rate
end

"""Query the stored interpolant, or explicitly evaluate a local cross section.

The physical s-domain cut always applies. `direct=true` allows a diagnostic
query inside that domain but outside the discrete support certified by the
cache, without changing the cached table, rate or validity.
"""
function get_sigma(cache::ThermalWeightedCache,s::Float64;direct::Bool=false)
    isfinite(s) || throw(ArgumentError("s must be finite"))
    ctx=cache.context
    first(ctx.bounds)<s<last(ctx.bounds) || return 0.0
    direct || return cache.table(s)
    value=total_cross_section(cache.process,s,ctx.quark_params,ctx.thermo_params,ctx.K_coeffs;
        n_points=ctx.n_sigma_points)
    isfinite(value) && value>=0 || throw(DomainError((s,value),"invalid direct cross section"))
    return value
end

export ThermalSamplingConfig, ThermalWeightedCache, build_thermal_weighted_cache
