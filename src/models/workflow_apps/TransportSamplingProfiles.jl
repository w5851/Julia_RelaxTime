module TransportSamplingProfiles

using TOML, SHA
using ..Models: create_model, PNJLModel
using ..TransportWorkflow: WorkflowCache, build_equilibrium_params, solve_transport_from_equilibrium,
    _transport_inputs_from_equilibrium

const RT = Main.RelaxTime
const AS = Main.AverageScatteringRate
const PROFILE_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "pnjl", "transport_sampling"))

struct SamplingProfile{B,A,Q,T}
    id::String
    source_path::String
    source_sha256::String
    background::B
    a_integral::A
    tau_integral::Q
    terminal_integral::T
    sampling::AS.ThermalSamplingConfig
    physics::NamedTuple{(:a_method,:b0_direction),Tuple{Symbol,Symbol}}
end

function _section(data, name, int_keys, real_keys)
    raw = data[name]
    Set(keys(raw)) == Set(String.((int_keys..., real_keys...))) ||
        throw(ArgumentError("unsupported or missing fields in sampling profile section $name"))
    integers = map(int_keys) do key
        value = raw[String(key)]
        value isa Integer && !(value isa Bool) && value > 0 ||
            throw(ArgumentError("$name.$key must be a positive integer"))
        Int(value)
    end
    reals = map(real_keys) do key
        value = Float64(raw[String(key)])
        isfinite(value) && value > 0 || throw(ArgumentError("$name.$key must be finite and positive"))
        value
    end
    return NamedTuple{(int_keys..., real_keys...)}((integers..., reals...))
end

function load_profile(name_or_path::AbstractString)
    path = isfile(name_or_path) ? abspath(name_or_path) : joinpath(PROFILE_ROOT, String(name_or_path) * ".toml")
    isfile(path) || throw(ArgumentError("transport sampling profile does not exist: $name_or_path"))
    bytes = read(path)
    data = TOML.parse(String(copy(bytes)))
    required = Set(("schema_version", "id", "algorithm", "background", "a_integral",
        "tau_integral", "terminal_integral", "sampling"))
    fields = Set(keys(data))
    required ⊆ fields && fields ⊆ union(required, Set(("physics",))) ||
        throw(ArgumentError("unsupported sampling profile fields"))
    data["schema_version"] == 1 || throw(ArgumentError("unsupported sampling profile schema"))
    data["algorithm"] == "thermal_weighted" || throw(ArgumentError("sampling profile requires thermal_weighted"))
    id = String(data["id"])
    isempty(strip(id)) && throw(ArgumentError("sampling profile id must not be empty"))
    background = _section(data, "background", (:p_nodes, :cos_nodes), (:thermal_p_max_inv_fm, :residual_norm_max))
    a = _section(data, "a_integral", (:p_nodes, :cos_nodes), (:scale_inv_fm,))
    q = _section(data, "tau_integral", (:p_nodes, :cos_nodes, :phi_nodes, :t_nodes), (:scale_inv_fm,))
    terminal = _section(data, "terminal_integral", (:p_nodes, :cos_nodes), (:scale_inv_fm,))
    sampling_values = Dict{Symbol,Any}(Symbol(k) => v for (k,v) in data["sampling"])
    if haskey(sampling_values, :relative_velocity)
        sampling_values[:relative_velocity] = Symbol(sampling_values[:relative_velocity])
    end
    sampling = AS.ThermalSamplingConfig(; sampling_values...)
    raw_physics = get(data, "physics", Dict("a_method"=>"quadrature", "b0_direction"=>"aligned"))
    Set(keys(raw_physics)) == Set(("a_method", "b0_direction")) ||
        throw(ArgumentError("physics requires a_method and b0_direction"))
    a_method = RT.AFieldBuilder._normalize_a_build_method(Symbol(raw_physics["a_method"]))
    a_method in (:quadrature, :ward_matched) ||
        throw(ArgumentError("sampling profile requires quadrature or ward_matched A"))
    b0_direction = RT.OneLoopIntegralsCorrection.normalize_b0_direction(Symbol(raw_physics["b0_direction"]))
    physics = (; a_method, b0_direction)
    return SamplingProfile(id, abspath(path), bytes2hex(sha256(bytes)), background, a, q, terminal, sampling, physics)
end

function profile_model(profile::SamplingProfile; template=create_model(:PNJL),
    thermal_p_max_inv_fm::Real=profile.background.thermal_p_max_inv_fm)
    isfinite(thermal_p_max_inv_fm) && thermal_p_max_inv_fm>0 ||
        throw(ArgumentError("thermal_p_max_inv_fm must be positive and finite"))
    p = template.params
    values = NamedTuple{fieldnames(typeof(p))}(Tuple(getfield(p, k) for k in fieldnames(typeof(p))))
    return PNJLModel(typeof(p)(; merge(values,
        (thermal_p_max_inv_fm=Float64(thermal_p_max_inv_fm),))...))
end

function effective_config(profile::SamplingProfile)
    return (id=profile.id, source_path=profile.source_path, source_sha256=profile.source_sha256,
        algorithm="thermal_weighted", background=profile.background, a_integral=profile.a_integral,
        tau_integral=profile.tau_integral, terminal_integral=profile.terminal_integral,
        sampling=AS._thermal_sampling_identity(profile.sampling),
        a_polarization=profile.physics.a_method===:ward_matched ? "ward_matched_RS" : "full_RS",
        a_coupling=profile.physics.a_method===:ward_matched ? "ward_matched_RS" : "full_RS",
        b0_direction=String(profile.physics.b0_direction),
        sigma_cutoff_inv_fm=Main.Constants_PNJL.Λ_inv_fm,
        internal_density=(p_nodes=AS.DEFAULT_SEMI_INF_NODES, scale_inv_fm=AS.DEFAULT_SEMI_INF_SCALE,
            cos_nodes=profile.tau_integral.cos_nodes),
        external_density=(p_nodes=64, cos_nodes=8, p_max_inv_fm=profile.background.thermal_p_max_inv_fm),
        cache_lifecycle="REQUIRED_RATE_KEYS_FOR_TAU order; reset polarization cache before each channel")
end

function build_runtime(profile::SamplingProfile)
    model = profile_model(profile)
    cache = WorkflowCache()
    cache.model_cache[:PNJL] = model
    q = profile.tau_integral
    p, pw = AS._cached_semi_infinite_momentum_grid(q.p_nodes, q.scale_inv_fm)
    c, cw = AS._cached_interval_grid(-1.0, 1.0, q.cos_nodes)
    phi, phiw = AS._cached_interval_grid(0.0, 2pi, q.phi_nodes)
    t = profile.terminal_integral
    tp, tw = AS._cached_semi_infinite_momentum_grid(t.p_nodes, t.scale_inv_fm)
    tc, tcw = AS._cached_interval_grid(-1.0, 1.0, t.cos_nodes)
    integration = RT.TransportCoefficients.TransportIntegrationConfig(p_grid=tp,p_w=tw,cos_grid=tc,cos_w=tcw)
    tau = (p_nodes=q.p_nodes,angle_nodes=q.cos_nodes,phi_nodes=q.phi_nodes,n_sigma_points=q.t_nodes,
        p_grid=p,p_w=pw,cos_grid=c,cos_w=cw,phi_grid=phi,phi_w=phiw,
        sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,sampling_config=profile.sampling)
    return (profile=profile,model=model,workflow_cache=cache,tau_kwargs=tau,transport_config=integration,
        p_grid=p,p_w=pw,cos_grid=c,cos_w=cw,phi_grid=phi,phi_w=phiw,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm)
end

function point_state(profile::SamplingProfile, equilibrium, T_inv_fm, mu_inv_fm; xi=0.0,
    propagator_xi_policy::Symbol=:match_thermo)
    propagator_xi_policy in (:match_thermo, :isotropic) || throw(ArgumentError("invalid propagator_xi_policy"))
    parameters = build_equilibrium_params(equilibrium,T_inv_fm,mu_inv_fm;xi=xi)
    masses = NamedTuple{(:u,:d,:s)}(Tuple(Float64.(equilibrium.masses)))
    all(m->isfinite(m) && m>0,masses) || throw(ArgumentError("equilibrium masses must be positive and finite"))
    qp = merge(Main.ParameterAdapters.normalize_quark_input(parameters.quark_params),(m=masses,))
    tp = Main.ParameterAdapters.normalize_thermo_input(parameters.thermo_params)
    return fixed_point_state(profile, qp, tp; propagator_xi_policy)
end

"""Construct loop inputs on a supplied background without solving mean field."""
function fixed_point_state(profile::SamplingProfile, qp::NamedTuple, tp::NamedTuple;
    propagator_xi_policy::Symbol=:match_thermo)
    AS._normalize_propagator_xi_policy(propagator_xi_policy)
    if profile.physics.b0_direction !== :aligned
        tp = merge(tp, (b0_direction=profile.physics.b0_direction,))
    end
    a = profile.a_integral
    p,w = AS._cached_semi_infinite_momentum_grid(a.p_nodes,a.scale_inv_fm)
    av = RT.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=a.cos_nodes,use_aniso=true,
        a_method=profile.physics.a_method)
    iso = RT.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=a.cos_nodes,use_aniso=false)
    gu = RT.EffectiveCouplings.calculate_G_from_A(av.u,qp.m.u)
    gs = RT.EffectiveCouplings.calculate_G_from_A(av.s,qp.m.s)
    kc = RT.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    return (quark_params=merge(qp,(A=av,)),thermo_params=tp,K_coeffs=kc,
        propagator_quark_params=merge(qp,(A=propagator_xi_policy===:isotropic ? iso : av,)))
end

"""
    solve_frozen_transport(runtime, qp, tp, densities, entropy; bulk_coeffs, on_channel)

Recompute all 17 scattering channels and terminal transport at a frozen
background. Densities, entropy and isentropic derivatives are supplied by the
caller and never replaced by a new phase or mean-field solve.
"""
function solve_frozen_transport(runtime, qp::NamedTuple, tp::NamedTuple,
    densities::NamedTuple, entropy::Real; bulk_coeffs=nothing, on_channel=nothing,
    existing_rates::AbstractDict=Dict{Symbol,Float64}())
    all(v -> isfinite(v) && v > 0, values(densities)) ||
        throw(ArgumentError("frozen species densities must be positive and finite"))
    isfinite(entropy) && entropy > 0 || throw(ArgumentError("frozen entropy must be positive and finite"))
    state = fixed_point_state(runtime.profile, qp, tp)
    rates = Dict{Symbol,Float64}()
    diagnostics = NamedTuple[]
    for process in RT.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
        if haskey(existing_rates, process)
            value = Float64(existing_rates[process])
            isfinite(value) && value >= 0 || throw(ArgumentError("invalid resumed channel rate: $process"))
            rates[process] = value
            continue
        end
        GC.gc(false)
        RT.PolarizationCache.reset_cache!()
        measured = @timed AS.build_thermal_weighted_cache(process, state.quark_params,
            state.thermo_params, state.K_coeffs; runtime.tau_kwargs...)
        cache = measured.value
        rates[process] = cache.rate
        row = merge((channel=String(process),seconds=measured.time,rate=cache.rate,
            estimated_error=cache.estimated_error), cache.statistics)
        push!(diagnostics, row)
        on_channel === nothing || on_channel(row, cache)
    end
    inverse = RT.RelaxationTime.relaxation_rates(densities, rates)
    tau = map(inv, inverse)
    config = runtime.transport_config
    eta = RT.TransportCoefficients.shear_viscosity(qp,tp; tau,config)
    sigma = RT.TransportCoefficients.electric_conductivity(qp,tp; tau,config)
    zeta = bulk_coeffs === nothing ? NaN : RT.TransportCoefficients.bulk_viscosity_isentropic(
        qp,tp; bulk_coeffs_isentropic=bulk_coeffs,tau,config)
    all(v -> isfinite(v) && v > 0, (values(tau)..., eta, sigma)) ||
        error("frozen transport returned invalid observables")
    bulk_coeffs === nothing || (isfinite(zeta) && zeta >= 0) ||
        error("frozen transport returned invalid bulk viscosity")
    return (; tau, tau_inv=inverse, rates, transport=(; eta,sigma,zeta),
        entropy=Float64(entropy), sampling_diagnostics=diagnostics,
        sampling_profile=effective_config(runtime.profile), loop_state=state)
end

function solve_profile_transport(runtime, equilibrium, T_inv_fm, mu_inv_fm; xi=0.0,
    compute_bulk::Bool=true, propagator_xi_policy::Symbol=:match_thermo, on_channel=nothing)
    profile = runtime.profile
    state = point_state(profile,equilibrium,T_inv_fm,mu_inv_fm;xi=xi,propagator_xi_policy=propagator_xi_policy)
    caches = Dict{Symbol,AS.ThermalWeightedCache}()
    diagnostics = NamedTuple[]
    options = merge(runtime.tau_kwargs,(propagator_xi_policy=propagator_xi_policy,
        propagator_quark_params=state.propagator_quark_params))
    for process in RT.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
        GC.gc()
        RT.PolarizationCache.reset_cache!()
        measured = @timed AS.build_thermal_weighted_cache(process,state.quark_params,state.thermo_params,
            state.K_coeffs;options...)
        cache = measured.value
        caches[process] = cache
        row = merge((channel=String(process),seconds=measured.time,rate=cache.rate,
            estimated_error=cache.estimated_error),cache.statistics)
        push!(diagnostics,row)
        on_channel === nothing || on_channel(row)
    end
    rates = Dict(process=>cache.rate for (process,cache) in caches)
    result = if profile.physics === (a_method=:quadrature,b0_direction=:aligned)
        kwargs = merge(options,(sigma_cache_policy=:thermal_weighted,interpolation_mode=:auto,cs_caches=caches))
        solve_transport_from_equilibrium(equilibrium,T_inv_fm,mu_inv_fm;xi=xi,
            compute_tau=true,compute_bulk=compute_bulk,p_num=profile.background.p_nodes,
            t_num=profile.background.cos_nodes,K_coeffs=state.K_coeffs,tau_kwargs=kwargs,
            transport_config=runtime.transport_config,workflow_cache=runtime.workflow_cache)
    else
        inputs = _transport_inputs_from_equilibrium(runtime.workflow_cache,
            equilibrium,T_inv_fm,mu_inv_fm;xi=xi,p_num=profile.background.p_nodes,
            t_num=profile.background.cos_nodes)
        inverse = RT.RelaxationTime.relaxation_rates(inputs.densities, rates)
        tau = map(inv, inverse)
        computed = solve_transport_from_equilibrium(equilibrium,T_inv_fm,mu_inv_fm;xi=xi,
            compute_tau=false,tau=tau,compute_bulk=compute_bulk,p_num=profile.background.p_nodes,
            t_num=profile.background.cos_nodes,K_coeffs=state.K_coeffs,
            transport_config=runtime.transport_config,workflow_cache=runtime.workflow_cache)
        merge(computed,(tau_inv=inverse,))
    end
    compute_bulk && result.bulk_coeffs===nothing && error("profile transport requires converged bulk derivatives")
    values = (Tuple(result.tau)...,result.transport.eta,result.transport.sigma,result.thermo_background.entropy)
    all(v->isfinite(v) && v>0,values) || error("profile transport returned invalid observables")
    compute_bulk && !(isfinite(result.transport.zeta) && result.transport.zeta>=0) &&
        error("profile transport returned invalid bulk viscosity")
    rate_values = NamedTuple{RT.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU}(
        Tuple(rates[k] for k in RT.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
    return merge(result,(rates=rate_values,sampling_diagnostics=diagnostics,
        sampling_profile=effective_config(profile)))
end

end # module TransportSamplingProfiles
