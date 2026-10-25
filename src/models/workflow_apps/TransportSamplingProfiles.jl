module TransportSamplingProfiles

using TOML, SHA
using ..Models: create_model, PNJLModel
using ..TransportWorkflow: WorkflowCache, build_equilibrium_params, solve_transport_from_equilibrium

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
    Set(keys(data)) == Set(("schema_version", "id", "algorithm", "background", "a_integral",
        "tau_integral", "terminal_integral", "sampling")) || throw(ArgumentError("unsupported sampling profile fields"))
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
    return SamplingProfile(id, abspath(path), bytes2hex(sha256(bytes)), background, a, q, terminal, sampling)
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
        a_polarization="full_RS", a_coupling="full_RS", sigma_cutoff_inv_fm=Main.Constants_PNJL.Λ_inv_fm,
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
    a = profile.a_integral
    p,w = AS._cached_semi_infinite_momentum_grid(a.p_nodes,a.scale_inv_fm)
    av = RT.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=a.cos_nodes,use_aniso=true)
    iso = RT.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=a.cos_nodes,use_aniso=false)
    gu = RT.EffectiveCouplings.calculate_G_from_A(av.u,qp.m.u)
    gs = RT.EffectiveCouplings.calculate_G_from_A(av.s,qp.m.s)
    kc = RT.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    return (quark_params=merge(qp,(A=av,)),thermo_params=tp,K_coeffs=kc,
        propagator_quark_params=merge(qp,(A=propagator_xi_policy===:isotropic ? iso : av,)))
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
    kwargs = merge(options,(sigma_cache_policy=:thermal_weighted,interpolation_mode=:auto,cs_caches=caches))
    result = solve_transport_from_equilibrium(equilibrium,T_inv_fm,mu_inv_fm;xi=xi,
        compute_tau=true,compute_bulk=compute_bulk,p_num=profile.background.p_nodes,
        t_num=profile.background.cos_nodes,K_coeffs=state.K_coeffs,tau_kwargs=kwargs,
        transport_config=runtime.transport_config,workflow_cache=runtime.workflow_cache)
    compute_bulk && result.bulk_coeffs===nothing && error("profile transport requires converged bulk derivatives")
    values = (Tuple(result.tau)...,result.transport.eta,result.transport.sigma,result.thermo_background.entropy)
    all(v->isfinite(v) && v>0,values) || error("profile transport returned invalid observables")
    compute_bulk && !(isfinite(result.transport.zeta) && result.transport.zeta>=0) &&
        error("profile transport returned invalid bulk viscosity")
    return merge(result,(sampling_diagnostics=diagnostics,sampling_profile=effective_config(profile)))
end

end # module TransportSamplingProfiles
