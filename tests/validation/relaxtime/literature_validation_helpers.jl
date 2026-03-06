const RELAXTIME_LITERATURE_VALIDATION_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const RELAXTIME_LITERATURE_VALIDATION_MODELS_PATH = joinpath(
    RELAXTIME_LITERATURE_VALIDATION_PROJECT_ROOT,
    "src",
    "models",
    "Models.jl",
)
const RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW_PATH = joinpath(
    RELAXTIME_LITERATURE_VALIDATION_PROJECT_ROOT,
    "src",
    "models",
    "workflows",
    "TransportWorkflow.jl",
)

if !isdefined(Main, :Models)
    Base.include(Main, RELAXTIME_LITERATURE_VALIDATION_MODELS_PATH)
end

using .Models

const RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW = Models.transport_workflow_module()
const RELAXTIME_LITERATURE_VALIDATION_PNJL = Models.pnjl_module()
const RELAXTIME_LITERATURE_VALIDATION_MODEL = Models.create_model(:PNJL)
const RELAXTIME_LITERATURE_VALIDATION_DEFAULT_MOMENTUM_NODES = getproperty(
    getproperty(RELAXTIME_LITERATURE_VALIDATION_PNJL, :Integrals),
    :DEFAULT_MOMENTUM_NODES,
)
const RELAXTIME_LITERATURE_VALIDATION_DEFAULT_MOMENTUM_WEIGHTS = getproperty(
    getproperty(RELAXTIME_LITERATURE_VALIDATION_PNJL, :Integrals),
    :DEFAULT_MOMENTUM_WEIGHTS,
)
const RELAXTIME_LITERATURE_VALIDATION_HADRON_SEED_5 = getproperty(
    RELAXTIME_LITERATURE_VALIDATION_PNJL,
    :HADRON_SEED_5,
)

const RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm
const RELAXTIME_LITERATURE_VALIDATION_G_FM2 = Main.Constants_PNJL.G_fm2
const RELAXTIME_LITERATURE_VALIDATION_K_FM5 = Main.Constants_PNJL.K_fm5
const RELAXTIME_LITERATURE_VALIDATION_MB_PER_FM2 = 10.0
const RELAXTIME_LITERATURE_VALIDATION_VALIDATION_P_NUM = 12
const RELAXTIME_LITERATURE_VALIDATION_VALIDATION_T_NUM = 6
const RELAXTIME_LITERATURE_VALIDATION_VALIDATION_MAX_ITER = 40
const RELAXTIME_LITERATURE_VALIDATION_SIGMA_N_POINTS = 64

const RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_CACHE = Dict{Tuple{Float64, Float64}, NamedTuple}()
const RELAXTIME_LITERATURE_VALIDATION_SIGMA_PARAM_CACHE = Dict{Tuple{Float64, Float64}, NamedTuple}()
const RELAXTIME_LITERATURE_VALIDATION_SIGMA_VALUE_CACHE = Dict{Tuple{Float64, Float64, Symbol, Float64}, Float64}()

function _solve_relaxtime_literature_validation_equilibrium(T_fm::Float64, muq_fm::Float64, xi::Float64)
    try
        return RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:models,
            p_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_P_NUM,
            t_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_T_NUM,
            seed_state=RELAXTIME_LITERATURE_VALIDATION_HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )
    catch
        return RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:legacy,
            p_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_P_NUM,
            t_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_T_NUM,
            seed_state=nothing,
            solver_kwargs=(iterations=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_MAX_ITER,),
        )
    end
end

function _build_relaxtime_literature_validation_K_coeffs(
    T_fm::Float64,
    muq_fm::Float64,
    masses::NamedTuple,
    Phi::Float64,
    Phibar::Float64,
)
    A_u = Main.OneLoopIntegrals.A(
        masses.u,
        muq_fm,
        T_fm,
        Phi,
        Phibar,
        RELAXTIME_LITERATURE_VALIDATION_DEFAULT_MOMENTUM_NODES,
        RELAXTIME_LITERATURE_VALIDATION_DEFAULT_MOMENTUM_WEIGHTS,
    )
    A_s = Main.OneLoopIntegrals.A(
        masses.s,
        muq_fm,
        T_fm,
        Phi,
        Phibar,
        RELAXTIME_LITERATURE_VALIDATION_DEFAULT_MOMENTUM_NODES,
        RELAXTIME_LITERATURE_VALIDATION_DEFAULT_MOMENTUM_WEIGHTS,
    )
    G_u = Main.EffectiveCouplings.calculate_G_from_A(A_u, masses.u)
    G_s = Main.EffectiveCouplings.calculate_G_from_A(A_s, masses.s)
    return (
        K_coeffs=Main.EffectiveCouplings.calculate_effective_couplings(
            RELAXTIME_LITERATURE_VALIDATION_G_FM2,
            RELAXTIME_LITERATURE_VALIDATION_K_FM5,
            G_u,
            G_s,
        ),
        A_vals=(u=A_u, d=A_u, s=A_s),
    )
end

function _compute_relaxtime_literature_transport_point(T_MeV::Float64, muB_MeV::Float64)
    key = (T_MeV, muB_MeV)
    haskey(RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_CACHE, key) &&
        return RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_CACHE[key]

    T_fm = T_MeV / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    xi = 0.0

    equilibrium = _solve_relaxtime_literature_validation_equilibrium(T_fm, muq_fm, xi)
    Bool(equilibrium.converged) || error("equilibrium solve failed at T=$(T_MeV) MeV, muB=$(muB_MeV) MeV")

    Phi = Float64(equilibrium.x_state[4])
    Phibar = Float64(equilibrium.x_state[5])
    masses = (
        u=Float64(equilibrium.masses[1]),
        d=Float64(equilibrium.masses[2]),
        s=Float64(equilibrium.masses[3]),
    )

    ktmp = _build_relaxtime_literature_validation_K_coeffs(T_fm, muq_fm, masses, Phi, Phibar)

    result = RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW.solve_gap_and_transport(
        T_fm,
        muq_fm;
        xi=xi,
        equilibrium=equilibrium,
        compute_tau=true,
        K_coeffs=ktmp.K_coeffs,
        compute_bulk=false,
        p_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_P_NUM,
        t_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_T_NUM,
        seed_state=Vector(equilibrium.x_state),
        solver_kwargs=(iterations=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_MAX_ITER,),
    )

    _, _, s_fm3inv, _ = Models.model_thermo(
        RELAXTIME_LITERATURE_VALIDATION_MODEL,
        equilibrium.x_state,
        equilibrium.mu_vec,
        T_fm;
        p_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_P_NUM,
        t_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_T_NUM,
        xi=xi,
    )
    eta_over_s = result.transport.eta / s_fm3inv

    RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_CACHE[key] = (
        eta_over_s=eta_over_s,
        transport=result.transport,
        equilibrium=equilibrium,
        masses=masses,
        K_coeffs=ktmp.K_coeffs,
    )
    return RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_CACHE[key]
end

function _build_relaxtime_literature_sigma_params(T_MeV::Float64, muB_MeV::Float64)
    key = (T_MeV, muB_MeV)
    haskey(RELAXTIME_LITERATURE_VALIDATION_SIGMA_PARAM_CACHE, key) &&
        return RELAXTIME_LITERATURE_VALIDATION_SIGMA_PARAM_CACHE[key]

    T_fm = T_MeV / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    xi = 0.0

    equilibrium = _solve_relaxtime_literature_validation_equilibrium(T_fm, muq_fm, xi)
    Bool(equilibrium.converged) || error("equilibrium solve failed at T=$(T_MeV) MeV, muB=$(muB_MeV) MeV")

    Phi = Float64(equilibrium.x_state[4])
    Phibar = Float64(equilibrium.x_state[5])
    masses = (
        u=Float64(equilibrium.masses[1]),
        d=Float64(equilibrium.masses[2]),
        s=Float64(equilibrium.masses[3]),
    )

    ktmp = _build_relaxtime_literature_validation_K_coeffs(T_fm, muq_fm, masses, Phi, Phibar)
    quark_params = (
        m=masses,
        μ=(u=muq_fm, d=muq_fm, s=muq_fm),
        A=ktmp.A_vals,
    )
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)

    RELAXTIME_LITERATURE_VALIDATION_SIGMA_PARAM_CACHE[key] = (
        quark_params=quark_params,
        thermo_params=thermo_params,
        K_coeffs=ktmp.K_coeffs,
    )
    return RELAXTIME_LITERATURE_VALIDATION_SIGMA_PARAM_CACHE[key]
end

function _compute_relaxtime_literature_sigma_mb(
    T_MeV::Float64,
    muB_MeV::Float64,
    process::Symbol,
    sqrt_s_GeV::Float64,
)
    key = (T_MeV, muB_MeV, process, sqrt_s_GeV)
    haskey(RELAXTIME_LITERATURE_VALIDATION_SIGMA_VALUE_CACHE, key) &&
        return RELAXTIME_LITERATURE_VALIDATION_SIGMA_VALUE_CACHE[key]

    params = _build_relaxtime_literature_sigma_params(T_MeV, muB_MeV)
    sqrt_s_fm = (sqrt_s_GeV * 1000.0) / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    s_fm2 = sqrt_s_fm^2
    sigma_fm2 = Main.TotalCrossSection.total_cross_section(
        process,
        s_fm2,
        params.quark_params,
        params.thermo_params,
        params.K_coeffs;
        n_points=RELAXTIME_LITERATURE_VALIDATION_SIGMA_N_POINTS,
    )
    sigma_mb = sigma_fm2 * RELAXTIME_LITERATURE_VALIDATION_MB_PER_FM2

    RELAXTIME_LITERATURE_VALIDATION_SIGMA_VALUE_CACHE[key] = sigma_mb
    return sigma_mb
end