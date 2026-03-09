if !isdefined(Main, :validation_targets_path)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "common", "data_paths.jl"))
end

const RELAXTIME_TAU_VALIDATION_PROJECT_ROOT = VALIDATION_PROJECT_ROOT
const RELAXTIME_TAU_VALIDATION_MODELS_PATH = joinpath(
    RELAXTIME_TAU_VALIDATION_PROJECT_ROOT,
    "src",
    "models",
    "Models.jl",
)

if !isdefined(Main, :Models)
    Base.include(Main, RELAXTIME_TAU_VALIDATION_MODELS_PATH)
end

using .Models

const RELAXTIME_TAU_VALIDATION_TRANSPORT_WORKFLOW = Models.transport_workflow_module()
const RELAXTIME_TAU_VALIDATION_PNJL = Models.pnjl_module()
const RELAXTIME_TAU_VALIDATION_DEFAULT_MOMENTUM_NODES = getproperty(
    getproperty(RELAXTIME_TAU_VALIDATION_PNJL, :Integrals),
    :DEFAULT_MOMENTUM_NODES,
)
const RELAXTIME_TAU_VALIDATION_DEFAULT_MOMENTUM_WEIGHTS = getproperty(
    getproperty(RELAXTIME_TAU_VALIDATION_PNJL, :Integrals),
    :DEFAULT_MOMENTUM_WEIGHTS,
)
const RELAXTIME_TAU_VALIDATION_HADRON_SEED_5 = getproperty(
    RELAXTIME_TAU_VALIDATION_PNJL,
    :HADRON_SEED_5,
)

const RELAXTIME_TAU_VALIDATION_HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm
const RELAXTIME_TAU_VALIDATION_G_FM2 = Main.Constants_PNJL.G_fm2
const RELAXTIME_TAU_VALIDATION_K_FM5 = Main.Constants_PNJL.K_fm5
const RELAXTIME_TAU_VALIDATION_LAMBDA_INV_FM = Main.Constants_PNJL.Λ_inv_fm

const RELAXTIME_TAU_VALIDATION_RT_ASR = Main.AverageScatteringRate
const RELAXTIME_TAU_VALIDATION_RT_RELAX = Main.RelaxationTime
const RELAXTIME_TAU_VALIDATION_N_SIGMA_POINTS = 4
const RELAXTIME_TAU_VALIDATION_P_NODES = 12
const RELAXTIME_TAU_VALIDATION_ANGLE_NODES = 3
const RELAXTIME_TAU_VALIDATION_PHI_NODES = 4
const RELAXTIME_TAU_VALIDATION_SIGMA_GRID_N = 24
const RELAXTIME_TAU_VALIDATION_W0CDF_P_NODES = 8
const RELAXTIME_TAU_VALIDATION_W0CDF_ANGLE_NODES = 3
const RELAXTIME_TAU_VALIDATION_W0CDF_PHI_NODES = 4
const RELAXTIME_TAU_VALIDATION_P_NUM = 12
const RELAXTIME_TAU_VALIDATION_T_NUM = 6
const RELAXTIME_TAU_VALIDATION_MAX_ITER = 40
const RELAXTIME_TAU_VALIDATION_P_GRID, RELAXTIME_TAU_VALIDATION_P_W = Main.GaussLegendre.gauleg(0.0, 15.0, RELAXTIME_TAU_VALIDATION_P_NODES)
const RELAXTIME_TAU_VALIDATION_COS_GRID, RELAXTIME_TAU_VALIDATION_COS_W = Main.GaussLegendre.gauleg(-1.0, 1.0, RELAXTIME_TAU_VALIDATION_ANGLE_NODES)
const RELAXTIME_TAU_VALIDATION_PHI_GRID, RELAXTIME_TAU_VALIDATION_PHI_W = Main.GaussLegendre.gauleg(0.0, 2π, RELAXTIME_TAU_VALIDATION_PHI_NODES)

const RELAXTIME_TAU_VALIDATION_POINT_CACHE = Dict{Tuple{Float64, Float64}, NamedTuple}()

function _load_tau_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 7 || length(cols) == 9 || error("invalid tau validation target row: $line")
        if length(cols) == 9
            push!(rows, (
                record_id=strip(cols[1]),
                series=strip(cols[2]),
                muB_MeV=parse(Float64, strip(cols[3])),
                flavor=strip(cols[4]),
                antiparticle=parse(Bool, lowercase(strip(cols[5]))),
                T_MeV=parse(Float64, strip(cols[6])),
                expected_tau_fm=parse(Float64, strip(cols[7])),
                rtol=parse(Float64, strip(cols[8])),
                source=strip(cols[9]),
            ))
        else
            push!(rows, (
                record_id=strip(cols[1]),
                T_MeV=parse(Float64, strip(cols[2])),
                muB_MeV=parse(Float64, strip(cols[3])),
                flavor=strip(cols[4]),
                antiparticle=parse(Bool, lowercase(strip(cols[5]))),
                expected_tau_fm=parse(Float64, strip(cols[6])),
                rtol=parse(Float64, strip(cols[7])),
                source="legacy_transport_fixedpoints_v1",
            ))
        end
    end
    return rows
end

function _build_tau_validation_K_coeffs(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Phi::Float64, Phibar::Float64)
    A_u = Main.OneLoopIntegrals.A(
        masses.u,
        muq_fm,
        T_fm,
        Phi,
        Phibar,
        RELAXTIME_TAU_VALIDATION_DEFAULT_MOMENTUM_NODES,
        RELAXTIME_TAU_VALIDATION_DEFAULT_MOMENTUM_WEIGHTS,
    )
    A_s = Main.OneLoopIntegrals.A(
        masses.s,
        muq_fm,
        T_fm,
        Phi,
        Phibar,
        RELAXTIME_TAU_VALIDATION_DEFAULT_MOMENTUM_NODES,
        RELAXTIME_TAU_VALIDATION_DEFAULT_MOMENTUM_WEIGHTS,
    )
    G_u = Main.EffectiveCouplings.calculate_G_from_A(A_u, masses.u)
    G_s = Main.EffectiveCouplings.calculate_G_from_A(A_s, masses.s)
    return (
        K_coeffs=Main.EffectiveCouplings.calculate_effective_couplings(
            RELAXTIME_TAU_VALIDATION_G_FM2,
            RELAXTIME_TAU_VALIDATION_K_FM5,
            G_u,
            G_s,
        ),
        A_vals=(u=A_u, d=A_u, s=A_s),
    )
end

function _tau_validation_kwargs()
    return (
        p_nodes=RELAXTIME_TAU_VALIDATION_P_NODES,
        angle_nodes=RELAXTIME_TAU_VALIDATION_ANGLE_NODES,
        phi_nodes=RELAXTIME_TAU_VALIDATION_PHI_NODES,
        n_sigma_points=RELAXTIME_TAU_VALIDATION_N_SIGMA_POINTS,
        p_grid=RELAXTIME_TAU_VALIDATION_P_GRID,
        p_w=RELAXTIME_TAU_VALIDATION_P_W,
        cos_grid=RELAXTIME_TAU_VALIDATION_COS_GRID,
        cos_w=RELAXTIME_TAU_VALIDATION_COS_W,
        phi_grid=RELAXTIME_TAU_VALIDATION_PHI_GRID,
        phi_w=RELAXTIME_TAU_VALIDATION_PHI_W,
        sigma_cutoff=RELAXTIME_TAU_VALIDATION_LAMBDA_INV_FM,
    )
end

function _build_tau_validation_cs_caches(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple,
)
    caches = Dict{Symbol, RELAXTIME_TAU_VALIDATION_RT_ASR.CrossSectionCache}()
    for process in RELAXTIME_TAU_VALIDATION_RT_RELAX.REQUIRED_PROCESSES
        caches[process] = RELAXTIME_TAU_VALIDATION_RT_ASR.build_w0cdf_pchip_cache(
            process,
            quark_params,
            thermo_params,
            K_coeffs;
            N=RELAXTIME_TAU_VALIDATION_SIGMA_GRID_N,
            design_p_nodes=RELAXTIME_TAU_VALIDATION_W0CDF_P_NODES,
            design_angle_nodes=RELAXTIME_TAU_VALIDATION_W0CDF_ANGLE_NODES,
            design_phi_nodes=RELAXTIME_TAU_VALIDATION_W0CDF_PHI_NODES,
            p_cutoff=RELAXTIME_TAU_VALIDATION_LAMBDA_INV_FM,
            n_sigma_points=RELAXTIME_TAU_VALIDATION_N_SIGMA_POINTS,
        )
    end
    return caches
end

function _solve_tau_validation_equilibrium(T_fm::Float64, muq_fm::Float64, xi::Float64)
    try
        return RELAXTIME_TAU_VALIDATION_TRANSPORT_WORKFLOW.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:models,
            p_num=RELAXTIME_TAU_VALIDATION_P_NUM,
            t_num=RELAXTIME_TAU_VALIDATION_T_NUM,
            seed_state=RELAXTIME_TAU_VALIDATION_HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )
    catch
        return RELAXTIME_TAU_VALIDATION_TRANSPORT_WORKFLOW.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:legacy,
            p_num=RELAXTIME_TAU_VALIDATION_P_NUM,
            t_num=RELAXTIME_TAU_VALIDATION_T_NUM,
            seed_state=nothing,
            solver_kwargs=(iterations=RELAXTIME_TAU_VALIDATION_MAX_ITER,),
        )
    end
end

function _compute_relaxtime_tau_point(T_MeV::Float64, muB_MeV::Float64)
    key = (T_MeV, muB_MeV)
    haskey(RELAXTIME_TAU_VALIDATION_POINT_CACHE, key) && return RELAXTIME_TAU_VALIDATION_POINT_CACHE[key]

    T_fm = T_MeV / RELAXTIME_TAU_VALIDATION_HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / RELAXTIME_TAU_VALIDATION_HBARC_MEV_FM
    xi = 0.0

    base = _solve_tau_validation_equilibrium(T_fm, muq_fm, xi)
    Bool(base.converged) || error("equilibrium solve failed at T=$(T_MeV) MeV, muB=$(muB_MeV) MeV")

    Phi = Float64(base.x_state[4])
    Phibar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

    ktmp = _build_tau_validation_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Phi, Phibar)
    quark_params = (
        m=masses,
        μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)),
        A=ktmp.A_vals,
    )
    thermo_params = (T=Float64(T_fm), Φ=Phi, Φbar=Phibar, ξ=xi)
    cs_caches = _build_tau_validation_cs_caches(quark_params, thermo_params, ktmp.K_coeffs)

    tau_res = RELAXTIME_TAU_VALIDATION_TRANSPORT_WORKFLOW.solve_gap_and_transport(
        T_fm,
        muq_fm;
        xi=xi,
        equilibrium=base,
        compute_tau=true,
        K_coeffs=ktmp.K_coeffs,
        tau=nothing,
        compute_bulk=false,
        p_num=RELAXTIME_TAU_VALIDATION_P_NUM,
        t_num=RELAXTIME_TAU_VALIDATION_T_NUM,
        seed_state=Vector(base.x_state),
        solver_kwargs=(iterations=RELAXTIME_TAU_VALIDATION_MAX_ITER,),
        tau_kwargs=(; _tau_validation_kwargs()..., cs_caches=cs_caches),
    )

    RELAXTIME_TAU_VALIDATION_POINT_CACHE[key] = tau_res.tau
    return tau_res.tau
end

function _tau_species_key(row)
    if row.flavor == "u"
        return row.antiparticle ? :ubar : :u
    elseif row.flavor == "s"
        return row.antiparticle ? :sbar : :s
    end
    error("unsupported flavor in tau validation target: $(row.flavor)")
end