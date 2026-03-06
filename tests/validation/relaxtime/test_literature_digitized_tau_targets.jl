using Test

const RELAXTIME_VALIDATION_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RELAXTIME_VALIDATION_DATA_PATH = joinpath(RELAXTIME_VALIDATION_PROJECT_ROOT, "tests", "validation", "data", "relaxtime_tau_literature_targets_v1.csv")

const RELAXTIME_VALIDATION_MODELS_PATH = joinpath(RELAXTIME_VALIDATION_PROJECT_ROOT, "src", "models", "Models.jl")

if !isdefined(Main, :Models)
    Base.include(Main, RELAXTIME_VALIDATION_MODELS_PATH)
end

using .Models

const TransportWorkflow = Models.transport_workflow_module()
const PNJL = Models.pnjl_module()
const DEFAULT_MOMENTUM_NODES = getproperty(getproperty(PNJL, :Integrals), :DEFAULT_MOMENTUM_NODES)
const DEFAULT_MOMENTUM_WEIGHTS = getproperty(getproperty(PNJL, :Integrals), :DEFAULT_MOMENTUM_WEIGHTS)
const HADRON_SEED_5 = getproperty(PNJL, :HADRON_SEED_5)

const RELAXTIME_VALIDATION_HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm
const G_FM2 = Main.Constants_PNJL.G_fm2
const K_FM5 = Main.Constants_PNJL.K_fm5
const LAMBDA_INV_FM = Main.Constants_PNJL.Λ_inv_fm

const RT_ASR = Main.AverageScatteringRate
const RT_TCS = Main.TotalCrossSection
const N_SIGMA_POINTS = RT_TCS.DEFAULT_T_INTEGRAL_POINTS
const TAU_P_NODES = RT_ASR.DEFAULT_P_NODES
const TAU_ANGLE_NODES = RT_ASR.DEFAULT_ANGLE_NODES
const TAU_PHI_NODES = RT_ASR.DEFAULT_PHI_NODES
const VALIDATION_P_NUM = 12
const VALIDATION_T_NUM = 6
const VALIDATION_MAX_ITER = 40
const P_GRID, P_W = Main.GaussLegendre.gauleg(0.0, 15.0, TAU_P_NODES)
const COS_GRID, COS_W = Main.GaussLegendre.gauleg(-1.0, 1.0, TAU_ANGLE_NODES)
const PHI_GRID, PHI_W = Main.GaussLegendre.gauleg(0.0, 2π, TAU_PHI_NODES)

const POINT_CACHE = Dict{Tuple{Float64, Float64}, NamedTuple}()

function _load_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            series=strip(cols[2]),
            muB_MeV=parse(Float64, strip(cols[3])),
            flavor=strip(cols[4]),
            antiparticle=parse(Bool, lowercase(strip(cols[5]))),
            T_MeV=parse(Float64, strip(cols[6])),
            expected_tau_fm=parse(Float64, strip(cols[7])),
            rtol=parse(Float64, strip(cols[8])),
            source=strip(cols[9]),
        ))
    end
    return rows
end

function _build_K_coeffs(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Phi::Float64, Phibar::Float64)
    A_u = Main.OneLoopIntegrals.A(masses.u, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = Main.OneLoopIntegrals.A(masses.s, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = Main.EffectiveCouplings.calculate_G_from_A(A_u, masses.u)
    G_s = Main.EffectiveCouplings.calculate_G_from_A(A_s, masses.s)
    return (
        K_coeffs=Main.EffectiveCouplings.calculate_effective_couplings(G_FM2, K_FM5, G_u, G_s),
        A_vals=(u=A_u, d=A_u, s=A_s),
    )
end

function _validation_tau_kwargs()
    return (
        p_nodes=TAU_P_NODES,
        angle_nodes=TAU_ANGLE_NODES,
        phi_nodes=TAU_PHI_NODES,
        n_sigma_points=N_SIGMA_POINTS,
        p_grid=P_GRID,
        p_w=P_W,
        cos_grid=COS_GRID,
        cos_w=COS_W,
        phi_grid=PHI_GRID,
        phi_w=PHI_W,
        sigma_cutoff=LAMBDA_INV_FM,
    )
end

function _solve_validation_equilibrium(T_fm::Float64, muq_fm::Float64, xi::Float64)
    try
        return TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:models,
            p_num=VALIDATION_P_NUM,
            t_num=VALIDATION_T_NUM,
            seed_state=HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )
    catch
        return TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:legacy,
            p_num=VALIDATION_P_NUM,
            t_num=VALIDATION_T_NUM,
            seed_state=nothing,
            solver_kwargs=(iterations=VALIDATION_MAX_ITER,),
        )
    end
end

function _compute_tau_point(T_MeV::Float64, muB_MeV::Float64)
    key = (T_MeV, muB_MeV)
    haskey(POINT_CACHE, key) && return POINT_CACHE[key]

    T_fm = T_MeV / RELAXTIME_VALIDATION_HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / RELAXTIME_VALIDATION_HBARC_MEV_FM
    xi = 0.0

    base = _solve_validation_equilibrium(T_fm, muq_fm, xi)
    Bool(base.converged) || error("equilibrium solve failed at T=$(T_MeV) MeV, muB=$(muB_MeV) MeV")

    Phi = Float64(base.x_state[4])
    Phibar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

    ktmp = _build_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Phi, Phibar)

    tau_res = TransportWorkflow.solve_gap_and_transport(
        T_fm,
        muq_fm;
        xi=xi,
        equilibrium=base,
        compute_tau=true,
        K_coeffs=ktmp.K_coeffs,
        tau=nothing,
        compute_bulk=false,
        p_num=VALIDATION_P_NUM,
        t_num=VALIDATION_T_NUM,
        seed_state=Vector(base.x_state),
        solver_kwargs=(iterations=VALIDATION_MAX_ITER,),
        tau_kwargs=_validation_tau_kwargs(),
    )

    POINT_CACHE[key] = tau_res.tau
    return tau_res.tau
end

function _species_key(row)
    if row.flavor == "u"
        return row.antiparticle ? :ubar : :u
    elseif row.flavor == "s"
        return row.antiparticle ? :sbar : :s
    end
    error("unsupported flavor in validation target: $(row.flavor)")
end

@testset "RelaxTime literature digitized tau targets" begin
    targets = _load_targets(RELAXTIME_VALIDATION_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        tau = _compute_tau_point(row.T_MeV, row.muB_MeV)
        species = _species_key(row)
        actual_tau_fm = getproperty(tau, species)

        @test isfinite(actual_tau_fm)
        @test isapprox(actual_tau_fm, row.expected_tau_fm; rtol=row.rtol, atol=0.0)
    end
end