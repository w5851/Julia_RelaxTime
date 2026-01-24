# 收敛性分析：RelaxationTime 模块的数值积分参数回归测试
#
# 本脚本用途：在若干代表点上对影响较大的数值参数做系统化扫描并检查收敛性。
# 具体测试项：
#  - `tau_p_nodes`（平均散射率中对动量 p 的 Gauss-Legendre 节点数）对 τ、η、σ 的影响；
#  - `tau_n_sigma_points`（计算 σ(s) 时 t 积分点数）对上述量的影响；
#  -（可扩展）角度/φ 节点、σ 预计算网格大小等。
#
# 运行：
#   julia --project=. tests/analysis/convergence/run_convergence.jl
#
# 可选环境变量：
#   CONV_TOL=0.05                 # 通用相对收敛容忍度（默认 0.08）
#   CONV_TOL_ANGLEPHI=0.12        # 角度/φ 扫掠容忍度（默认 0.12）
#   CONV_TOL_SIGMA_GRID=0.20      # sigma_grid_n 扫掠容忍度（默认 0.20）
#   CONV_VERBOSE=1                # 输出所有原始结果

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, gauleg

const RT_ASR = TransportWorkflow.RelaxationTime.AverageScatteringRate
const RT_TCS = TransportWorkflow.RelaxationTime.TotalCrossSection
const REQUIRED_PROCESSES = TransportWorkflow.RelaxationTime.REQUIRED_PROCESSES

const DEFAULT_P_NUM = 12
const DEFAULT_T_NUM = 6
const DEFAULT_MAX_ITER = 40

struct ConvCase
    name::String
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    sigma_grid_n::Int
    integration_mode::Symbol
    tr_p_nodes::Int
    tr_p_max_fm::Float64
end

function build_K_data(T_fm::Float64, mu_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    nodes = DEFAULT_MOMENTUM_NODES
    weights = DEFAULT_MOMENTUM_WEIGHTS
    A_u = A(masses.u, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    A_s = A(masses.s, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function integration_grids(mode::Symbol, p_nodes::Int)
    if mode == :finite_15
        pg, pw = gauleg(0.0, 15.0, p_nodes)
        return (pg, pw, Λ_inv_fm)
    elseif mode == :finite_lambda
        pg, pw = gauleg(0.0, Λ_inv_fm, p_nodes)
        return (pg, pw, Λ_inv_fm)
    else
        return (nothing, nothing, nothing)
    end
end

function safe_total_cross_section(process::Symbol, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int, max_tries::Int=4)
    s_try = s
    last_err = nothing
    for _ in 1:max_tries
        try
            σ = RT_TCS.total_cross_section(process, s_try, quark_params, thermo_params, K_coeffs; n_points=n_points)
            isfinite(σ) && return σ
        catch err
            last_err = err
        end
        s_try = s_try * (1.0 + 1e-6) + 1e-10
    end
    @warn "failed to compute sigma; returning NaN" process=process s=s last_error=last_err
    return NaN
end

function build_sigma_caches(processes::Tuple, quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_sigma_points::Int, sigma_grid_n::Int)
    caches = Dict{Symbol,RT_ASR.CrossSectionCache}()
    for process in processes
        s_grid = RT_ASR.design_w0cdf_s_grid(process, quark_params, thermo_params;
            N=sigma_grid_n, p_cutoff=Λ_inv_fm)

        cache = RT_ASR.CrossSectionCache(process)
        n_ok, n_bad = 0, 0
        for s in s_grid
            σ = safe_total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
            if !isfinite(σ)
                n_bad += 1
                continue
            end
            RT_ASR.insert_sigma!(cache, s, σ)
            n_ok += 1
        end
        n_bad > 0 && @warn "sigma grid had non-finite points" process=process n_ok=n_ok n_bad=n_bad
        n_ok >= 2 || error("sigma cache has too few valid points for $process (n_ok=$n_ok)")
        caches[process] = cache
    end
    return caches
end

function solve_point(T_mev::Float64, muB_mev::Float64, xi::Float64, c::ConvCase)
    T_fm = T_mev / ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / ħc_MeV_fm

    base = TransportWorkflow.PNJL.solve(TransportWorkflow.PNJL.FixedMu(), T_fm, muq_fm;
        xi=xi,
        p_num=DEFAULT_P_NUM,
        t_num=DEFAULT_T_NUM,
        seed_strategy=TransportWorkflow.PNJL.MultiSeed(),
        iterations=DEFAULT_MAX_ITER,
    )

    base.converged || error("PNJL equilibrium did not converge at T=$(T_mev), muB=$(muB_mev), xi=$(xi)")

    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

    ktmp = build_K_data(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)
    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi))
    quark_params = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)), A=ktmp.A_vals)

    cs_caches = build_sigma_caches(REQUIRED_PROCESSES, quark_params, thermo_params, ktmp.K_coeffs;
        n_sigma_points=c.tau_n_sigma_points, sigma_grid_n=c.sigma_grid_n)

    p_grid, p_w, sigma_cutoff = integration_grids(c.integration_mode, c.tau_p_nodes)
    cos_grid, cos_w = gauleg(-1.0, 1.0, c.tau_angle_nodes)
    phi_grid, phi_w = gauleg(0.0, 2 * pi, c.tau_phi_nodes)

    res = solve_gap_and_transport(
        T_fm,
        muq_fm;
        xi=xi,
        equilibrium=base,
        compute_tau=true,
        K_coeffs=ktmp.K_coeffs,
        tau=nothing,
        compute_bulk=false,
        p_num=DEFAULT_P_NUM,
        t_num=DEFAULT_T_NUM,
        seed_state=Vector(base.solution),
        solver_kwargs=(iterations=DEFAULT_MAX_ITER,),
        tau_kwargs=(
            p_nodes=c.tau_p_nodes,
            angle_nodes=c.tau_angle_nodes,
            phi_nodes=c.tau_phi_nodes,
            n_sigma_points=c.tau_n_sigma_points,
            cs_caches=cs_caches,
            p_grid=p_grid,
            p_w=p_w,
            cos_grid=cos_grid,
            cos_w=cos_w,
            phi_grid=phi_grid,
            phi_w=phi_w,
            sigma_cutoff=sigma_cutoff,
        ),
        transport_config=TransportWorkflow.TransportIntegrationConfig(p_nodes=c.tr_p_nodes, p_max=c.tr_p_max_fm),
    )

    eq = res.equilibrium
    tr = res.transport
    tau = res.tau
    s_fm3inv = eq.entropy

    eta_over_s = (isfinite(tr.eta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.eta / s_fm3inv) : NaN

    return Dict(
        :tau_u => tau.u,
        :tau_s => tau.s,
        :eta => tr.eta,
        :sigma => tr.sigma,
        :eta_over_s => eta_over_s,
    )
end

function rel_diff(a::Float64, b::Float64)
    if !isfinite(a) || !isfinite(b)
        return NaN
    end
    denom = max(abs(b), 1e-12)
    return abs(a - b) / denom
end

function run_sweep(label::String, T_mev::Float64, muB_mev::Float64, xi::Float64, cases::Vector{ConvCase})
    println("\n== Sweep: $label (T=$(T_mev) MeV, muB=$(muB_mev) MeV, xi=$(xi)) ==")

    results = Dict{String,Dict{Symbol,Float64}}()
    for c in cases
        println("  running: $(c.name)")
        results[c.name] = solve_point(T_mev, muB_mev, xi, c)
    end

    ref = results[cases[end].name]
    metrics = [:tau_u, :tau_s, :eta, :sigma, :eta_over_s]

    for c in cases
        vals = results[c.name]
        diffs = Dict(m => rel_diff(vals[m], ref[m]) for m in metrics)
        max_diff = maximum(values(diffs))
        println("  $(c.name): max_rel_diff=$(round(max_diff, sigdigits=4))")
        for m in metrics
            println("    $(m): value=$(vals[m]), rel_diff=$(round(diffs[m], sigdigits=4))")
        end
    end

    return results
end

const tol = parse(Float64, get(ENV, "CONV_TOL", "0.08"))
const tol_angle_phi = parse(Float64, get(ENV, "CONV_TOL_ANGLEPHI", "0.12"))
const tol_sigma_grid = parse(Float64, get(ENV, "CONV_TOL_SIGMA_GRID", "0.20"))
const verbose = get(ENV, "CONV_VERBOSE", "0") in ("1", "true", "TRUE", "yes", "YES")

@testset "Convergence (RelaxTime)" begin
    # 要测试的代表点列表：
    # - 低 μB、无各向异性： (T=200, μB=0, ξ=0)
    # - 高 μB（近 CSV 中示例），并测试各向异性 ξ=-0.4/+0.4
    sample_points = [
        (200.0, 0.0, 0.0),
        (200.0, 800.0, 0.0),
        (200.0, 800.0, -0.4),
        (200.0, 800.0, 0.4),
    ]

    # 定义要扫的参数阵列
    base = ConvCase("base", 20, 4, 8, 6, 60, :finite_15, 24, 8.0)

    sweep_p = [
        ConvCase("p16", 16, 4, 8, 6, 60, :finite_15, 24, 8.0),
        ConvCase("p20", 20, 4, 8, 6, 60, :finite_15, 24, 8.0),
        ConvCase("p24", 24, 4, 8, 6, 60, :finite_15, 24, 8.0),
    ]

    sweep_sigma = [
        ConvCase("sigma6", 20, 4, 8, 6, 60, :finite_15, 24, 8.0),
        ConvCase("sigma8", 20, 4, 8, 8, 60, :finite_15, 24, 8.0),
        ConvCase("sigma10", 20, 4, 8, 10, 60, :finite_15, 24, 8.0),
    ]

    sweep_angle_phi = [
        ConvCase("ang4_phi8", 20, 4, 8, 6, 60, :finite_15, 24, 8.0),
        ConvCase("ang6_phi12", 20, 6, 12, 6, 60, :finite_15, 24, 8.0),
    ]

    sweep_sigma_grid = [
        ConvCase("grid40", 20, 4, 8, 6, 40, :finite_15, 24, 8.0),
        ConvCase("grid60", 20, 4, 8, 6, 60, :finite_15, 24, 8.0),
        ConvCase("grid80", 20, 4, 8, 6, 80, :finite_15, 24, 8.0),
    ]

    # 为每个代表点运行相同的扫掠，并做断言
    for (T_mev, muB_mev, xi) in sample_points
        results_p = run_sweep("tau_p_nodes", T_mev, muB_mev, xi, sweep_p)
        results_sigma = run_sweep("tau_n_sigma_points", T_mev, muB_mev, xi, sweep_sigma)
        results_angle_phi = run_sweep("tau_angle_phi_nodes", T_mev, muB_mev, xi, sweep_angle_phi)
        results_sigma_grid = run_sweep("sigma_grid_n", T_mev, muB_mev, xi, sweep_sigma_grid)

        # Basic convergence assertions vs highest resolution in each sweep
        ref_p = results_p[sweep_p[end].name]
        for (name, vals) in results_p
            name == sweep_p[end].name && continue
            for key in keys(ref_p)
                @test rel_diff(vals[key], ref_p[key]) <= tol
            end
        end

        ref_sigma = results_sigma[sweep_sigma[end].name]
        for (name, vals) in results_sigma
            name == sweep_sigma[end].name && continue
            for key in keys(ref_sigma)
                @test rel_diff(vals[key], ref_sigma[key]) <= tol
            end
        end

        ref_angle_phi = results_angle_phi[sweep_angle_phi[end].name]
        for (name, vals) in results_angle_phi
            name == sweep_angle_phi[end].name && continue
            for key in keys(ref_angle_phi)
                @test rel_diff(vals[key], ref_angle_phi[key]) <= tol_angle_phi
            end
        end

        ref_sigma_grid = results_sigma_grid[sweep_sigma_grid[end].name]
        for (name, vals) in results_sigma_grid
            name == sweep_sigma_grid[end].name && continue
            for key in keys(ref_sigma_grid)
                @test rel_diff(vals[key], ref_sigma_grid[key]) <= tol_sigma_grid
            end
        end

        if verbose
            println("\nRaw results (tau_p_nodes) for T=$(T_mev), muB=$(muB_mev), xi=$(xi):")
            for c in sweep_p
                println("  $(c.name): $(results_p[c.name])")
            end
            println("\nRaw results (tau_n_sigma_points) for T=$(T_mev), muB=$(muB_mev), xi=$(xi):")
            for c in sweep_sigma
                println("  $(c.name): $(results_sigma[c.name])")
            end
            println("\nRaw results (tau_angle_phi_nodes) for T=$(T_mev), muB=$(muB_mev), xi=$(xi):")
            for c in sweep_angle_phi
                println("  $(c.name): $(results_angle_phi[c.name])")
            end
            println("\nRaw results (sigma_grid_n) for T=$(T_mev), muB=$(muB_mev), xi=$(xi):")
            for c in sweep_sigma_grid
                println("  $(c.name): $(results_sigma_grid[c.name])")
            end
        end
    end
end
