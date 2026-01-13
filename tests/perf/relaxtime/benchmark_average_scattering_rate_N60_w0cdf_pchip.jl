#!/usr/bin/env julia

# Benchmark: AverageScatteringRate with the finalized production method
# - σ(s) cache: w0cdf grid (default N=60) + PCHIP interpolation
# - ω integral: AverageScatteringRate defaults (p=20, angle=4, phi=8)
#
# Run:
#   julia --project=. --eval 'include("tests/perf/relaxtime/benchmark_average_scattering_rate_N60_w0cdf_pchip.jl")'
#
# Optional env knobs (defaults match production):
#   PROCESS=udbar_to_udbar
#   T_MEV=150  MUB_MEV=800  XI=0.8  PHI=0.5  PHIBAR=0.5
#   FAST=1                         (quick sanity run; overrides nodes unless explicitly set)
#   DESIGN_P_NODES=14  DESIGN_ANGLE_NODES=4  DESIGN_PHI_NODES=8
#   OMEGA_P_NODES=20   OMEGA_ANGLE_NODES=4   OMEGA_PHI_NODES=8
#   N_SIGMA_POINTS=<int>           (t-integral points for σ(s))

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings

const ASR = AverageScatteringRate

"""Build minimal parameter bundle required by TotalCrossSection/ASR."""
function build_params(; T_MeV::Float64, muB_MeV::Float64,
    mu_s_MeV::Float64=0.0, m_u_MeV::Float64=300.0, m_d_MeV::Float64=300.0, m_s_MeV::Float64=500.0,
    phi::Float64=0.5, phibar::Float64=0.5, xi::Float64=0.0)

    T = T_MeV / ħc_MeV_fm
    μ_q = (muB_MeV / 3.0) / ħc_MeV_fm
    μ_u = μ_q
    μ_d = μ_q
    μ_s = mu_s_MeV / ħc_MeV_fm

    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    A_u = A(m_u, μ_u, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_d = A(m_d, μ_d, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(m_s, μ_s, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)

    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ_u, d=μ_d, s=μ_s), A=(u=A_u, d=A_d, s=A_s))
    thermo_params = (T=T, Φ=phi, Φbar=phibar, ξ=xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function main()
    process = Symbol(get(ENV, "PROCESS", "udbar_to_udbar"))

    T_MeV = parse(Float64, get(ENV, "T_MEV", "150.0"))
    muB_MeV = parse(Float64, get(ENV, "MUB_MEV", "800.0"))
    phi = parse(Float64, get(ENV, "PHI", "0.5"))
    phibar = parse(Float64, get(ENV, "PHIBAR", "0.5"))
    # Default xi>0 is more representative for anisotropy-sensitive integration
    xi = parse(Float64, get(ENV, "XI", "0.8"))

    fast = get(ENV, "FAST", "0") in ("1", "true", "TRUE", "yes", "YES")

    design_p_nodes = parse(Int, get(ENV, "DESIGN_P_NODES", string(ASR.DEFAULT_W0CDF_P_NODES)))
    design_angle_nodes = parse(Int, get(ENV, "DESIGN_ANGLE_NODES", string(ASR.DEFAULT_W0CDF_ANGLE_NODES)))
    design_phi_nodes = parse(Int, get(ENV, "DESIGN_PHI_NODES", string(ASR.DEFAULT_W0CDF_PHI_NODES)))

    omega_p_nodes = parse(Int, get(ENV, "OMEGA_P_NODES", string(ASR.DEFAULT_P_NODES)))
    omega_angle_nodes = parse(Int, get(ENV, "OMEGA_ANGLE_NODES", string(ASR.DEFAULT_ANGLE_NODES)))
    omega_phi_nodes = parse(Int, get(ENV, "OMEGA_PHI_NODES", string(ASR.DEFAULT_PHI_NODES)))

    n_sigma_points_env = get(ENV, "N_SIGMA_POINTS", "")
    n_sigma_points = isempty(n_sigma_points_env) ? nothing : parse(Int, n_sigma_points_env)

    # FAST mode picks smaller node counts for a quick regression run.
    if fast
        if !haskey(ENV, "DESIGN_P_NODES")
            design_p_nodes = min(design_p_nodes, 6)
        end
        if !haskey(ENV, "DESIGN_ANGLE_NODES")
            design_angle_nodes = min(design_angle_nodes, 2)
        end
        if !haskey(ENV, "DESIGN_PHI_NODES")
            design_phi_nodes = min(design_phi_nodes, 2)
        end
        if !haskey(ENV, "OMEGA_P_NODES")
            omega_p_nodes = min(omega_p_nodes, 6)
        end
        if !haskey(ENV, "OMEGA_ANGLE_NODES")
            omega_angle_nodes = min(omega_angle_nodes, 2)
        end
        if !haskey(ENV, "OMEGA_PHI_NODES")
            omega_phi_nodes = min(omega_phi_nodes, 2)
        end
        if n_sigma_points === nothing
            n_sigma_points = 4
        end
    end

    params = build_params(T_MeV=T_MeV, muB_MeV=muB_MeV, phi=phi, phibar=phibar, xi=xi)

    @printf("process=%s\n", string(process))
    @printf("thermo: T=%.1f MeV  muB=%.1f MeV  xi=%.3f  Phi=%.3f  Phibar=%.3f\n", T_MeV, muB_MeV, xi, phi, phibar)
    @printf("sigma_grid_N=%d\n", ASR.DEFAULT_SIGMA_GRID_N)
    @printf("design nodes (p,angle,phi)=(%d,%d,%d)\n", design_p_nodes, design_angle_nodes, design_phi_nodes)
    @printf("omega nodes  (p,angle,phi)=(%d,%d,%d)\n", omega_p_nodes, omega_angle_nodes, omega_phi_nodes)
    if n_sigma_points === nothing
        @printf("sigma t-integral points: DEFAULT (TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)\n\n")
    else
        @printf("sigma t-integral points: %d\n\n", n_sigma_points)
    end

    # 1) w0cdf grid design (no σ evaluation)
    s_grid = Float64[]
    t_design = @elapsed begin
        s_grid = ASR.design_w0cdf_s_grid(
            process,
            params.quark_params,
            params.thermo_params;
            N=ASR.DEFAULT_SIGMA_GRID_N,
            p_nodes=design_p_nodes,
            angle_nodes=design_angle_nodes,
            phi_nodes=design_phi_nodes,
        )
    end

    # 2) σ(s) precompute on that grid + PCHIP slope build
    cache = ASR.CrossSectionCache(process)
    t_sigma = @elapsed begin
        if n_sigma_points === nothing
            ASR.precompute_cross_section!(cache, s_grid, params.quark_params, params.thermo_params, params.K_coeffs)
        else
            ASR.precompute_cross_section!(cache, s_grid, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_sigma_points)
        end
        # Ensure PCHIP slopes are materialized for fair ω timing.
        ASR._ensure_pchip_slopes!(cache)
    end

    # 3) ω integral (AverageScatteringRate defaults)
    ω = NaN
    t_omega = @elapsed begin
        kwargs = (
            p_nodes=omega_p_nodes,
            angle_nodes=omega_angle_nodes,
            phi_nodes=omega_phi_nodes,
            cs_cache=cache,
        )
        if n_sigma_points === nothing
            ω = ASR.average_scattering_rate(process, params.quark_params, params.thermo_params, params.K_coeffs; kwargs...)
        else
            ω = ASR.average_scattering_rate(process, params.quark_params, params.thermo_params, params.K_coeffs; kwargs..., n_sigma_points=n_sigma_points)
        end
    end

    @printf("w0cdf design:        %.3fs (N=%d)\n", t_design, length(s_grid))
    @printf("sigma precompute:    %.3fs (cache_points=%d)\n", t_sigma, length(cache.s_vals))
    @printf("omega integral:      %.3fs\n", t_omega)
    @printf("total:              %.3fs\n", t_design + t_sigma + t_omega)
    @printf("omega value:         %.6e\n", ω)

    isfinite(ω) || error("ω is not finite; check inputs")
    return nothing
end

main()
