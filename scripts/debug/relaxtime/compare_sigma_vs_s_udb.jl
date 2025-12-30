#!/usr/bin/env julia

using Printf
using StaticArrays

function find_project_root(start_dir::AbstractString)
    dir = abspath(start_dir)
    while true
        if isfile(joinpath(dir, "Project.toml"))
            return dir
        end
        parent = dirname(dir)
        parent == dir && error("Could not find Project.toml from: $start_dir")
        dir = parent
    end
end

const PROJECT_ROOT = find_project_root(@__DIR__)

push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src", "relaxtime"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section, DEFAULT_T_INTEGRAL_POINTS

@inline relerr(a::Float64, b::Float64) = abs(a - b) / max(1e-12, abs(b))

@inline function _parse_numbers(line::AbstractString)
    matches = eachmatch(r"(?<![A-Za-z0-9_])[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
    return [parse(Float64, m.match) for m in matches]
end

function read_cpp_sigma_vs_s(path::AbstractString)
    rows = NamedTuple{(:sqrt_s_MeV, :sigma_cpp), Tuple{Float64, Float64}}[]
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        startswith(s, "#") && continue
        nums = _parse_numbers(s)
        length(nums) < 3 && continue
        push!(rows, (sqrt_s_MeV=nums[1], sigma_cpp=nums[3]))
    end
    isempty(rows) && error("No data rows parsed from: $path")
    return rows
end

function build_params(; T_MeV::Float64, muB_MeV::Float64, xi::Float64)
    muq_MeV = muB_MeV / 3.0
    T = T_MeV / ħc_MeV_fm
    μ = muq_MeV / ħc_MeV_fm

    gap_res = PNJL.solve(PNJL.FixedMu(), T, Float64(μ); xi=Float64(xi))
    gap_res.converged || error("Gap solver did not converge (T=$T_MeV MeV, μq=$muq_MeV MeV, xi=$xi)")

    x = gap_res.solution
    ϕ = SVector{3, Float64}(x[1], x[2], x[3])
    Φ = Float64(x[4])
    Φbar = Float64(x[5])
    m_u, m_d, m_s = gap_res.masses[1], gap_res.masses[2], gap_res.masses[3]

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS
    A_u = A(m_u, μ, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ, T, Φ, Φbar, nodes_p, weights_p)
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ, d=μ, s=μ), A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=xi)
    return (quark_params=quark_params, thermo_params=thermo_params, Kc=Kc,
        m_u=m_u, m_d=m_d, m_s=m_s, Phi=Φ, Phibar=Φbar, muq_MeV=muq_MeV)
end

function summarize_errors(label::AbstractString, errs::Vector{Float64}, sqrt_s_MeV::Vector{Float64}, sigma_cpp::Vector{Float64}, sigma_j::Vector{Float64})
    max_err = maximum(errs)
    i = argmax(errs)
    @printf("%s: max relerr = %.3e at sqrt_s=%.3f MeV  (cpp=%.6g julia=%.6g)\n",
        label, max_err, sqrt_s_MeV[i], sigma_cpp[i], sigma_j[i])
end

function main()
    cpp_path = normpath(joinpath(PROJECT_ROOT, "..", "20250413备份", "results", "debug_sigma_vs_s_udb.txt"))
    rows = read_cpp_sigma_vs_s(cpp_path)

    T_MeV = 150.0
    muB_MeV = 800.0
    xi = 0.0
    process = :udbar_to_udbar

    params = build_params(T_MeV=T_MeV, muB_MeV=muB_MeV, xi=xi)
    @printf("Julia gap: muq=%.3f MeV  m_u=%.6f MeV  Phi=%.6f  Phibar=%.6f\n",
        params.muq_MeV, params.m_u * ħc_MeV_fm, params.Phi, params.Phibar)

    sqrt_s_MeV_all = [r.sqrt_s_MeV for r in rows]
    s_vals_all = [(r.sqrt_s_MeV / ħc_MeV_fm)^2 for r in rows]  # fm^-2
    sigma_cpp_all = [r.sigma_cpp for r in rows]

    n_default = DEFAULT_T_INTEGRAL_POINTS
    n_ref = parse(Int, get(ENV, "N_REF", "64"))

    sample_count = parse(Int, get(ENV, "SAMPLE_COUNT", "12"))
    full_scan = lowercase(get(ENV, "FULL_SCAN", "0")) in ("1", "true", "yes")

    idxs = full_scan ? collect(eachindex(rows)) : unique(round.(Int, range(1, length(rows); length=sample_count)))
    sqrt_s_MeV = sqrt_s_MeV_all[idxs]
    s_vals = s_vals_all[idxs]
    sigma_cpp = sigma_cpp_all[idxs]

    sigma_j_default = Vector{Float64}(undef, length(idxs))
    sigma_j_ref = Vector{Float64}(undef, length(idxs))

    for (j, s) in enumerate(s_vals)
        sigma_j_default[j] = total_cross_section(process, s, params.quark_params, params.thermo_params, params.Kc; n_points=n_default)
        sigma_j_ref[j] = total_cross_section(process, s, params.quark_params, params.thermo_params, params.Kc; n_points=n_ref)
        if full_scan && (j % 10 == 0)
            @printf("progress: %d/%d\n", j, length(idxs))
        end
    end

    errs_default = [relerr(sigma_j_default[i], sigma_cpp[i]) for i in eachindex(sigma_cpp)]
    errs_ref = [relerr(sigma_j_ref[i], sigma_cpp[i]) for i in eachindex(sigma_cpp)]

    println("=== Compare σ(s): Julia vs C++ (udb, T=150, muB=800) ===")
    if full_scan
        @printf("Points: %d (FULL_SCAN=1)  (C++ grid from %s)\n", length(rows), cpp_path)
    else
        @printf("Points: %d/%d (sampled; set FULL_SCAN=1 for all)  (C++ grid from %s)\n", length(idxs), length(rows), cpp_path)
    end
    summarize_errors(@sprintf("Julia n_points=%d (DEFAULT)", n_default), errs_default, sqrt_s_MeV, sigma_cpp, sigma_j_default)
    summarize_errors(@sprintf("Julia n_points=%d", n_ref), errs_ref, sqrt_s_MeV, sigma_cpp, sigma_j_ref)

    abs_err_default = maximum(abs.(sigma_j_default .- sigma_cpp))
    abs_err_ref = maximum(abs.(sigma_j_ref .- sigma_cpp))
    @printf("Max abs err (DEFAULT): %.6g\n", abs_err_default)
    @printf("Max abs err (n_points=%d): %.6g\n", n_ref, abs_err_ref)

    # Trimmed relerr: ignore points too close to threshold and where σ is tiny.
    # This is useful because (a) Julia/C++ gap solvers may produce slightly different m_u,
    # shifting the exact threshold; (b) relative error over-emphasizes extremely small σ.
    sqrt_s_lo = parse(Float64, get(ENV, "TRIM_SQRT_S_MIN_MEV", "250"))
    sqrt_s_hi = parse(Float64, get(ENV, "TRIM_SQRT_S_MAX_MEV", "1100"))
    sigma_min = parse(Float64, get(ENV, "TRIM_SIGMA_MIN", "0.05"))
    keep = [i for i in eachindex(sigma_cpp) if (sqrt_s_MeV[i] >= sqrt_s_lo && sqrt_s_MeV[i] <= sqrt_s_hi && sigma_cpp[i] >= sigma_min)]
    if !isempty(keep)
        max_rel_trim_def = maximum(errs_default[keep])
        max_rel_trim_ref = maximum(errs_ref[keep])
        @printf("Trimmed max relerr (sqrt_s in [%.0f,%.0f] MeV and sigma_cpp>=%.3g):\n", sqrt_s_lo, sqrt_s_hi, sigma_min)
        @printf("  DEFAULT n_points=%d: %.3e\n", n_default, max_rel_trim_def)
        @printf("  n_points=%d: %.3e\n", n_ref, max_rel_trim_ref)
    else
        @printf("Trimmed window empty; adjust TRIM_* env vars if needed.\n")
    end

    println("\nSample rows (first 8):")
    @printf("%10s %12s %12s %12s %12s\n", "sqrt_s(MeV)", "sigma_cpp", "sigma_j_def", "sigma_j_64", "relerr64")
    for i in 1:min(8, length(idxs))
        @printf("%10.3f %12.6g %12.6g %12.6g %12.3e\n",
            sqrt_s_MeV[i], sigma_cpp[i], sigma_j_default[i], sigma_j_ref[i], errs_ref[i])
    end
end

main()
