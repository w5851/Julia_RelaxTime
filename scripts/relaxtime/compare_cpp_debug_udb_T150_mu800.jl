#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src", "relaxtime"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "DifferentialCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalPropagator.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "PolarizationCache.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScatteringAmplitude: scattering_amplitude_squared
using .DifferentialCrossSection: differential_cross_section
using .TotalCrossSection: total_cross_section
using .TotalPropagator: calculate_all_propagators_by_channel, calculate_cms_momentum
using .PolarizationCache: polarization_aniso_cached

# -----------------------------
# parsing helpers
# -----------------------------

@inline function _parse_numbers(line::AbstractString)
    # captures floats like -13.0321, 1.23e-4
    # avoid capturing digits embedded in identifiers like k00/k88/k08
    matches = eachmatch(r"(?<![A-Za-z0-9_])[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
    return [parse(Float64, m.match) for m in matches]
end

function parse_cpp_point_summary(path::AbstractString)
    lines = readlines(path)
    T_MeV = nothing
    muB_MeV = nothing
    Phi = nothing
    Phibar = nothing
    m_u_MeV = nothing
    m_d_MeV = nothing
    m_s_MeV = nothing

    for ln in lines
        if startswith(ln, "T(MeV)")
            nums = _parse_numbers(ln)
            T_MeV = nums[end]
        elseif startswith(ln, "muB(MeV)")
            nums = _parse_numbers(ln)
            muB_MeV = nums[end]
        elseif startswith(ln, "Phi1 Phi2")
            nums = _parse_numbers(ln)
            Phi = nums[end-1]
            Phibar = nums[end]
        elseif startswith(ln, "m_u m_d m_s")
            nums = _parse_numbers(ln)
            m_u_MeV, m_d_MeV, m_s_MeV = nums[end-2], nums[end-1], nums[end]
        end
    end

    any(x -> x === nothing, (T_MeV, muB_MeV, Phi, Phibar, m_u_MeV, m_d_MeV, m_s_MeV)) &&
        error("Failed to parse some fields from C++ point summary: $path")

    return (T_MeV=T_MeV, muB_MeV=muB_MeV, Phi=Phi, Phibar=Phibar, m_u_MeV=m_u_MeV, m_d_MeV=m_d_MeV, m_s_MeV=m_s_MeV)
end

function parse_cpp_amp_vs_t(path::AbstractString)
    lines = readlines(path)

    s = nothing
    t_min = nothing
    t_max = nothing
    denom = nothing
    block = nothing

    # header physics
    s_pij_P = nothing; s_pij_S = nothing
    s_PiP = nothing; s_PiS = nothing
    t_PiP = nothing; t_PiS = nothing

    t_k00_P = nothing; t_k88_P = nothing; t_k08_P = nothing; t_detK_P = nothing
    t_k00_S = nothing; t_k88_S = nothing; t_k08_S = nothing; t_detK_S = nothing
    t_Pi_uu_P = nothing; t_Pi_ss_P = nothing
    t_Pi_uu_S = nothing; t_Pi_ss_S = nothing

    pick_t = nothing
    prop_s_P = nothing; prop_s_S = nothing
    prop_t_P = nothing; prop_t_S = nothing

    # t table
    t_rows = Float64[]
    w_rows = Float64[]
    amp2_rows = Float64[]
    wamp2_rows = Float64[]

    in_table = false
    for ln in lines
        if occursin(" pick_sqrt_s", ln) && occursin(" s=", ln)
            # contains s=...
            nums = _parse_numbers(ln)
            # last number is s
            s = nums[end]
        elseif occursin("t_min=", ln) && occursin("t_max=", ln)
            nums = _parse_numbers(ln)
            # expect t_min, t_max
            t_min = nums[end-1]
            t_max = nums[end]
        elseif occursin("block=", ln) && occursin("denom", ln)
            m = match(r"block=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+denom\([^\)]*\)=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", ln)
            m === nothing && error("Failed to parse block/denom line: $ln")
            block = parse(Float64, m.captures[1])
            denom = parse(Float64, m.captures[2])
        elseif occursin("# s-channel light", ln)
            nums = _parse_numbers(ln)
            # expected tail: ... p_ij(P) p_ij(S) PiP_re PiP_im PiS_re PiS_im
            s_pij_P = nums[end-5]
            s_pij_S = nums[end-4]
            s_PiP = (re=nums[end-3], im=nums[end-2])
            s_PiS = (re=nums[end-1], im=nums[end])
        elseif occursin("# t-channel light", ln)
            nums = _parse_numbers(ln)
            t_PiP = (re=nums[end-3], im=nums[end-2])
            t_PiS = (re=nums[end-1], im=nums[end])
        elseif occursin("# t-channel heavy(P)", ln)
            nums = _parse_numbers(ln)
            # expected: k00 k88 k08 detK Pi_uu_re Pi_uu_im Pi_ss_re Pi_ss_im
            t_k00_P, t_k88_P, t_k08_P, t_detK_P = nums[1], nums[2], nums[3], nums[4]
            t_Pi_uu_P = (re=nums[5], im=nums[6])
            t_Pi_ss_P = (re=nums[7], im=nums[8])
        elseif occursin("# t-channel heavy(S)", ln)
            nums = _parse_numbers(ln)
            # expected: k00 k88 k08 detK Pi_uu_re Pi_uu_im Pi_ss_re Pi_ss_im
            t_k00_S, t_k88_S, t_k08_S, t_detK_S = nums[1], nums[2], nums[3], nums[4]
            t_Pi_uu_S = (re=nums[5], im=nums[6])
            t_Pi_ss_S = (re=nums[7], im=nums[8])
        elseif occursin("# pick_t=", ln)
            nums = _parse_numbers(ln)
            pick_t = nums[1]
        elseif occursin("# s-channel k0=", ln) && occursin("P(re,im)=", ln) && occursin("S(re,im)=", ln)
            nums = _parse_numbers(ln)
            # ... k0 k fac P_re P_im S_re S_im
            prop_s_P = (re=nums[end-3], im=nums[end-2])
            prop_s_S = (re=nums[end-1], im=nums[end])
        elseif occursin("# t-channel k0=", ln) && occursin("P(re,im)=", ln) && occursin("S(re,im)=", ln)
            nums = _parse_numbers(ln)
            prop_t_P = (re=nums[end-3], im=nums[end-2])
            prop_t_S = (re=nums[end-1], im=nums[end])
        elseif startswith(ln, "# t((MeV/hc)^2)")
            in_table = true
        elseif in_table
            isempty(strip(ln)) && continue
            startswith(strip(ln), "#") && continue
            nums = _parse_numbers(ln)
            if length(nums) >= 5
                push!(t_rows, nums[1])
                push!(w_rows, nums[3])
                push!(amp2_rows, nums[4])
                push!(wamp2_rows, nums[5])
            end
        end
    end

    any(x -> x === nothing, (s, t_min, t_max, denom, block, pick_t)) &&
        error("Failed to parse required fields from C++ amp file: $path")

    return (
        s=s, t_min=t_min, t_max=t_max, denom=denom, block=block, pick_t=pick_t,
        s_pij_P=s_pij_P, s_pij_S=s_pij_S, s_PiP=s_PiP, s_PiS=s_PiS,
        t_PiP=t_PiP, t_PiS=t_PiS,
        t_k00_P=t_k00_P, t_k88_P=t_k88_P, t_k08_P=t_k08_P, t_detK_P=t_detK_P,
        t_k00_S=t_k00_S, t_k88_S=t_k88_S, t_k08_S=t_k08_S, t_detK_S=t_detK_S,
        t_Pi_uu_P=t_Pi_uu_P, t_Pi_ss_P=t_Pi_ss_P,
        t_Pi_uu_S=t_Pi_uu_S, t_Pi_ss_S=t_Pi_ss_S,
        prop_s_P=prop_s_P, prop_s_S=prop_s_S,
        prop_t_P=prop_t_P, prop_t_S=prop_t_S,
        t_rows=t_rows, w_rows=w_rows, amp2_rows=amp2_rows, wamp2_rows=wamp2_rows
    )
end

@inline relerr(a::Float64, b::Float64) = abs(a - b) / max(1e-12, abs(b))

function main()
    cpp_results_dir = normpath(joinpath(PROJECT_ROOT, "..", "20250413备份", "results"))
    cpp_point = joinpath(cpp_results_dir, "debug_point_summary_udb.txt")
    cpp_amp = joinpath(cpp_results_dir, "debug_amp_vs_t_udb.txt")

    @info "Reading C++ debug files" cpp_point cpp_amp
    point = parse_cpp_point_summary(cpp_point)
    amp = parse_cpp_amp_vs_t(cpp_amp)

    T = point.T_MeV / ħc_MeV_fm
    muq = (point.muB_MeV / 3.0) / ħc_MeV_fm

    Phi = point.Phi
    Phibar = point.Phibar

    m_u = point.m_u_MeV / ħc_MeV_fm
    m_d = point.m_d_MeV / ħc_MeV_fm
    m_s = point.m_s_MeV / ħc_MeV_fm

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS

    A_u = A(m_u, muq, T, Phi, Phibar, nodes_p, weights_p)
    A_s = A(m_s, muq, T, Phi, Phibar, nodes_p, weights_p)

    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=muq, d=muq, s=muq), A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=T, Φ=Phi, Φbar=Phibar, ξ=0.0)

    process = :udbar_to_udbar
    s = amp.s

    println("=== Julia vs C++ compare (udb, T=150, muB=800) ===")
    @printf("Inputs from C++ summary: Phi=%.6f Phibar=%.6f  m_u=%.6fMeV m_s=%.6fMeV\n",
        Phi, Phibar, point.m_u_MeV, point.m_s_MeV)
    @printf("Parsed from C++ amp file: s=%.6f  t_min=%.6f  t_max=%.6f  denom=%.3f  block=%.6f\n",
        amp.s, amp.t_min, amp.t_max, amp.denom, amp.block)

    # denom cross-check
    splus = s - (m_u + m_d)^2
    sminus = s - (m_u - m_d)^2
    denom_julia = 16π * splus * sminus
    @printf("denom check: C++ %.6f vs Julia %.6f (relerr=%.3e)\n", amp.denom, denom_julia, relerr(denom_julia, amp.denom))

    # Pi & propagators at pick_t
    pick_t = amp.pick_t
    cms_s = calculate_cms_momentum(process, s, pick_t, :s, quark_params)
    cms_t = calculate_cms_momentum(process, s, pick_t, :t, quark_params)

    PiP_s_re, PiP_s_im = polarization_aniso_cached(:P, cms_s.k0, cms_s.k, m_u, m_d, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)
    PiS_s_re, PiS_s_im = polarization_aniso_cached(:S, cms_s.k0, cms_s.k, m_u, m_d, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)

    PiP_t_uu_re, PiP_t_uu_im = polarization_aniso_cached(:P, cms_t.k0, cms_t.k, m_u, m_u, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)
    PiS_t_uu_re, PiS_t_uu_im = polarization_aniso_cached(:S, cms_t.k0, cms_t.k, m_u, m_u, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)

    PiP_t_ss_re, PiP_t_ss_im = polarization_aniso_cached(:P, cms_t.k0, cms_t.k, m_s, m_s, muq, muq, T, Phi, Phibar, 0.0, A_s, A_s, 0)
    PiS_t_ss_re, PiS_t_ss_im = polarization_aniso_cached(:S, cms_t.k0, cms_t.k, m_s, m_s, muq, muq, T, Phi, Phibar, 0.0, A_s, A_s, 0)

    props_s = calculate_all_propagators_by_channel(process, cms_s.k0, cms_s.k, quark_params, thermo_params, Kc)
    props_t = calculate_all_propagators_by_channel(process, cms_t.k0, cms_t.k, quark_params, thermo_params, Kc)

    println("\n--- Polarization Π (at C++ pick_t) ---")
    @printf("s ΠP Julia %.6f,%.6f  | C++ %.6f,%.6f\n", PiP_s_re, PiP_s_im, amp.s_PiP.re, amp.s_PiP.im)
    @printf("s ΠS Julia %.6f,%.6f  | C++ %.6f,%.6f\n", PiS_s_re, PiS_s_im, amp.s_PiS.re, amp.s_PiS.im)
    @printf("t(uu) ΠP Julia %.6f,%.6f | C++ %.6f,%.6f\n", PiP_t_uu_re, PiP_t_uu_im, amp.t_PiP.re, amp.t_PiP.im)
    @printf("t(uu) ΠS Julia %.6f,%.6f | C++ %.6f,%.6f\n", PiS_t_uu_re, PiS_t_uu_im, amp.t_PiS.re, amp.t_PiS.im)
    if amp.t_Pi_ss_P !== nothing
        @printf("t(ss) ΠP Julia %.6f,%.6f | C++ %.6f,%.6f\n", PiP_t_ss_re, PiP_t_ss_im, amp.t_Pi_ss_P.re, amp.t_Pi_ss_P.im)
    end
    if amp.t_Pi_ss_S !== nothing
        @printf("t(ss) ΠS Julia %.6f,%.6f | C++ %.6f,%.6f\n", PiS_t_ss_re, PiS_t_ss_im, amp.t_Pi_ss_S.re, amp.t_Pi_ss_S.im)
    end

    println("\n--- Propagators D (at C++ pick_t) ---")
    @printf("s D_P Julia %.6f%+.6fim | C++ %.6f%+.6fim\n", real(props_s.s_P), imag(props_s.s_P), amp.prop_s_P.re, amp.prop_s_P.im)
    @printf("s D_S Julia %.6f%+.6fim | C++ %.6f%+.6fim\n", real(props_s.s_S), imag(props_s.s_S), amp.prop_s_S.re, amp.prop_s_S.im)
    @printf("t D_P Julia %.6f%+.6fim | C++ %.6f%+.6fim\n", real(props_t.t_P), imag(props_t.t_P), amp.prop_t_P.re, amp.prop_t_P.im)
    @printf("t D_S Julia %.6f%+.6fim | C++ %.6f%+.6fim\n", real(props_t.t_S), imag(props_t.t_S), amp.prop_t_S.re, amp.prop_t_S.im)

    # Msq / dsdt table compare
    println("\n--- |M|^2 and dσ/dt over Gauss-Legendre t nodes ---")
    n = length(amp.t_rows)
    @printf("Parsed %d t-nodes from C++ file\n", n)

    max_rel_M = 0.0
    max_rel_dsdt = 0.0
    worst_idx_M = 0
    worst_idx_dsdt = 0

    for i in 1:n
        t = amp.t_rows[i]
        amp2_cpp = amp.amp2_rows[i]
        dsdt_cpp = amp2_cpp / amp.denom

        Msq = scattering_amplitude_squared(process, s, t, quark_params, thermo_params, Kc)
        dsdt = differential_cross_section(splus, sminus, Msq)

        rM = relerr(Msq, amp2_cpp)
        rD = relerr(dsdt, dsdt_cpp)
        if rM > max_rel_M
            max_rel_M = rM
            worst_idx_M = i
        end
        if rD > max_rel_dsdt
            max_rel_dsdt = rD
            worst_idx_dsdt = i
        end
    end

    @printf("max relerr(|M|^2) = %.3e at i=%d (t=%.6f)\n", max_rel_M, worst_idx_M, worst_idx_M == 0 ? NaN : amp.t_rows[worst_idx_M])
    @printf("max relerr(dσ/dt) = %.3e at i=%d (t=%.6f)\n", max_rel_dsdt, worst_idx_dsdt, worst_idx_dsdt == 0 ? NaN : amp.t_rows[worst_idx_dsdt])

    # sigma compare
    sigma_cpp = amp.block / amp.denom * sum(amp.wamp2_rows)
    sigma_julia = total_cross_section(process, s, quark_params, thermo_params, Kc; n_points=n)

    println("\n--- σ_total comparison (using C++ table weights) ---")
    @printf("sigma_cpp_from_table = %.12f\n", sigma_cpp)
    @printf("sigma_julia(n_points=%d) = %.12f\n", n, sigma_julia)
    @printf("relerr(sigma) = %.3e\n", relerr(sigma_julia, sigma_cpp))

    println("\nDone.")
end

main()
