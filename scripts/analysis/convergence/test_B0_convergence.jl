# B0 函数收敛性测试
# 运行：julia --project=. tests/analysis/convergence/test_B0_convergence.jl

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

using .OneLoopIntegrals: B0, EPS_K, EPS_SEGMENT, PV_GAP_REL,
    energy_cutoff, internal_momentum, distribution_value_b0, distribution_integral_b0,
    singularity_k_zero, singularity_k_positive,
    real_integrand_k_zero, real_integrand_k_positive,
    integrate_hybrid_interval, SING_NONE, SING_LEFT, SING_RIGHT, SING_BOTH

function rel_diff(a::Float64, b::Float64)
    denom = max(abs(b), 1e-12)
    return abs(a - b) / denom
end

const TOL = parse(Float64, get(ENV, "B0_CONV_TOL", "5e-4"))

# 选取包含虚部的代表点（参数单位：fm）
const CASES = [
    (λ=3.0, k=0.0, m1=0.04, m2=0.04, μ1=0.0, μ2=0.0, T=1.5, Φ=0.84, Φbar=0.84),
    (λ=3.3, k=0.0, m1=0.04, m2=1.03, μ1=0.0, μ2=0.0, T=1.5, Φ=0.84, Φbar=0.84),
    (λ=2.5, k=0.4, m1=0.30, m2=0.30, μ1=0.0, μ2=0.0, T=1.0, Φ=0.5, Φbar=0.5),
]

function integrate_piecewise_hybrid_ref(integrand_fun, Emin, Emax, intervals; n_smooth::Int, n_sing::Int)
    if isempty(intervals)
        return integrate_hybrid_interval(integrand_fun, Emin, Emax, SING_NONE; n=n_smooth)
    end

    if length(intervals) == 1
        E1, E2 = intervals[1]
        total = 0.0
        if E1 > Emin
            total += integrate_hybrid_interval(integrand_fun, Emin, E1, SING_RIGHT; n=n_sing)
        end
        if E2 > E1
            total += integrate_hybrid_interval(integrand_fun, E1, E2, SING_BOTH; n=n_sing)
        end
        if Emax > E2
            total += integrate_hybrid_interval(integrand_fun, E2, Emax, SING_LEFT; n=n_smooth)
        end
        return total
    elseif length(intervals) == 2
        a1, b1 = intervals[1]
        a2, b2 = intervals[2]
        tol = 1e-12 * max(1.0, Emax - Emin)

        E1 = NaN
        E2 = NaN
        if abs(a1 - Emin) <= tol
            E1 = b1
        elseif abs(a2 - Emin) <= tol
            E1 = b2
        end

        if abs(b1 - Emax) <= tol
            E2 = a1
        elseif abs(b2 - Emax) <= tol
            E2 = a2
        end

        if isfinite(E1) && isfinite(E2)
            total = 0.0
            if E1 > Emin
                total += integrate_hybrid_interval(integrand_fun, Emin, E1, SING_RIGHT; n=n_sing)
            end
            if E2 > E1
                total += integrate_hybrid_interval(integrand_fun, E1, E2, SING_BOTH; n=n_sing)
            end
            if Emax > E2
                total += integrate_hybrid_interval(integrand_fun, E2, Emax, SING_LEFT; n=n_smooth)
            end
            return total
        end
    end

    pts = Float64[Emin]
    for (a, b) in intervals
        push!(pts, a)
        push!(pts, b)
    end
    push!(pts, Emax)
    pts = sort(unique(pts))

    total = 0.0
    for i in 1:(length(pts) - 1)
        a = pts[i]
        b = pts[i + 1]
        if !(isfinite(a) && isfinite(b)) || b <= a
            continue
        end
        sing_pos = (i == 1) ? SING_RIGHT : (i == length(pts) - 1 ? SING_LEFT : SING_BOTH)
        n_nodes = (i == length(pts) - 1) ? n_smooth : n_sing
        total += integrate_hybrid_interval(integrand_fun, a, b, sing_pos; n=n_nodes)
    end
    return total
end

function tilde_B0_k_zero_ref(sign_flag::Symbol, λ::Float64, m::Float64, m_prime::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64; n_smooth::Int, n_sing::Int)
    m_pos = max(m, 0.0)
    m_prime_pos = max(m_prime, 0.0)
    Emin = m_pos
    Emax = energy_cutoff(m_pos)
    denominator_term = (λ ^ 2 + m_pos ^ 2 - m_prime_pos ^ 2) / 2.0
    singularity = singularity_k_zero(λ, Emin, Emax, denominator_term)
    integrand_fun(E) = real_integrand_k_zero(sign_flag, λ, m_pos, denominator_term, μ, T, Φ, Φbar, E)

    imag_part = 0.0
    if isempty(singularity)
        real_part = integrate_hybrid_interval(integrand_fun, Emin, Emax, SING_NONE; n=n_smooth)
        return real_part * 2.0, imag_part / λ
    end

    E0 = singularity[1]

    @inline function numer(E::Float64)
        p = internal_momentum(E, m_pos)
        dist = distribution_value_b0(sign_flag, E, μ, T, Φ, Φbar)
        if !isfinite(p) || !isfinite(dist)
            return 0.0
        end
        val = p * dist
        return isfinite(val) ? val : 0.0
    end

    f0 = numer(E0)
    endpoint_tol = 1e-10 * max(1.0, Emax - Emin)
    if abs(E0 - Emin) <= endpoint_tol || abs(Emax - E0) <= endpoint_tol
        gap = max(PV_GAP_REL * (Emax - Emin), 64 * eps(E0))
        left_end = max(Emin, E0 - gap)
        right_start = min(Emax, E0 + gap)
        pv_integral = 0.0
        if left_end > Emin
            pv_integral += integrate_hybrid_interval(integrand_fun, Emin, left_end, SING_NONE; n=n_smooth)
        end
        if right_start < Emax
            pv_integral += integrate_hybrid_interval(integrand_fun, right_start, Emax, SING_NONE; n=n_smooth)
        end
    else
        @inline function regular_integrand(E::Float64)
            d = E - E0
            if abs(d) < 64 * eps(E0)
                δ = max(1e-7 * (Emax - Emin), 64 * eps(E0))
                El = max(Emin, E0 - δ)
                Er = min(Emax, E0 + δ)
                if Er > El
                    return (numer(Er) - numer(El)) / (Er - El)
                else
                    return 0.0
                end
            end
            return (numer(E) - f0) / d
        end

        regular = integrate_hybrid_interval(regular_integrand, Emin, Emax, SING_NONE; n=n_sing)
        logterm = f0 * (log(abs(Emax - E0)) - log(abs(Emin - E0)))
        pv_integral = (regular + logterm) / λ
    end

    p0 = internal_momentum(E0, m_pos)
    imag_part = 2.0 * π * p0 * distribution_value_b0(sign_flag, E0, μ, T, Φ, Φbar)

    return pv_integral * 2.0, imag_part / λ
end

function tilde_B0_k_positive_ref(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64; n_smooth::Int, n_sing::Int)
    m_pos = max(m, 0.0)
    m_prime_pos = max(m_prime, 0.0)
    Emin = m_pos
    Emax = energy_cutoff(m_pos)
    integrand_fun(E) = real_integrand_k_positive(sign_flag, λ, k, m_pos, m_prime_pos, μ, T, Φ, Φbar, E)
    intervals, _ = singularity_k_positive(λ, k, m_pos, m_prime_pos, Emin, Emax)

    imag_part = 0.0
    real_part = integrate_piecewise_hybrid_ref(integrand_fun, Emin, Emax, intervals; n_smooth=n_smooth, n_sing=n_sing)

    if !isempty(intervals)
        for (E1, E2) in intervals
            imag_part += distribution_integral_b0(sign_flag, E1, E2, μ, T, Φ, Φbar)
        end
        imag_part *= π * sign(λ)
    end

    return real_part / k, imag_part / k
end

function tilde_B0_ref(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64; n_smooth::Int, n_sing::Int)
    if abs(k) < EPS_K
        return tilde_B0_k_zero_ref(sign_flag, λ, m, m_prime, μ, T, Φ, Φbar; n_smooth=n_smooth, n_sing=n_sing)
    else
        return tilde_B0_k_positive_ref(sign_flag, λ, k, m, m_prime, μ, T, Φ, Φbar; n_smooth=n_smooth, n_sing=n_sing)
    end
end

function B0_ref(λ::Float64, k::Float64, m1::Float64, μ1::Float64, m2::Float64, μ2::Float64, T::Float64;
    Φ::Float64=0.0, Φbar::Float64=0.0, n_smooth::Int=64, n_sing::Int=128)
    term1 = tilde_B0_ref(:plus, -λ, k, m1, m2, μ1, T, Φ, Φbar; n_smooth=n_smooth, n_sing=n_sing)
    term2 = tilde_B0_ref(:minus, λ, k, m1, m2, μ1, T, Φ, Φbar; n_smooth=n_smooth, n_sing=n_sing)
    term3 = tilde_B0_ref(:plus, λ, k, m2, m1, μ2, T, Φ, Φbar; n_smooth=n_smooth, n_sing=n_sing)
    term4 = tilde_B0_ref(:minus, -λ, k, m2, m1, μ2, T, Φ, Φbar; n_smooth=n_smooth, n_sing=n_sing)
    real_part = term1[1] - term2[1] + term3[1] - term4[1]
    imag_part = term1[2] - term2[2] + term3[2] - term4[2]
    return real_part, imag_part
end

@testset "B0 convergence" begin
    imag_nonzero = 0
    for c in CASES
        B0_re_ref, B0_im_ref = B0_ref(c.λ, c.k, c.m1, c.μ1, c.m2, c.μ2, c.T; Φ=c.Φ, Φbar=c.Φbar)
        B0_re_def, B0_im_def = B0(c.λ, c.k, c.m1, c.μ1, c.m2, c.μ2, c.T; Φ=c.Φ, Φbar=c.Φbar)

        @test rel_diff(B0_re_def, B0_re_ref) < TOL
        @test rel_diff(B0_im_def, B0_im_ref) < TOL

        if abs(B0_im_ref) > 1e-8
            imag_nonzero += 1
        end
    end
    @test imag_nonzero >= 1
end
