"""
# MottTransition.jl

Mott 相变判据与差值计算。
"""
module MottTransition

export mott_threshold_mass, mott_gap, is_mott_point
export mott_threshold_masses, mott_gaps

@inline function _is_mixed_meson(meson::Symbol)::Bool
    return (meson == :eta) || (meson == :eta_prime) || (meson == :sigma) || (meson == :sigma_prime)
end

"""Return all relevant quark-pair thresholds for a meson.

For non-mixed mesons, returns a NamedTuple with a single field `threshold`.
For mixed mesons (η/η′ and σ/σ′), returns thresholds for uu and ss channels.
"""
@inline function mott_threshold_masses(meson::Symbol, quark_params::NamedTuple)
    if _is_mixed_meson(meson)
        thr_uu = 2.0 * quark_params.m.u
        thr_ss = 2.0 * quark_params.m.s
        return (uu=thr_uu, ss=thr_ss, min=min(thr_uu, thr_ss))
    end

    return (threshold=mott_threshold_mass(meson, quark_params),)
end

@inline function mott_threshold_mass(meson::Symbol, quark_params::NamedTuple)::Float64
    if meson == :pi
        return quark_params.m.u + quark_params.m.d
    elseif meson == :K
        return quark_params.m.u + quark_params.m.s
    elseif meson == :sigma_pi
        return quark_params.m.u + quark_params.m.d
    elseif meson == :sigma_K
        return quark_params.m.u + quark_params.m.s
    elseif _is_mixed_meson(meson)
        # Mixed mesons can decay into uu or ss; keep the legacy API by using the lower threshold.
        thr = mott_threshold_masses(meson, quark_params)
        return thr.min
    else
        error("Unknown meson type: $meson. Use :pi, :K, :sigma_pi, :sigma_K, :eta, :eta_prime, :sigma, or :sigma_prime")
    end
end

"""
    mott_gap(meson, meson_mass, quark_params) -> Float64

返回介子质量与组分夸克阈值之差：Δ = M_mes - (M_{q1}+M_{q2})。
"""
@inline function mott_gap(meson::Symbol, meson_mass::Float64, quark_params::NamedTuple)::Float64
    return meson_mass - mott_threshold_mass(meson, quark_params)
end

"""Return Mott gaps for all relevant thresholds.

For non-mixed mesons, returns `(gap=...)`.
For mixed mesons, returns `(uu=..., ss=..., min=...)`.
"""
@inline function mott_gaps(meson::Symbol, meson_mass::Float64, quark_params::NamedTuple)
    thr = mott_threshold_masses(meson, quark_params)
    if hasproperty(thr, :threshold)
        return (gap=meson_mass - thr.threshold,)
    end
    gap_uu = meson_mass - thr.uu
    gap_ss = meson_mass - thr.ss
    return (uu=gap_uu, ss=gap_ss, min=meson_mass - thr.min)
end

"""
    is_mott_point(meson, meson_mass, quark_params; atol=1e-6) -> Bool

判断是否满足 Mott 相变条件。
"""
@inline function is_mott_point(meson::Symbol, meson_mass::Float64, quark_params::NamedTuple;
                               atol::Float64=1e-6)::Bool
    return abs(mott_gap(meson, meson_mass, quark_params)) <= atol
end

end # module MottTransition