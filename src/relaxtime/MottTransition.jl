"""
# MottTransition.jl

Mott 相变判据与差值计算。
"""
module MottTransition

export mott_threshold_mass, mott_gap, is_mott_point

@inline function mott_threshold_mass(meson::Symbol, quark_params::NamedTuple)::Float64
    if meson == :pi
        return quark_params.m.u + quark_params.m.d
    elseif meson == :K
        return quark_params.m.u + quark_params.m.s
    elseif meson == :sigma_pi
        return quark_params.m.u + quark_params.m.d
    elseif meson == :sigma_K
        return quark_params.m.u + quark_params.m.s
    else
        error("Unknown meson type: $meson. Use :pi, :K, :sigma_pi, or :sigma_K")
    end
end

"""
    mott_gap(meson, meson_mass, quark_params) -> Float64

返回介子质量与组分夸克阈值之差：Δ = M_mes - (M_{q1}+M_{q2})。
"""
@inline function mott_gap(meson::Symbol, meson_mass::Float64, quark_params::NamedTuple)::Float64
    return meson_mass - mott_threshold_mass(meson, quark_params)
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