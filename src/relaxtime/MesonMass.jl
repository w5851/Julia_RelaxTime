"""
# MesonMass.jl

介子质量计算模块（π、K、σ_π、σ_K；以及混合介子 η/η′、σ/σ′）。
采用实数 B0 积分 + 宽度耦合的极点条件：
    1 - 4K * Π(p0) = 0,  p0 = M + iΓ/2
"""
module MesonMass

include("../integration/GaussLegendre.jl")
include("../Constants_PNJL.jl")
include("OneLoopIntegrals.jl")
include("OneLoopIntegralsAniso.jl")
include("EffectiveCouplings.jl")
include("PolarizationAniso.jl")

using .GaussLegendre: gauleg, DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS
using .Constants_PNJL: G_fm2, K_fm5
using .OneLoopIntegrals: A
using .OneLoopIntegralsCorrection: A_aniso
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, mixing_matrix_elements
using .PolarizationAniso: polarization_with_width
using NLsolve

export ensure_quark_params_has_A
export default_meson_mass_guess
export meson_mass_equation, solve_meson_mass

@inline _get_xi(thermo_params::NamedTuple)::Float64 = hasproperty(thermo_params, :ξ) ? Float64(thermo_params.ξ) : 0.0

@inline function _is_mixed_meson(meson::Symbol)::Bool
    return (meson == :eta) || (meson == :eta_prime) || (meson == :sigma) || (meson == :sigma_prime)
end

@inline function _mixed_channel(meson::Symbol)::Symbol
    if (meson == :eta) || (meson == :eta_prime)
        return :P
    elseif (meson == :sigma) || (meson == :sigma_prime)
        return :S
    else
        error("Not a mixed meson: $meson")
    end
end

@inline function _mixed_matrix_elements(channel::Symbol, k0::Float64, gamma::Float64, k_norm::Float64,
                                        quark_params::NamedTuple, thermo_params::NamedTuple,
                                        K_coeffs::NamedTuple)
    ξ = _get_xi(thermo_params)

    Πuu_re, Πuu_im = polarization_with_width(
        channel, k0, gamma, k_norm,
        quark_params.m.u, quark_params.m.u,
        quark_params.μ.u, quark_params.μ.u,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, ξ,
        quark_params.A.u, quark_params.A.u,
        0,
    )
    Πss_re, Πss_im = polarization_with_width(
        channel, k0, gamma, k_norm,
        quark_params.m.s, quark_params.m.s,
        quark_params.μ.s, quark_params.μ.s,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, ξ,
        quark_params.A.s, quark_params.A.s,
        0,
    )

    Πuu = ComplexF64(Πuu_re, Πuu_im)
    Πss = ComplexF64(Πss_re, Πss_im)

    elems = mixing_matrix_elements(Πuu, Πss, K_coeffs, channel)
    return (M00=elems.M00, M08=elems.M08, M88=elems.M88, Πuu=Πuu, Πss=Πss)
end

@inline function _mixed_inverse(channel::Symbol, which::Symbol, elems)::ComplexF64
    M00 = elems.M00
    M08 = elems.M08
    M88 = elems.M88
    sqrt_term = sqrt((M00 - M88)^2 + 4.0 * M08^2)
    if which == :light
        return (M00 + M88) - sqrt_term
    elseif which == :heavy
        return (M00 + M88) + sqrt_term
    else
        error("Unknown mixed eigenmode: $which")
    end
end

@inline function ensure_quark_params_has_A(quark_params::NamedTuple, thermo_params::NamedTuple;
                                           p_nodes::Int=16,
                                           p_max::Float64=20.0,
                                           cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
                                           use_aniso::Bool=true)
    if hasproperty(quark_params, :A)
        return quark_params
    end
    hasproperty(quark_params, :m) || error("quark_params is missing :m")
    hasproperty(quark_params, :μ) || error("quark_params is missing :μ")
    hasproperty(thermo_params, :T) || error("thermo_params is missing :T")
    hasproperty(thermo_params, :Φ) || error("thermo_params is missing :Φ")
    hasproperty(thermo_params, :Φbar) || error("thermo_params is missing :Φbar")

    nodes_p, weights_p = gauleg(0.0, p_max, p_nodes)
    ξ = hasproperty(thermo_params, :ξ) ? thermo_params.ξ : 0.0
    if use_aniso && abs(ξ) > 0.0
        nodes_cos, weights_cos = gauleg(-1.0, 1.0, cos_nodes)
        A_u = A_aniso(quark_params.m.u, quark_params.μ.u, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
                      ξ, nodes_p, weights_p, nodes_cos, weights_cos)
        A_d = A_aniso(quark_params.m.d, quark_params.μ.d, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
                      ξ, nodes_p, weights_p, nodes_cos, weights_cos)
        A_s = A_aniso(quark_params.m.s, quark_params.μ.s, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
                      ξ, nodes_p, weights_p, nodes_cos, weights_cos)
    else
        A_u = A(quark_params.m.u, quark_params.μ.u, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, nodes_p, weights_p)
        A_d = A(quark_params.m.d, quark_params.μ.d, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, nodes_p, weights_p)
        A_s = A(quark_params.m.s, quark_params.μ.s, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, nodes_p, weights_p)
    end

    return merge(quark_params, (A=(u=A_u, d=A_d, s=A_s),))
end

@inline function _ensure_K_coeffs(quark_params::NamedTuple, K_coeffs::Union{Nothing,NamedTuple})
    if K_coeffs !== nothing
        return K_coeffs
    end
    G_u = calculate_G_from_A(quark_params.A.u, quark_params.m.u)
    G_s = calculate_G_from_A(quark_params.A.s, quark_params.m.s)
    return calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

@inline function _meson_polarization_params(meson::Symbol, quark_params::NamedTuple)
    if meson == :pi
        return (channel = :P,
                m1 = quark_params.m.u, m2 = quark_params.m.u,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.u,
                A1 = quark_params.A.u, A2 = quark_params.A.u,
                num_s_quark = 0)
    elseif meson == :K
        return (channel = :P,
                m1 = quark_params.m.u, m2 = quark_params.m.s,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.s,
                A1 = quark_params.A.u, A2 = quark_params.A.s,
                num_s_quark = 1)
    elseif meson == :sigma_pi
        return (channel = :S,
                m1 = quark_params.m.u, m2 = quark_params.m.u,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.u,
                A1 = quark_params.A.u, A2 = quark_params.A.u,
                num_s_quark = 0)
    elseif meson == :sigma_K
        return (channel = :S,
                m1 = quark_params.m.u, m2 = quark_params.m.s,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.s,
                A1 = quark_params.A.u, A2 = quark_params.A.s,
                num_s_quark = 1)
    else
        error("Unknown meson type: $meson. Use :pi, :K, :sigma_pi, or :sigma_K")
    end
end

@inline function _meson_coupling(meson::Symbol, K_coeffs::NamedTuple)::Float64
    if meson == :pi
        return K_coeffs.K123_plus
    elseif meson == :K
        return K_coeffs.K4567_plus
    elseif meson == :sigma_pi
        return K_coeffs.K123_minus
    elseif meson == :sigma_K
        return K_coeffs.K4567_minus
    else
        error("Unknown meson type: $meson. Use :pi, :K, :sigma_pi, or :sigma_K")
    end
end

@inline function default_meson_mass_guess(meson::Symbol, quark_params::NamedTuple)::Float64
    if meson == :pi
        return quark_params.m.u + quark_params.m.d
    elseif meson == :K
        return quark_params.m.u + quark_params.m.s
    elseif meson == :sigma_pi
        return 2.0 * quark_params.m.u
    elseif meson == :sigma_K
        return quark_params.m.u + quark_params.m.s
    elseif meson == :eta
        return 2.0 * quark_params.m.u
    elseif meson == :eta_prime
        return 2.0 * quark_params.m.s
    elseif meson == :sigma
        return 2.0 * quark_params.m.u
    elseif meson == :sigma_prime
        return 2.0 * quark_params.m.s
    else
        error("Unknown meson type: $meson. Use :pi, :K, :sigma_pi, :sigma_K, :eta, :eta_prime, :sigma, or :sigma_prime")
    end
end

"""
    meson_mass_equation(meson, k0, gamma, k_norm, quark_params, thermo_params, K_coeffs) -> ComplexF64

返回介子极点方程的复数残差：f = 1 - 4K * Π。
"""
function meson_mass_equation(meson::Symbol, k0::Float64, gamma::Float64, k_norm::Float64,
                             quark_params::NamedTuple, thermo_params::NamedTuple,
                             K_coeffs::NamedTuple)
    ξ = _get_xi(thermo_params)

    if _is_mixed_meson(meson)
        channel = _mixed_channel(meson)
        elems = _mixed_matrix_elements(channel, k0, gamma, k_norm, quark_params, thermo_params, K_coeffs)
        which = ((meson == :eta) || (meson == :sigma)) ? :light : :heavy
        return _mixed_inverse(channel, which, elems)
    end

    pol_params = _meson_polarization_params(meson, quark_params)
    K = _meson_coupling(meson, K_coeffs)

    Π_real, Π_imag = polarization_with_width(
        pol_params.channel, k0, gamma, k_norm,
        pol_params.m1, pol_params.m2,
        pol_params.μ1, pol_params.μ2,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, ξ,
        pol_params.A1, pol_params.A2,
        pol_params.num_s_quark,
    )
    Π = ComplexF64(Π_real, Π_imag)
    return 1.0 - 4.0 * K * Π
end

@inline function _mass_residual!(fvec::Vector{Float64}, x::Vector{Float64},
                                 meson::Symbol, k_norm::Float64,
                                 quark_params::NamedTuple, thermo_params::NamedTuple,
                                 K_coeffs::NamedTuple)
    k0 = x[1]
    gamma = x[2]
    f = meson_mass_equation(meson, k0, gamma, k_norm, quark_params, thermo_params, K_coeffs)
    fvec[1] = real(f)
    fvec[2] = imag(f)
    return nothing
end

"""
    solve_meson_mass(meson, quark_params, thermo_params; K_coeffs=nothing,
                     k_norm=0.0, initial_mass=nothing, initial_gamma=0.0,
                     nlsolve_kwargs...) -> NamedTuple

求解介子极点方程，返回 (mass, gamma, converged, residual_norm, solution)。
"""
function solve_meson_mass(meson::Symbol, quark_params::NamedTuple, thermo_params::NamedTuple;
                          K_coeffs::Union{Nothing,NamedTuple}=nothing,
                          k_norm::Float64=0.0,
                          initial_mass::Union{Nothing,Float64}=nothing,
                          initial_gamma::Float64=0.0,
                          nlsolve_kwargs...)
    quark_params = ensure_quark_params_has_A(quark_params, thermo_params)
    K_coeffs = _ensure_K_coeffs(quark_params, K_coeffs)

    m0 = initial_mass === nothing ? default_meson_mass_guess(meson, quark_params) : initial_mass
    x0 = [m0, initial_gamma]

    result = nlsolve((f, x) -> _mass_residual!(f, x, meson, k_norm, quark_params, thermo_params, K_coeffs), x0;
                     nlsolve_kwargs...)

    mass = result.zero[1]
    gamma = result.zero[2]
    residual_norm = result.residual_norm

    converged = if hasproperty(result, :converged)
        result.converged
    else
        (hasproperty(result, :f_converged) ? result.f_converged : false) ||
            (hasproperty(result, :x_converged) ? result.x_converged : false)
    end

    return (mass=mass, gamma=gamma, converged=converged,
            residual_norm=residual_norm, solution=result)
end

end # module MesonMass