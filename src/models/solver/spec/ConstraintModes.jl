"""
    ConstraintMode

约束模式抽象基类型。
"""
abstract type ConstraintMode end

"""
    FixedMu <: ConstraintMode

固定化学势模式。
"""
struct FixedMu <: ConstraintMode end

"""
    FixedRho <: ConstraintMode

固定重子密度模式。
"""
struct FixedRho <: ConstraintMode
    rho_target::Float64
end

"""
    FixedAsymmetricRho <: ConstraintMode

固定非对称约束模式（nB, ρu/ρd, ρs）。
"""
struct FixedAsymmetricRho <: ConstraintMode
    rho_target::Float64
    ud_ratio_target::Float64
    s_target::Float64
end

"""
    FixedMuBConservedCharges <: ConstraintMode

固定重子化学势 `muB_fm`，并联立求解满足给定电荷/重子比和净奇异数密度的
三味化学势。`strangeness_density_target` 使用 fm^-3。
"""
struct FixedMuBConservedCharges <: ConstraintMode
    muB_fm::Float64
    charge_to_baryon_ratio::Float64
    strangeness_density_target::Float64

    function FixedMuBConservedCharges(
        muB_fm::Real,
        charge_to_baryon_ratio::Real=0.4,
        strangeness_density_target::Real=0.0,
    )
        values = Float64(muB_fm), Float64(charge_to_baryon_ratio), Float64(strangeness_density_target)
        all(isfinite, values) || throw(ArgumentError(
            "FixedMuBConservedCharges inputs must be finite, got $(values)",
        ))
        return new(values...)
    end
end

FixedMuBConservedCharges(
    muB_fm::Real;
    charge_to_baryon_ratio::Real=0.4,
    strangeness_density_target::Real=0.0,
) = FixedMuBConservedCharges(muB_fm, charge_to_baryon_ratio, strangeness_density_target)

"""Map conserved-charge chemical potentials `(mu_B, mu_Q, mu_S)` to flavor values."""
@inline function flavor_mu_from_bqs(mu_B::Real, mu_Q::Real, mu_S::Real)
    return (
        mu_u=mu_B / 3 + 2 * mu_Q / 3,
        mu_d=mu_B / 3 - mu_Q / 3,
        mu_s=mu_B / 3 - mu_Q / 3 - mu_S,
    )
end

"""Invert flavor chemical potentials to `(mu_B, mu_Q, mu_S)`."""
@inline function conserved_mu_from_flavor(mu_u::Real, mu_d::Real, mu_s::Real)
    return (
        mu_B=mu_u + 2 * mu_d,
        mu_Q=mu_u - mu_d,
        mu_S=mu_d - mu_s,
    )
end

"""Map flavor net densities `(rho_u,rho_d,rho_s)` to conserved B/Q/S densities."""
@inline function conserved_densities_from_flavor(rho_vec::AbstractVector)
    length(rho_vec) == 3 || throw(ArgumentError(
        "conserved_densities_from_flavor expects length 3, got $(length(rho_vec))",
    ))
    rho_u, rho_d, rho_s = rho_vec
    return (
        rho_B=(rho_u + rho_d + rho_s) / 3,
        rho_Q=(2 * rho_u - rho_d - rho_s) / 3,
        rho_S=-rho_s,
    )
end

"""
    FixedEntropy <: ConstraintMode

固定熵密度模式。
"""
struct FixedEntropy <: ConstraintMode
    s_target::Float64
end

"""
    FixedSigma <: ConstraintMode

固定比熵模式（σ = s / n_B）。
"""
struct FixedSigma <: ConstraintMode
    sigma_target::Float64
end

state_dim(::FixedMu) = 5
state_dim(::FixedRho) = 8
state_dim(::FixedAsymmetricRho) = 8
state_dim(::FixedMuBConservedCharges) = 8
state_dim(::FixedEntropy) = 8
state_dim(::FixedSigma) = 8

state_var_dim(::FixedMu) = 5
state_var_dim(::FixedRho) = 5
state_var_dim(::FixedAsymmetricRho) = 5
state_var_dim(::FixedMuBConservedCharges) = 5
state_var_dim(::FixedEntropy) = 5
state_var_dim(::FixedSigma) = 5

mu_var_dim(::FixedMu) = 0
mu_var_dim(::FixedRho) = 3
mu_var_dim(::FixedAsymmetricRho) = 3
mu_var_dim(::FixedMuBConservedCharges) = 3
mu_var_dim(::FixedEntropy) = 3
mu_var_dim(::FixedSigma) = 3

solution_dim(mode::ConstraintMode) = state_var_dim(mode) + mu_var_dim(mode)

param_dim(::FixedMu) = 2
param_dim(::FixedRho) = 1
param_dim(::FixedAsymmetricRho) = 1
param_dim(::FixedMuBConservedCharges) = 1
param_dim(::FixedEntropy) = 1
param_dim(::FixedSigma) = 1

constraint_description(::FixedMu) = "Fixed chemical potential μ"
constraint_description(m::FixedRho) = "Fixed baryon density ρ/ρ₀ = $(m.rho_target)"
constraint_description(m::FixedAsymmetricRho) = "Fixed asymmetric density constraints (ρ/ρ₀=$(m.rho_target), ρu/ρd=$(m.ud_ratio_target), ρs=$(m.s_target))"
constraint_description(m::FixedMuBConservedCharges) = "Fixed mu_B=$(m.muB_fm) fm^-1 with rho_Q/rho_B=$(m.charge_to_baryon_ratio) and rho_S=$(m.strangeness_density_target) fm^-3"
constraint_description(m::FixedEntropy) = "Fixed entropy density s = $(m.s_target) fm⁻³"
constraint_description(m::FixedSigma) = "Fixed specific entropy σ = s/n_B = $(m.sigma_target)"

Base.show(io::IO, ::FixedMu) = print(io, "FixedMu()")
Base.show(io::IO, m::FixedRho) = print(io, "FixedRho(ρ/ρ₀=$(m.rho_target))")
Base.show(io::IO, m::FixedAsymmetricRho) = print(io, "FixedAsymmetricRho(ρ/ρ₀=$(m.rho_target), ρu/ρd=$(m.ud_ratio_target), ρs=$(m.s_target))")
Base.show(io::IO, m::FixedMuBConservedCharges) = print(io, "FixedMuBConservedCharges(muB_fm=$(m.muB_fm), rhoQ/rhoB=$(m.charge_to_baryon_ratio), rhoS=$(m.strangeness_density_target))")
Base.show(io::IO, m::FixedEntropy) = print(io, "FixedEntropy(s=$(m.s_target))")
Base.show(io::IO, m::FixedSigma) = print(io, "FixedSigma(σ=$(m.sigma_target))")
