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
state_dim(::FixedEntropy) = 8
state_dim(::FixedSigma) = 8

param_dim(::FixedMu) = 2
param_dim(::FixedRho) = 1
param_dim(::FixedAsymmetricRho) = 1
param_dim(::FixedEntropy) = 1
param_dim(::FixedSigma) = 1

constraint_description(::FixedMu) = "Fixed chemical potential μ"
constraint_description(m::FixedRho) = "Fixed baryon density ρ/ρ₀ = $(m.rho_target)"
constraint_description(m::FixedAsymmetricRho) = "Fixed asymmetric density constraints (ρ/ρ₀=$(m.rho_target), ρu/ρd=$(m.ud_ratio_target), ρs=$(m.s_target))"
constraint_description(m::FixedEntropy) = "Fixed entropy density s = $(m.s_target) fm⁻³"
constraint_description(m::FixedSigma) = "Fixed specific entropy σ = s/n_B = $(m.sigma_target)"

Base.show(io::IO, ::FixedMu) = print(io, "FixedMu()")
Base.show(io::IO, m::FixedRho) = print(io, "FixedRho(ρ/ρ₀=$(m.rho_target))")
Base.show(io::IO, m::FixedAsymmetricRho) = print(io, "FixedAsymmetricRho(ρ/ρ₀=$(m.rho_target), ρu/ρd=$(m.ud_ratio_target), ρs=$(m.s_target))")
Base.show(io::IO, m::FixedEntropy) = print(io, "FixedEntropy(s=$(m.s_target))")
Base.show(io::IO, m::FixedSigma) = print(io, "FixedSigma(σ=$(m.sigma_target))")
