"""
    ConstraintComponents

约束组件契约与默认组件映射。
"""

abstract type AbstractConstraintComponent end

struct StationarityComponent <: AbstractConstraintComponent end
struct EqualMuComponent <: AbstractConstraintComponent end
struct FixedBaryonDensityComponent <: AbstractConstraintComponent end
struct FixedEntropyComponent <: AbstractConstraintComponent end
struct FixedSigmaComponent <: AbstractConstraintComponent end
struct AsymmetricDensityComponent <: AbstractConstraintComponent end
struct FixedMuBComponent <: AbstractConstraintComponent end
struct ConservedChargeDensityComponent <: AbstractConstraintComponent end

@inline constraint_name(::StationarityComponent) = :stationarity
@inline constraint_name(::EqualMuComponent) = :equal_mu
@inline constraint_name(::FixedBaryonDensityComponent) = :fixed_baryon_density
@inline constraint_name(::FixedEntropyComponent) = :fixed_entropy
@inline constraint_name(::FixedSigmaComponent) = :fixed_sigma
@inline constraint_name(::AsymmetricDensityComponent) = :asymmetric_density
@inline constraint_name(::FixedMuBComponent) = :fixed_muB
@inline constraint_name(::ConservedChargeDensityComponent) = :conserved_charge_density

@inline constraint_dim(::StationarityComponent) = 5
@inline constraint_dim(::EqualMuComponent) = 2
@inline constraint_dim(::FixedBaryonDensityComponent) = 1
@inline constraint_dim(::FixedEntropyComponent) = 1
@inline constraint_dim(::FixedSigmaComponent) = 1
@inline constraint_dim(::AsymmetricDensityComponent) = 3
@inline constraint_dim(::FixedMuBComponent) = 1
@inline constraint_dim(::ConservedChargeDensityComponent) = 2

@inline function build_constraint_components(::FixedMu)
    return AbstractConstraintComponent[
        StationarityComponent(),
    ]
end

@inline function build_constraint_components(::FixedRho)
    return AbstractConstraintComponent[
        StationarityComponent(),
        EqualMuComponent(),
        FixedBaryonDensityComponent(),
    ]
end

@inline function build_constraint_components(::FixedAsymmetricRho)
    return AbstractConstraintComponent[
        StationarityComponent(),
        AsymmetricDensityComponent(),
    ]
end

@inline function build_constraint_components(::FixedMuBConservedCharges)
    return AbstractConstraintComponent[
        StationarityComponent(),
        FixedMuBComponent(),
        ConservedChargeDensityComponent(),
    ]
end

@inline function build_constraint_components(::FixedEntropy)
    return AbstractConstraintComponent[
        StationarityComponent(),
        EqualMuComponent(),
        FixedEntropyComponent(),
    ]
end

@inline function build_constraint_components(::FixedSigma)
    return AbstractConstraintComponent[
        StationarityComponent(),
        EqualMuComponent(),
        FixedSigmaComponent(),
    ]
end

@inline constraint_total_dim(components::AbstractVector{<:AbstractConstraintComponent}) = sum(constraint_dim, components)
