"""
    MesonRPA

Pure algebra for the neutral three-flavor RPA backend.  The module consumes a
`FullKMTInteraction` produced on a fixed mean-field background and three
flavor-resolved quark-bubble polarizations.  It does not solve the PNJL gap
equations, evaluate the bubbles, or alter the legacy propagator interfaces.

The matrix convention follows Tian et al., Phys. Rev. D 114, 034012 (2026),
Eqs. (20)-(22):

```text
M = 2 K * inv(I - 2 K * Pi)
```

The order of multiplication is intentional because the full matrices need
not commute away from the isospin-symmetric limit.
"""
module MesonRPA

using LinearAlgebra: det, I
using StaticArrays
using ..MesonInteractionKernel: FullKMTInteraction, neutral_coupling_matrix

export neutral_polarization_matrix
export neutral_rpa_inverse_matrix, neutral_rpa_propagator, neutral_rpa_determinant

const _NEUTRAL_SIZE = 3

@inline _real_type(x::Real) = typeof(x)
@inline _real_type(x::Complex) = typeof(real(x))

@inline function _complex_finite(x::Number)
    return isfinite(real(x)) && isfinite(imag(x))
end

function _coerce_flavor_bubbles(bubbles::SVector{3,<:Number})
    return bubbles
end

function _coerce_flavor_bubbles(bubbles::NTuple{3,<:Number})
    return SVector{3}(bubbles)
end

function _coerce_flavor_bubbles(bubbles::AbstractVector{<:Number})
    length(bubbles) == 3 || throw(ArgumentError("bubbles must contain exactly three entries (u, d, s)"))
    return SVector{3}(bubbles)
end

function _coerce_flavor_bubbles(bubbles::NamedTuple)
    all(name -> hasproperty(bubbles, name), (:u, :d, :s)) ||
        throw(ArgumentError("bubbles NamedTuple must provide fields :u, :d, and :s"))
    return _coerce_flavor_bubbles((bubbles.u, bubbles.d, bubbles.s))
end

_coerce_flavor_bubbles(bubbles) = throw(ArgumentError(
    "bubbles must be a length-three real/complex vector/tuple or a NamedTuple with fields :u, :d, :s",
))

function _promote_bubble_type(bubbles::SVector{3,<:Number})
    T = promote_type((_real_type(value) for value in bubbles)...)
    T <: Real || throw(ArgumentError("bubble values must have real-valued real/imaginary parts"))
    return T
end

@inline function _validate_bubbles(bubbles::SVector{3,<:Number})
    all(_complex_finite, bubbles) || throw(ArgumentError("bubble values must be finite"))
    return nothing
end

"""
    neutral_polarization_matrix(bubbles)

Build the pseudoscalar/scalar neutral polarization matrix in the `(0, 3, 8)`
basis from flavor-diagonal bubbles `(Pi_u, Pi_d, Pi_s)`.  The coefficients
follow Tian et al. (2026), Eq. (26), and are valid for real or complex bubble
values.
"""
function neutral_polarization_matrix(bubbles)
    values = _coerce_flavor_bubbles(bubbles)
    _validate_bubbles(values)
    T = _promote_bubble_type(values)
    u, d, s = (Complex{T}(value) for value in values)

    pi0 = (T(2) / T(3)) * (u + d + s)
    pi3 = u + d
    pi8 = (one(T) / T(3)) * (u + d + T(4) * s)
    pi03 = (sqrt(T(6)) / T(3)) * (u - d)
    pi08 = (sqrt(T(2)) / T(3)) * (u + d - T(2) * s)
    pi38 = (sqrt(T(3)) / T(3)) * (u - d)

    # SMatrix is column-major: columns are (Pi0,Pi03,Pi08),
    # (Pi03,Pi3,Pi38), and (Pi08,Pi38,Pi8).
    return SMatrix{3,3,Complex{T},9}((
        pi0, pi03, pi08,
        pi03, pi3, pi38,
        pi08, pi38, pi8,
    ))
end

function _coerce_polarization_matrix(polarization::SMatrix{3,3,<:Number,9})
    return polarization
end

function _coerce_polarization_matrix(polarization::AbstractMatrix{<:Number})
    size(polarization) == (_NEUTRAL_SIZE, _NEUTRAL_SIZE) ||
        throw(ArgumentError("polarization must be a 3x3 matrix in the (0, 3, 8) basis"))
    return SMatrix{3,3}(polarization)
end

_coerce_polarization_matrix(polarization) = throw(ArgumentError(
    "polarization must be a 3x3 real/complex matrix in the (0, 3, 8) basis",
))

@inline function _validate_polarization(polarization::SMatrix{3,3,<:Number,9})
    all(_complex_finite, polarization) || throw(ArgumentError("polarization entries must be finite"))
    return nothing
end

@inline function _coerce_channel_matrix(kernel::FullKMTInteraction, channel::Symbol)
    K = neutral_coupling_matrix(kernel, channel)
    T = promote_type(eltype(K), Float64)
    return SMatrix{3,3,T,9}(K)
end

function _promoted_matrix_type(K::SMatrix{3,3,<:Real,9}, Pi::SMatrix{3,3,<:Number,9})
    T = promote_type(eltype(K), (_real_type(value) for value in Pi)...)
    T <: Real || throw(ArgumentError("matrix values must have real-valued real/imaginary parts"))
    return T
end

"""
    neutral_rpa_inverse_matrix(kernel, polarization; channel=:P)

Return the RPA denominator matrix `I - 2 K Pi` for the neutral `(0, 3, 8)`
sector.  The matrix order is the one used by the source equation and is kept
explicit for later pole/phase-shift consumers.
"""
function neutral_rpa_inverse_matrix(
    kernel::FullKMTInteraction,
    polarization;
    channel::Symbol=:P,
)
    K_real = _coerce_channel_matrix(kernel, channel)
    Pi_raw = _coerce_polarization_matrix(polarization)
    _validate_polarization(Pi_raw)
    T = _promoted_matrix_type(K_real, Pi_raw)
    K = SMatrix{3,3,Complex{T},9}(K_real)
    Pi = SMatrix{3,3,Complex{T},9}(Pi_raw)
    identity = SMatrix{3,3,Complex{T},9}(I)
    result = identity - T(2) * K * Pi
    all(_complex_finite, result) || throw(ArgumentError("RPA denominator matrix is not finite"))
    return result
end

"""
    neutral_rpa_propagator(kernel, polarization; channel=:P)

Compute the full neutral RPA matrix
`2 K * inv(I - 2 K * Pi)` in the `(0, 3, 8)` basis.  This is a diagnostic
backend: it does not locate poles, unwrap phase shifts, or feed meson pressure
back into the mean-field equations.
"""
function neutral_rpa_propagator(
    kernel::FullKMTInteraction,
    polarization;
    channel::Symbol=:P,
)
    K_real = _coerce_channel_matrix(kernel, channel)
    Pi_raw = _coerce_polarization_matrix(polarization)
    _validate_polarization(Pi_raw)
    T = _promoted_matrix_type(K_real, Pi_raw)
    K = SMatrix{3,3,Complex{T},9}(K_real)
    denominator = neutral_rpa_inverse_matrix(kernel, Pi_raw; channel=channel)
    result = T(2) * K * inv(denominator)
    all(_complex_finite, result) || throw(ArgumentError("RPA propagator matrix is not finite"))
    return result
end

"""Return `det(I - 2 K Pi)` for the neutral RPA denominator."""
function neutral_rpa_determinant(
    kernel::FullKMTInteraction,
    polarization;
    channel::Symbol=:P,
)
    denominator = neutral_rpa_inverse_matrix(kernel, polarization; channel=channel)
    result = det(denominator)
    _complex_finite(result) || throw(ArgumentError("RPA denominator determinant is not finite"))
    return result
end

end # module MesonRPA
