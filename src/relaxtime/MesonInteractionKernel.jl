"""
    MesonInteractionKernel

Phase-1/1.5 pure algebra backend for the three-flavor KMT meson interaction
kernel.  The backend is deliberately parallel to `EffectiveCouplings`: it
does not alter the PNJL mean-field solver or the legacy 0/8 propagator API.

The input condensates use the convention
`phi_f = <bar(q_f) q_f>`.  The P/S signs and the neutral `(0, 3, 8)` basis are
recorded in `KernelConvention` so that a later RPA implementation can audit
normalization and sign conversions before using the matrices in production.
"""
module MesonInteractionKernel

using StaticArrays
using Main.Constants_PNJL: G_fm2, K_fm5

export KernelConvention, ChargedKMTChannels, FullKMTInteraction
export default_kernel_convention
export build_full_kmt_interaction, build_full_kmt_interaction_from_equilibrium
export neutral_coupling_matrix, charged_couplings, charged_coupling, get_Kab

const _NEUTRAL_BASIS = (:lambda0, :lambda3, :lambda8)

"""Metadata that fixes the convention used to generate a KMT interaction kernel."""
struct KernelConvention
    name::Symbol
    lambda0_normalization::Symbol
    condensate_definition::Symbol
    channel_sign_convention::Symbol
    source::Symbol
end

"""Return the source-aligned convention used by the full-KMT builder."""
@inline function default_kernel_convention()
    return KernelConvention(
        :full_kmt_tian2026_eq3,
        :sqrt_2_over_3_identity,
        :phi_equals_sigma,
        :pseudoscalar_plus_scalar_minus,
        :tian_et_al_2026_eq3,
    )
end

"""Flavor-resolved charged/kaon channel couplings in one P/S sector."""
struct ChargedKMTChannels{T<:Real}
    K12::T
    K45::T
    K67::T
end

"""Full candidate KMT interaction kernel on a fixed mean-field background.

`neutral_P` and `neutral_S` use the row/column order `(0, 3, 8)`.  `charged_P`
and `charged_S` retain the three flavor channels `K12`, `K45`, and `K67`.
"""
struct FullKMTInteraction{T<:Real}
    G::T
    K::T
    phi::SVector{3,T}
    neutral_P::SMatrix{3,3,T,9}
    neutral_S::SMatrix{3,3,T,9}
    charged_P::ChargedKMTChannels{T}
    charged_S::ChargedKMTChannels{T}
    neutral_basis::NTuple{3,Symbol}
    convention::KernelConvention
end

@inline function _coerce_phi(phi::SVector{3,T}) where {T<:Real}
    return phi
end

function _coerce_phi(phi::AbstractVector{<:Real})
    length(phi) == 3 || throw(ArgumentError("phi must contain exactly three entries (u, d, s)"))
    T = promote_type(typeof(phi[1]), typeof(phi[2]), typeof(phi[3]))
    return SVector{3,T}(phi[1], phi[2], phi[3])
end

function _coerce_phi(phi::NTuple{3,<:Real})
    T = promote_type(typeof(phi[1]), typeof(phi[2]), typeof(phi[3]))
    return SVector{3,T}(phi[1], phi[2], phi[3])
end

function _coerce_phi(phi::NamedTuple)
    all(name -> hasproperty(phi, name), (:u, :d, :s)) ||
        throw(ArgumentError("phi NamedTuple must provide fields :u, :d, and :s"))
    return _coerce_phi((phi.u, phi.d, phi.s))
end

_coerce_phi(phi) = throw(ArgumentError("phi must be a length-three real vector/tuple or a NamedTuple with fields :u, :d, :s"))

@inline function _validate_inputs(phi::SVector{3,T}, G::T, K::T) where {T<:Real}
    all(isfinite, phi) || throw(ArgumentError("phi entries must be finite"))
    isfinite(G) || throw(ArgumentError("G must be finite"))
    isfinite(K) || throw(ArgumentError("K must be finite"))
    return nothing
end

@inline function _channel_sign(channel::Symbol)
    if channel === :P || channel === :pseudoscalar || channel === :plus
        return 1.0
    elseif channel === :S || channel === :scalar || channel === :minus
        return -1.0
    end
    throw(ArgumentError("unknown channel $(channel); use :P/:S (or :pseudoscalar/:scalar)"))
end

@inline function _symmetric_neutral_matrix(
    ::Type{T},
    K0::T,
    K3::T,
    K8::T,
    K03::T,
    K08::T,
    K38::T,
) where {T<:Real}
    # SMatrix is column-major: columns are (K0,K03,K08), (K03,K3,K38),
    # and (K08,K38,K8), respectively.
    return SMatrix{3,3,T,9}((K0, K03, K08, K03, K3, K38, K08, K38, K8))
end

function _build_channel(
    phi::SVector{3,T},
    G::T,
    K::T,
    channel::Symbol,
) where {T<:Real}
    u, d, strange = phi
    sign = T(_channel_sign(channel))

    # Source-aligned convention from Tian et al., Phys. Rev. D 114, 034012
    # (2026), Eq. (3), after identifying phi_f with sigma_f.
    K0 = G + sign * (K / T(3)) * (u + d + strange)
    K3 = G - sign * (K / T(2)) * strange
    K45 = G - sign * (K / T(2)) * d
    K67 = G - sign * (K / T(2)) * u
    K8 = G - sign * (K / T(6)) * (T(2) * u + T(2) * d - strange)
    K03 = -sign * (K / (T(2) * sqrt(T(6)))) * (u - d)
    K08 = -sign * (sqrt(T(2)) * K / T(12)) * (u + d - T(2) * strange)
    K38 = sign * (K / (T(2) * sqrt(T(3)))) * (u - d)

    neutral = _symmetric_neutral_matrix(T, K0, K3, K8, K03, K08, K38)
    charged = ChargedKMTChannels{T}(K3, K45, K67)
    return neutral, charged
end

"""Build the full candidate KMT kernel from three scalar condensates.

The returned object is immutable and allocation-free after conversion to
`StaticArrays`.  `phi` may be an `SVector`, a length-three real vector/tuple,
or a `NamedTuple` with fields `u`, `d`, and `s`.
"""
function build_full_kmt_interaction(
    phi;
    G::Real=G_fm2,
    K::Real=K_fm5,
    convention::KernelConvention=default_kernel_convention(),
)
    phi_vec = _coerce_phi(phi)
    promoted_type = promote_type(eltype(phi_vec), typeof(G), typeof(K))
    # The neutral off-diagonal coefficients contain square roots.  Promote
    # an all-integer input to Float64 instead of attempting to store those
    # coefficients in an integer SMatrix; retain non-integer numeric types.
    T = (promoted_type <: Integer || promoted_type <: Rational) ? Float64 : promoted_type
    phi_T = SVector{3,T}(phi_vec)
    G_T = T(G)
    K_T = T(K)
    _validate_inputs(phi_T, G_T, K_T)

    neutral_P, charged_P = _build_channel(phi_T, G_T, K_T, :P)
    neutral_S, charged_S = _build_channel(phi_T, G_T, K_T, :S)
    return FullKMTInteraction{T}(
        G_T,
        K_T,
        phi_T,
        neutral_P,
        neutral_S,
        charged_P,
        charged_S,
        _NEUTRAL_BASIS,
        convention,
    )
end

"""Build a full kernel by reading the first three scalar fields of `x_state`.

This adapter is intentionally solver-free: it accepts the `equilibrium`
NamedTuple returned by the current equilibrium facade, extracts
`equilibrium.x_state[1:3]`, and delegates to `build_full_kmt_interaction`.
"""
function build_full_kmt_interaction_from_equilibrium(
    equilibrium;
    G::Real=G_fm2,
    K::Real=K_fm5,
    convention::KernelConvention=default_kernel_convention(),
)
    hasproperty(equilibrium, :x_state) ||
        throw(ArgumentError("equilibrium must provide an x_state field"))
    x_state = equilibrium.x_state
    length(x_state) >= 3 || throw(ArgumentError("equilibrium.x_state must contain at least three scalar condensates"))
    return build_full_kmt_interaction((x_state[1], x_state[2], x_state[3]); G=G, K=K, convention=convention)
end

@inline function _channel_matrix(kernel::FullKMTInteraction, channel::Symbol)
    if channel === :P || channel === :pseudoscalar || channel === :plus
        return kernel.neutral_P
    elseif channel === :S || channel === :scalar || channel === :minus
        return kernel.neutral_S
    end
    throw(ArgumentError("unknown channel $(channel); use :P/:S (or :pseudoscalar/:scalar)"))
end

@inline function _channel_charged(kernel::FullKMTInteraction, channel::Symbol)
    if channel === :P || channel === :pseudoscalar || channel === :plus
        return kernel.charged_P
    elseif channel === :S || channel === :scalar || channel === :minus
        return kernel.charged_S
    end
    throw(ArgumentError("unknown channel $(channel); use :P/:S (or :pseudoscalar/:scalar)"))
end

"""Return the neutral `(0,3,8)` P/S coupling matrix."""
@inline neutral_coupling_matrix(kernel::FullKMTInteraction, channel::Symbol=:P) =
    _channel_matrix(kernel, channel)

"""Return the three flavor-resolved charged/kaon P/S couplings."""
@inline charged_couplings(kernel::FullKMTInteraction, channel::Symbol=:P) =
    _channel_charged(kernel, channel)

"""Read one of `K12`, `K45`, or `K67` from a P/S channel."""
@inline function charged_coupling(
    kernel::FullKMTInteraction,
    pair::Symbol,
    channel::Symbol=:P,
)
    pair === :K12 || pair === :K45 || pair === :K67 ||
        throw(ArgumentError("unknown charged pair $(pair); use :K12, :K45, or :K67"))
    return getproperty(_channel_charged(kernel, channel), pair)
end

@inline function _neutral_index(label::Integer)
    label == 0 && return 1
    label == 3 && return 2
    label == 8 && return 3
    throw(ArgumentError("neutral K_ab labels must be in (0, 3, 8), got $(label)"))
end

"""Read a neutral `K_ab` matrix element for `a,b ∈ (0,3,8)`."""
@inline function get_Kab(
    kernel::FullKMTInteraction,
    a::Integer,
    b::Integer,
    channel::Symbol=:P,
)
    matrix = _channel_matrix(kernel, channel)
    return matrix[_neutral_index(a), _neutral_index(b)]
end

end # module MesonInteractionKernel
