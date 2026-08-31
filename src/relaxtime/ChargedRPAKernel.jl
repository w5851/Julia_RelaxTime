"""
    ChargedRPAKernel

Contract and pure algebra for the charged single-channel RPA backend.

The module records the flavor ordering and KMT channel used by each charged
`pi^±`/`K^±` mode, then evaluates a scalar RPA numerator/denominator without
computing the quark bubble.  It is intentionally parallel to
`MesonInteractionKernel` and `MesonRPA`: the caller supplies a fixed
mean-field background, a charged coupling, and a polarization value.

The normalization is fixed by the charged-ladder trace and Goldstone gate in
the formula-route closure package: `D = 2K/(1 - 4K Pi)`.  This is a diagnostic
Phase-B backend.  It does not provide a retarded bubble,
locate poles, construct phase shifts, perform Beth-Uhlenbeck integration, or
modify the legacy `MesonPropagator`/`MesonDensity` entrypoints.
"""
module ChargedRPAKernel

import ..MesonInteractionKernel
using ..MesonInteractionKernel: FullKMTInteraction

export charged_rpa_spec, charged_rpa_coupling
export charged_rpa_inverse, charged_rpa_propagator

const _CHARGED_MESONS = (:pi_plus, :pi_minus, :K_plus, :K_minus)
const _VALID_KERNEL_PAIRS = (:K12, :K45)
const _NUMERATOR_FACTOR = 2.0
const _DENOMINATOR_FACTOR = 4.0
const _NORMALIZATION_SOURCE = :charged_ladder_goldstone
const _RETARDED_CONVENTION = :retarded_e_minus_iwt

"""Return the canonical P/S label used by the charged-RPA contract."""
@inline function _canonical_channel(channel::Symbol)
    if channel === :P || channel === :pseudoscalar || channel === :plus
        return :P
    elseif channel === :S || channel === :scalar || channel === :minus
        return :S
    end
    throw(ArgumentError("unknown channel $(channel); use :P/:S (or :pseudoscalar/:scalar)"))
end

"""Metadata for a charged meson flavor ordering and KMT pair."""
@inline function _charged_metadata(meson::Symbol)
    if meson === :pi_plus
        return ((:u, :d), :K12)
    elseif meson === :pi_minus
        return ((:d, :u), :K12)
    elseif meson === :K_plus
        return ((:u, :s), :K45)
    elseif meson === :K_minus
        return ((:s, :u), :K45)
    end
    throw(ArgumentError(
        "unknown charged meson $(meson); use :pi_plus, :pi_minus, :K_plus, or :K_minus",
    ))
end

@inline function _finite_number(x::Number)
    return isfinite(real(x)) && isfinite(imag(x))
end

"""
    ChargedRPAKernelSpec

Immutable metadata for one charged scalar/pseudoscalar RPA mode.

Fields:

- `meson`: one of `:pi_plus`, `:pi_minus`, `:K_plus`, `:K_minus`;
- `pair`: ordered quark flavors supplied to the polarization provider;
- `channel`: canonical `:P` or `:S`;
- `kernel_pair`: charged KMT pair (`:K12` or `:K45`);
- `retarded_convention`: fixed metadata only; no analytic continuation is performed;
- `normalization_source`: provenance label for the chosen scalar convention.

Use [`charged_rpa_spec`](@ref) for the validated standard mappings.
"""
struct ChargedRPAKernelSpec
    meson::Symbol
    pair::NTuple{2,Symbol}
    channel::Symbol
    kernel_pair::Symbol
    retarded_convention::Symbol
    normalization_source::Symbol

    function ChargedRPAKernelSpec(
        meson::Symbol,
        pair::NTuple{2,Symbol},
        channel::Symbol,
        kernel_pair::Symbol,
    )
        return new(
            meson,
            pair,
            channel,
            kernel_pair,
            _RETARDED_CONVENTION,
            _NORMALIZATION_SOURCE,
        )
    end
end

function _validated_spec(
    meson::Symbol,
    pair::NTuple{2,Symbol},
    channel::Symbol,
    kernel_pair::Symbol,
)
    meson in _CHARGED_MESONS || throw(ArgumentError("meson is not a supported charged mode"))
    pair[1] in (:u, :d, :s) && pair[2] in (:u, :d, :s) ||
        throw(ArgumentError("pair entries must be one of :u, :d, or :s"))
    pair[1] != pair[2] || throw(ArgumentError("charged pair must contain two distinct flavors"))
    channel === :P || channel === :S || throw(ArgumentError("channel must be canonical :P or :S"))
    kernel_pair in _VALID_KERNEL_PAIRS ||
        throw(ArgumentError("kernel_pair must be :K12 or :K45"))
    return ChargedRPAKernelSpec(
        meson,
        pair,
        channel,
        kernel_pair,
    )
end

"""
    charged_rpa_spec(meson; channel=:P)

Return the validated charged-RPA contract for a charged pion or kaon.

The standard flavor ordering is `(u,d)` for `pi_plus`, `(d,u)` for
`pi_minus`, `(u,s)` for `K_plus`, and `(s,u)` for `K_minus`.  Both charged
kaons use `K45`.  The normalization is fixed to
`D = 2K/(1 - 4K Pi)`; the neutral-matrix `1 - 2K Pi_matrix` form is equivalent
only after the charged-ladder relation `Pi_matrix = 2Pi` is applied.
"""
function charged_rpa_spec(
    meson::Symbol;
    channel::Symbol=:P,
)
    pair, kernel_pair = _charged_metadata(meson)
    canonical_channel = _canonical_channel(channel)
    return _validated_spec(
        meson,
        pair,
        canonical_channel,
        kernel_pair,
    )
end

"""Read the coupling selected by a charged-RPA spec from a full KMT kernel."""
@inline function charged_rpa_coupling(
    kernel::FullKMTInteraction,
    spec::ChargedRPAKernelSpec,
)
    return MesonInteractionKernel.charged_coupling(kernel, spec.kernel_pair, spec.channel)
end

@inline charged_rpa_coupling(spec::ChargedRPAKernelSpec, kernel::FullKMTInteraction) =
    charged_rpa_coupling(kernel, spec)

"""
    charged_rpa_inverse(spec, K_a, Pi_a)

Evaluate the charged single-channel RPA inverse propagator
`1 - 4 * K_a * Pi_a`.

`K_a` has units `fm^2`, `Pi_a` has units `fm^-2`, and the returned inverse is
dimensionless.  No pole regularization is applied: an exact zero is returned
as zero so callers can record the pole condition explicitly.
"""
function charged_rpa_inverse(
    spec::ChargedRPAKernelSpec,
    K_a::Real,
    Pi_a::Number,
)
    isfinite(K_a) || throw(ArgumentError("charged coupling K_a must be finite"))
    _finite_number(Pi_a) || throw(ArgumentError("charged polarization Pi_a must be finite"))
    denominator = one(Pi_a) - _DENOMINATOR_FACTOR * K_a * Pi_a
    _finite_number(denominator) || throw(ArgumentError("charged RPA inverse is not finite"))
    return denominator
end

"""
    charged_rpa_propagator(spec, K_a, Pi_a)

Evaluate the charged single-channel propagator
`2 * K_a / charged_rpa_inverse(spec, K_a, Pi_a)`.

The function deliberately throws `DomainError` at an exact pole instead of
adding the legacy epsilon regulator.  Retarded `eta`, pole finding, and
phase-shift continuation belong to later charged-RPA phases.
"""
function charged_rpa_propagator(
    spec::ChargedRPAKernelSpec,
    K_a::Real,
    Pi_a::Number,
)
    denominator = charged_rpa_inverse(spec, K_a, Pi_a)
    denominator == zero(denominator) &&
        throw(DomainError(denominator, "charged RPA propagator is singular at an exact pole"))
    result = _NUMERATOR_FACTOR * K_a / denominator
    _finite_number(result) || throw(ArgumentError("charged RPA propagator is not finite"))
    return result
end

end # module ChargedRPAKernel
