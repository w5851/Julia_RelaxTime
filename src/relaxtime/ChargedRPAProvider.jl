"""
    ChargedRPAProvider

Adapter contract for charged quark-antiquark polarization bubbles.

The provider keeps the ordered flavor pair from `ChargedRPAKernelSpec` and
calls the existing `PolarizationAniso` implementation with the same `A`/`B0`
regularization.  It returns a diagnostic `NamedTuple` containing the ordered
inputs, the complex polarization, and the continuation mode.  It does not
claim to implement a strict retarded analytic continuation, pole search, or
Beth-Uhlenbeck density.
"""
module ChargedRPAProvider

using ..ChargedRPAKernel: ChargedRPAKernelSpec
using ..PolarizationAniso: polarization_aniso, polarization_with_width

export charged_polarization

const _FLAVORS = (:u, :d, :s)
const _MODES = (:real_axis, :finite_width)

@inline function _finite_real(value, label::AbstractString)::Float64
    value isa Real || throw(ArgumentError("$(label) must be real"))
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError("$(label) must be finite"))
    return converted
end

@inline function _flavor_value(values, flavor::Symbol, label::AbstractString)::Float64
    all(candidate -> hasproperty(values, candidate), _FLAVORS) ||
        throw(ArgumentError("$(label) must provide fields :u, :d, and :s"))
    return _finite_real(getproperty(values, flavor), "$(label).$(flavor)")
end

@inline function _normalize_triplet(values, label::AbstractString)
    return (
        u=_flavor_value(values, :u, label),
        d=_flavor_value(values, :d, label),
        s=_flavor_value(values, :s, label),
    )
end

@inline function _normalize_thermo(thermo)
    hasproperty(thermo, :T) || throw(ArgumentError("thermo must provide :T"))
    hasproperty(thermo, :Φ) || throw(ArgumentError("thermo must provide :Φ"))
    hasproperty(thermo, :Φbar) || throw(ArgumentError("thermo must provide :Φbar"))
    return (
        T=_finite_real(thermo.T, "thermo.T"),
        Φ=_finite_real(thermo.Φ, "thermo.Φ"),
        Φbar=_finite_real(thermo.Φbar, "thermo.Φbar"),
        ξ=hasproperty(thermo, :ξ) ? _finite_real(thermo.ξ, "thermo.ξ") : 0.0,
    )
end

@inline function _ordered_num_s(spec::ChargedRPAKernelSpec, num_s_quark)
    if num_s_quark === nothing
        return (spec.pair[1] === :s || spec.pair[2] === :s) ? 1 : 0
    end
    num_s_quark isa Integer || throw(ArgumentError("num_s_quark must be an integer or nothing"))
    num_s_quark in (0, 1) || throw(ArgumentError("charged num_s_quark must be 0 or 1"))
    return Int(num_s_quark)
end

"""
    charged_polarization(spec, k0_inv_fm, q_inv_fm, masses, chemical_potentials,
                         thermo, A_values; kwargs...) -> NamedTuple

Evaluate one ordered charged bubble using the current `PolarizationAniso`
implementation.  `masses`, `chemical_potentials`, and `A_values` must each
provide `u`, `d`, and `s` fields in `fm^-1`, `fm^-1`, and `fm^-2`, respectively.
`thermo.T` is in `fm^-1`; `thermo.Φ`, `thermo.Φbar`, and optional `thermo.ξ`
are dimensionless.  The returned `value` is `Pi_a` in `fm^-2`.

`mode=:real_axis` calls `polarization_aniso` and requires `gamma_inv_fm=0`.
`mode=:finite_width` calls the existing real-axis finite-width algebra in
`polarization_with_width`; this mode is a compatibility proxy, not a proof of
the strict retarded analytic continuation.
"""
function charged_polarization(
    spec::ChargedRPAKernelSpec,
    k0_inv_fm::Real,
    q_inv_fm::Real,
    masses,
    chemical_potentials,
    thermo,
    A_values;
    gamma_inv_fm::Real=0.0,
    mode::Symbol=:real_axis,
    num_s_quark=nothing,
)
    mode in _MODES || throw(ArgumentError("unknown mode $(mode); use :real_axis or :finite_width"))
    k0 = _finite_real(k0_inv_fm, "k0_inv_fm")
    q = _finite_real(q_inv_fm, "q_inv_fm")
    q >= 0.0 || throw(ArgumentError("q_inv_fm must be non-negative"))
    gamma = _finite_real(gamma_inv_fm, "gamma_inv_fm")
    gamma >= 0.0 || throw(ArgumentError("gamma_inv_fm must be non-negative"))
    mode === :real_axis && gamma != 0.0 &&
        throw(ArgumentError("mode=:real_axis requires gamma_inv_fm=0"))

    mass_values = _normalize_triplet(masses, "masses")
    mu_values = _normalize_triplet(chemical_potentials, "chemical_potentials")
    A_triplet = _normalize_triplet(A_values, "A_values")
    thermo_values = _normalize_thermo(thermo)
    thermo_values.T > 0.0 || throw(ArgumentError("thermo.T must be positive"))
    num_s = _ordered_num_s(spec, num_s_quark)

    flavor1, flavor2 = spec.pair
    m1, m2 = mass_values[flavor1], mass_values[flavor2]
    μ1, μ2 = mu_values[flavor1], mu_values[flavor2]
    A1, A2 = A_triplet[flavor1], A_triplet[flavor2]

    Π_re, Π_im = if mode === :real_axis
        polarization_aniso(
            spec.channel,
            k0,
            q,
            m1,
            m2,
            μ1,
            μ2,
            thermo_values.T,
            thermo_values.Φ,
            thermo_values.Φbar,
            thermo_values.ξ,
            A1,
            A2,
            num_s,
        )
    else
        polarization_with_width(
            spec.channel,
            k0,
            gamma,
            q,
            m1,
            m2,
            μ1,
            μ2,
            thermo_values.T,
            thermo_values.Φ,
            thermo_values.Φbar,
            thermo_values.ξ,
            A1,
            A2,
            num_s,
        )
    end

    re = _finite_real(Π_re, "polarization real part")
    im = _finite_real(Π_im, "polarization imaginary part")
    return (
        spec=spec,
        meson=spec.meson,
        pair=spec.pair,
        channel=spec.channel,
        kernel_pair=spec.kernel_pair,
        value=ComplexF64(re, im),
        polarization_units=:fm_minus2,
        provider=:PolarizationAniso,
        mode=mode,
        retarded_convention=spec.retarded_convention,
        k0_inv_fm=k0,
        q_inv_fm=q,
        gamma_inv_fm=gamma,
        m1_inv_fm=m1,
        m2_inv_fm=m2,
        μ1_inv_fm=μ1,
        μ2_inv_fm=μ2,
        A1_inv_fm2=A1,
        A2_inv_fm2=A2,
        T_inv_fm=thermo_values.T,
        Φ=thermo_values.Φ,
        Φbar=thermo_values.Φbar,
        ξ=thermo_values.ξ,
        num_s_quark=num_s,
    )
end

end # module ChargedRPAProvider
