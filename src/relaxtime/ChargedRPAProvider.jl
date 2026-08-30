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

using ..ChargedRPAKernel: ChargedRPAKernelSpec, charged_rpa_inverse
using ..PolarizationAniso: polarization_aniso, polarization_with_width

export charged_polarization
export charged_pole_residual, charged_mott_diagnostic

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

"""
    charged_pole_residual(spec, K_a, Pi_a) -> NamedTuple

Record the complex inverse-propagator residual
`Delta_a = 1 - c_den*K_a*Pi_a` at a candidate pole.  This helper does not
search for a root and does not add the legacy propagator epsilon; it is a
solver-independent diagnostic record for later pole providers.
"""
function charged_pole_residual(
    spec::ChargedRPAKernelSpec,
    K_a::Real,
    Pi_a::Number,
)
    inverse = charged_rpa_inverse(spec, K_a, Pi_a)
    return (
        meson=spec.meson,
        pair=spec.pair,
        channel=spec.channel,
        kernel_pair=spec.kernel_pair,
        coupling_fm2=Float64(K_a),
        polarization=Pi_a,
        polarization_units=:fm_minus2,
        inverse= inverse,
        residual_complex=inverse,
        residual_real=Float64(real(inverse)),
        residual_imag=Float64(imag(inverse)),
        residual_norm=Float64(abs(inverse)),
        normalization_source=spec.normalization_source,
        retarded_convention=spec.retarded_convention,
        root_search=false,
    )
end

"""
    charged_mott_diagnostic(spec, pole_mass_inv_fm, masses; kwargs...) -> NamedTuple

Classify a candidate charged pole against the constituent threshold
`m1+m2`, using the ordered pair in `spec`.  The returned `status` is
`:bound` below threshold, `:at_threshold` within `atol`, or `:continuum` above
threshold.  This is a kinematic record only; it does not assert that a pole
was found or that the current finite-width proxy is retarded.
"""
function charged_mott_diagnostic(
    spec::ChargedRPAKernelSpec,
    pole_mass_inv_fm::Real,
    masses;
    pole_gamma_inv_fm::Real=0.0,
    atol::Real=1e-6,
)
    pole_mass = _finite_real(pole_mass_inv_fm, "pole_mass_inv_fm")
    pole_mass >= 0.0 || throw(ArgumentError("pole_mass_inv_fm must be non-negative"))
    pole_gamma = _finite_real(pole_gamma_inv_fm, "pole_gamma_inv_fm")
    pole_gamma >= 0.0 || throw(ArgumentError("pole_gamma_inv_fm must be non-negative"))
    threshold_atol = _finite_real(atol, "atol")
    threshold_atol >= 0.0 || throw(ArgumentError("atol must be non-negative"))

    mass_values = _normalize_triplet(masses, "masses")
    flavor1, flavor2 = spec.pair
    m1 = mass_values[flavor1]
    m2 = mass_values[flavor2]
    threshold = m1 + m2
    gap = pole_mass - threshold
    status = if abs(gap) <= threshold_atol
        :at_threshold
    elseif gap < 0.0
        :bound
    else
        :continuum
    end

    return (
        meson=spec.meson,
        pair=spec.pair,
        pole_mass_inv_fm=pole_mass,
        pole_gamma_inv_fm=pole_gamma,
        threshold_inv_fm=threshold,
        threshold_gap_inv_fm=gap,
        threshold_source=:constituent_mass_sum,
        status=status,
        is_mott=(status === :at_threshold),
        atol=threshold_atol,
        retarded_convention=spec.retarded_convention,
        root_search=false,
    )
end

end # module ChargedRPAProvider
