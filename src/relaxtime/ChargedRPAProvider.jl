"""
    ChargedRPAProvider

Adapter contract for ordered charged quark-antiquark polarization bubbles.

The provider keeps the ordered flavor pair from `ChargedRPAKernelSpec` and
calls the existing `PolarizationAniso` implementation with the same `A`/`B0`
regularization.  The default `:ordered_retarded` prescription keeps
`Pi_us`/`Pi_su` distinct with `num_s_quark=0`; the source-backed
`num_s_quark=1` average is available only as the explicit
`:legacy_symmetrized_B0` oracle.  The adapter does not locate poles, construct
phase shifts, or perform Beth-Uhlenbeck integration.
"""
module ChargedRPAProvider

using ..ChargedRPAKernel: ChargedRPAKernelSpec
using ..PolarizationAniso: polarization_aniso

export charged_polarization

const _FLAVORS = (:u, :d, :s)
const _PRESCRIPTIONS = (:ordered_retarded, :legacy_symmetrized_B0)

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

@inline function _num_s_quark(spec::ChargedRPAKernelSpec, prescription::Symbol)
    prescription in _PRESCRIPTIONS || throw(ArgumentError(
        "unknown prescription $(prescription); use :ordered_retarded or :legacy_symmetrized_B0",
    ))
    if prescription === :legacy_symmetrized_B0 &&
       (spec.pair[1] === :s || spec.pair[2] === :s)
        return 1
    end
    return 0
end

"""
    charged_polarization(spec, k0_inv_fm, q_inv_fm, masses, chemical_potentials,
                         thermo, A_values; kwargs...) -> NamedTuple

Evaluate one ordered charged bubble using the current `PolarizationAniso`
implementation.  `masses`, `chemical_potentials`, and `A_values` must each
provide `u`, `d`, and `s` fields in `fm^-1`, `fm^-1`, and `fm^-2`, respectively.
`thermo.T` is in `fm^-1`; `thermo.Φ`, `thermo.Φbar`, and optional `thermo.ξ`
are dimensionless.  The returned `value` is `Pi_a` in `fm^-2`.

`prescription=:ordered_retarded` is the strict-route real-axis candidate and
uses `num_s_quark=0`.  `:legacy_symmetrized_B0` reproduces the existing
Rehberg/Fortran/Cpp oracle for strange channels with `num_s_quark=1`.  Both are
real-axis adapters; finite-width and second-sheet pole semantics remain out of
scope.
"""
function charged_polarization(
    spec::ChargedRPAKernelSpec,
    k0_inv_fm::Real,
    q_inv_fm::Real,
    masses,
    chemical_potentials,
    thermo,
    A_values;
    prescription::Symbol=:ordered_retarded,
)
    k0 = _finite_real(k0_inv_fm, "k0_inv_fm")
    q = _finite_real(q_inv_fm, "q_inv_fm")
    q >= 0.0 || throw(ArgumentError("q_inv_fm must be non-negative"))

    mass_values = _normalize_triplet(masses, "masses")
    mu_values = _normalize_triplet(chemical_potentials, "chemical_potentials")
    A_triplet = _normalize_triplet(A_values, "A_values")
    thermo_values = _normalize_thermo(thermo)
    thermo_values.T > 0.0 || throw(ArgumentError("thermo.T must be positive"))
    num_s = _num_s_quark(spec, prescription)

    flavor1, flavor2 = spec.pair
    m1, m2 = mass_values[flavor1], mass_values[flavor2]
    μ1, μ2 = mu_values[flavor1], mu_values[flavor2]
    A1, A2 = A_triplet[flavor1], A_triplet[flavor2]

    Π_re, Π_im = polarization_aniso(
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
        prescription=prescription,
        analytic_scope=:real_axis,
        retarded_convention=spec.retarded_convention,
        k0_inv_fm=k0,
        q_inv_fm=q,
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
