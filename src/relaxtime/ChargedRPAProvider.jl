"""
    ChargedRPAProvider

Adapter contract for ordered charged quark-antiquark polarization bubbles.

The default `:ordered_retarded` prescription evaluates the ordered bubble at
`z=omega+i*eta` using the explicit complex `B0_retarded` backend. Historical
real-axis `B0` behavior remains available through `:ordered_legacy_B0`, while
the source-backed `num_s_quark=1` average remains the explicit
`:legacy_symmetrized_B0` oracle. The adapter does not locate poles, construct
phase shifts, or perform Beth-Uhlenbeck integration.
"""
module ChargedRPAProvider

using ..ChargedRPAKernel: ChargedRPAKernelSpec
using ..OneLoopIntegrals: B0_pv_cut, B0_retarded, EPS_SEGMENT
using ..PolarizationAniso: polarization_aniso
using Main.Constants_PNJL: N_color

export charged_polarization

const _FLAVORS = (:u, :d, :s)
const _PRESCRIPTIONS = (:ordered_retarded, :ordered_pv_cut, :ordered_legacy_B0, :legacy_symmetrized_B0)

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
        "unknown prescription $(prescription); use :ordered_retarded, :ordered_pv_cut, :ordered_legacy_B0, or :legacy_symmetrized_B0",
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

`prescription=:ordered_retarded` is the strict upper-half-plane probe and uses
the explicit `eta_inv_fm` and `energy_nodes` controls. `:ordered_pv_cut` is the
real-axis principal-value plus analytic-cut adapter with the retarded
`e^{-iωt}` imaginary-part sign. Both ordered strict prescriptions currently
support only `xi=0`. `:ordered_legacy_B0` preserves the earlier ordered adapter
with `num_s_quark=0`, while `:legacy_symmetrized_B0` reproduces the existing
Rehberg/Fortran/Cpp strange-channel oracle with `num_s_quark=1`. Finite-width
and second-sheet pole semantics remain out of scope.
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
    eta_inv_fm::Real=1.0e-3,
    energy_nodes::Integer=128,
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

    eta = _finite_real(eta_inv_fm, "eta_inv_fm")

    value, provider, analytic_scope, effective_eta, effective_nodes = if prescription === :ordered_retarded
        abs(thermo_values.ξ) <= EPS_SEGMENT || throw(ArgumentError(
            ":ordered_retarded currently supports only thermo.ξ=0; use a legacy prescription for diagnostic anisotropic evaluation",
        ))
        eta > 0.0 || throw(ArgumentError("eta_inv_fm must be positive for :ordered_retarded"))
        energy_nodes >= 4 || throw(ArgumentError("energy_nodes must be at least 4 for :ordered_retarded"))
        λ = ComplexF64(k0 + μ1 - μ2, eta)
        B0_value = B0_retarded(
            real(λ), q, m1, μ1, m2, μ2, thermo_values.T;
            Φ=thermo_values.Φ,
            Φbar=thermo_values.Φbar,
            eta_inv_fm=eta,
            energy_nodes=energy_nodes,
        )
        mass_term = spec.channel === :P ? (m1 - m2)^2 : (m1 + m2)^2
        prefactor = q^2 - λ^2 + mass_term
        Π = (-N_color / (8π^2)) * (A1 + A2 + prefactor * B0_value)
        (ComplexF64(Π), :OneLoopIntegralsRetarded, :upper_half_plane_probe, eta, Int(energy_nodes))
    elseif prescription === :ordered_pv_cut
        abs(thermo_values.ξ) <= EPS_SEGMENT || throw(ArgumentError(
            ":ordered_pv_cut currently supports only thermo.ξ=0; use a legacy prescription for diagnostic anisotropic evaluation",
        ))
        B0_value = B0_pv_cut(
            k0 + μ1 - μ2, q, m1, μ1, m2, μ2, thermo_values.T;
            Φ=thermo_values.Φ,
            Φbar=thermo_values.Φbar,
        )
        mass_term = spec.channel === :P ? (m1 - m2)^2 : (m1 + m2)^2
        prefactor = q^2 - (k0 + μ1 - μ2)^2 + mass_term
        Π = (-N_color / (8π^2)) * (A1 + A2 + prefactor * B0_value)
        (ComplexF64(Π), :OneLoopIntegralsPV, :real_axis_pv_cut, 0.0, 0)
    else
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
        Π = ComplexF64(
            _finite_real(Π_re, "polarization real part"),
            _finite_real(Π_im, "polarization imaginary part"),
        )
        (Π, :PolarizationAniso, :real_axis_legacy, 0.0, 0)
    end

    isfinite(real(value)) && isfinite(imag(value)) ||
        throw(ArgumentError("polarization value must be finite"))
    return (
        spec=spec,
        meson=spec.meson,
        pair=spec.pair,
        channel=spec.channel,
        kernel_pair=spec.kernel_pair,
        value=value,
        polarization_units=:fm_minus2,
        provider=provider,
        prescription=prescription,
        analytic_scope=analytic_scope,
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
        eta_inv_fm=effective_eta,
        energy_nodes=effective_nodes,
    )
end

end # module ChargedRPAProvider
