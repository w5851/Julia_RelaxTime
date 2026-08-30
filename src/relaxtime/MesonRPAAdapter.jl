"""
    MesonRPAAdapter

Diagnostic adapter from the current `PolarizationAniso` implementation to the
source-aligned neutral three-flavor RPA algebra in `MesonRPA`.

The adapter evaluates three flavor-diagonal `(f, f)` quark bubbles and passes
them to the pure `(0, 3, 8)` matrix backend.  It is intentionally separate
from the legacy propagator and meson-density paths: it does not solve a gap
equation, change the mean-field background, locate poles, or add meson
pressure feedback.
"""
module MesonRPAAdapter

using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input
using ..AFieldBuilder: ensure_quark_params_has_A
using ..GaussLegendre: DEFAULT_COSΘ_NODES
using ..MesonInteractionKernel: FullKMTInteraction
using ..MesonRPA: neutral_polarization_matrix, neutral_rpa_inverse_matrix,
    neutral_rpa_propagator, neutral_rpa_determinant
using ..PolarizationAniso: polarization_aniso, polarization_with_width

export neutral_flavor_bubbles, neutral_rpa_from_quark_params

const _FLAVORS = (:u, :d, :s)

@inline function _normalize_channel(channel::Symbol)::Symbol
    if channel === :P || channel === :pseudoscalar || channel === :plus
        return :P
    elseif channel === :S || channel === :scalar || channel === :minus
        return :S
    end
    throw(ArgumentError("unknown polarization channel $(channel); use :P or :S"))
end

@inline function _finite_real(value, label::AbstractString)::Float64
    value isa Real || throw(ArgumentError("$(label) must be real"))
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError("$(label) must be finite"))
    return converted
end

function _coerce_real_triplet(values, label::AbstractString)::NamedTuple{(:u, :d, :s),Tuple{Float64,Float64,Float64}}
    all(flavor -> hasproperty(values, flavor), _FLAVORS) ||
        throw(ArgumentError("$(label) must provide fields :u, :d, and :s"))
    converted = ntuple(3) do index
        flavor = _FLAVORS[index]
        _finite_real(getproperty(values, flavor), "$(label).$(flavor)")
    end
    return (u=converted[1], d=converted[2], s=converted[3])
end

function _coerce_num_s_quark(values)::NamedTuple{(:u, :d, :s),Tuple{Int,Int,Int}}
    raw = if values isa NamedTuple || all(flavor -> hasproperty(values, flavor), _FLAVORS)
        all(flavor -> hasproperty(values, flavor), _FLAVORS) ||
            throw(ArgumentError("num_s_quark must provide fields :u, :d, and :s"))
        (getproperty(values, :u), getproperty(values, :d), getproperty(values, :s))
    elseif values isa Tuple || values isa AbstractVector
        length(values) == 3 || throw(ArgumentError("num_s_quark must contain exactly three entries (u, d, s)"))
        (values[1], values[2], values[3])
    else
        throw(ArgumentError("num_s_quark must be a three-entry tuple/vector or a NamedTuple with fields :u, :d, :s"))
    end

    converted = ntuple(3) do index
        value = raw[index]
        value isa Integer || throw(ArgumentError("num_s_quark.$(_FLAVORS[index]) must be an integer"))
        value >= 0 || throw(ArgumentError("num_s_quark.$(_FLAVORS[index]) must be non-negative"))
        Int(value)
    end
    return (u=converted[1], d=converted[2], s=converted[3])
end

function _normalize_parameters(quark_params, thermo_params;
                               ensure_A::Bool,
                               a_p_nodes::Integer,
                               a_p_max::Float64,
                               a_cos_nodes::Integer,
                               a_use_aniso::Bool,
                               warn_on_auto_A::Bool)
    a_p_nodes > 0 || throw(ArgumentError("a_p_nodes must be positive"))
    a_p_max > 0.0 || throw(ArgumentError("a_p_max must be positive"))
    a_cos_nodes > 0 || throw(ArgumentError("a_cos_nodes must be positive"))

    q_raw = normalize_quark_input(quark_params)
    hasproperty(q_raw, :m) || throw(ArgumentError("quark_params is missing :m"))
    hasproperty(q_raw, :μ) || throw(ArgumentError("quark_params is missing :μ"))
    masses = _coerce_real_triplet(q_raw.m, "quark_params.m")
    chemical_potentials = _coerce_real_triplet(q_raw.μ, "quark_params.mu")
    q = merge(q_raw, (m=masses, μ=chemical_potentials))

    t_raw = normalize_thermo_input(thermo_params)
    hasproperty(t_raw, :T) || throw(ArgumentError("thermo_params is missing :T"))
    hasproperty(t_raw, :Φ) || throw(ArgumentError("thermo_params is missing :Φ"))
    hasproperty(t_raw, :Φbar) || throw(ArgumentError("thermo_params is missing :Φbar"))
    temperature = _finite_real(t_raw.T, "thermo_params.T")
    phi = _finite_real(t_raw.Φ, "thermo_params.Phi")
    phi_bar = _finite_real(t_raw.Φbar, "thermo_params.PhiBar")
    xi = hasproperty(t_raw, :ξ) ? _finite_real(t_raw.ξ, "thermo_params.xi") : 0.0
    t = (T=temperature, Φ=phi, Φbar=phi_bar, ξ=xi)

    had_A = hasproperty(q, :A)
    if had_A
        A_values = _coerce_real_triplet(q.A, "quark_params.A")
        return merge(q, (A=A_values,)), t, false
    end

    ensure_A || throw(ArgumentError("quark_params is missing :A; set ensure_A=true or provide A.u/A.d/A.s"))
    q_with_A = ensure_quark_params_has_A(
        q,
        t;
        p_nodes=Int(a_p_nodes),
        p_max=a_p_max,
        cos_nodes=Int(a_cos_nodes),
        use_aniso=a_use_aniso,
        warn_on_auto=warn_on_auto_A,
    )
    A_values = _coerce_real_triplet(q_with_A.A, "quark_params.A")
    return merge(q_with_A, (A=A_values,)), t, true
end

@inline function _bubble_for_flavor(channel::Symbol,
                                    k0_inv_fm::Float64,
                                    k_norm_inv_fm::Float64,
                                    gamma_inv_fm::Float64,
                                    mass::Float64,
                                    chemical_potential::Float64,
                                    temperature::Float64,
                                    phi::Float64,
                                    phi_bar::Float64,
                                    xi::Float64,
                                    A_value::Float64,
                                    num_s_quark::Int,
                                    with_width::Bool,
                                    flavor::Symbol)::ComplexF64
    real_part, imag_part = if with_width
        polarization_with_width(
            channel,
            k0_inv_fm,
            gamma_inv_fm,
            k_norm_inv_fm,
            mass,
            mass,
            chemical_potential,
            chemical_potential,
            temperature,
            phi,
            phi_bar,
            xi,
            A_value,
            A_value,
            num_s_quark,
        )
    else
        polarization_aniso(
            channel,
            k0_inv_fm,
            k_norm_inv_fm,
            mass,
            mass,
            chemical_potential,
            chemical_potential,
            temperature,
            phi,
            phi_bar,
            xi,
            A_value,
            A_value,
            num_s_quark,
        )
    end
    re = _finite_real(real_part, "$(flavor) polarization real part")
    im = _finite_real(imag_part, "$(flavor) polarization imaginary part")
    return ComplexF64(re, im)
end

"""
    neutral_flavor_bubbles(channel, k0_inv_fm, k_norm_inv_fm,
                           quark_params, thermo_params; kwargs...) -> NamedTuple

Evaluate the three flavor-diagonal `(u,u)`, `(d,d)`, and `(s,s)` bubbles
using the current `PolarizationAniso` implementation.

All momentum, mass, chemical-potential, and temperature values use the
project's internal natural units (`fm^-1`).  `num_s_quark` defaults to
`(u=0, d=0, s=0)`, which is the explicit same-flavor convention; callers may
set individual entries to reproduce an existing legacy label.  The returned
NamedTuple has top-level `u`, `d`, and `s` fields, so it can be passed directly
to `neutral_polarization_matrix`, and also carries the normalized inputs and
adapter settings for diagnostics.
"""
function neutral_flavor_bubbles(channel::Symbol,
                                k0_inv_fm::Real,
                                k_norm_inv_fm::Real,
                                quark_params,
                                thermo_params;
                                gamma_inv_fm::Real=0.0,
                                with_width::Bool=false,
                                num_s_quark=(u=0, d=0, s=0),
                                ensure_A::Bool=true,
                                a_p_nodes::Integer=16,
                                a_p_max::Real=20.0,
                                a_cos_nodes::Integer=length(DEFAULT_COSΘ_NODES),
                                a_use_aniso::Bool=true,
                                warn_on_auto_A::Bool=true)
    channel_norm = _normalize_channel(channel)
    k0 = _finite_real(k0_inv_fm, "k0_inv_fm")
    k_norm = _finite_real(k_norm_inv_fm, "k_norm_inv_fm")
    k_norm >= 0.0 || throw(ArgumentError("k_norm_inv_fm must be non-negative"))
    gamma = _finite_real(gamma_inv_fm, "gamma_inv_fm")
    gamma >= 0.0 || throw(ArgumentError("gamma_inv_fm must be non-negative"))
    with_width || gamma == 0.0 ||
        throw(ArgumentError("gamma_inv_fm is nonzero but with_width=false"))
    num_s = _coerce_num_s_quark(num_s_quark)
    a_p_max_f64 = _finite_real(a_p_max, "a_p_max")
    q, t, auto_A = _normalize_parameters(
        quark_params,
        thermo_params;
        ensure_A=ensure_A,
        a_p_nodes=a_p_nodes,
        a_p_max=a_p_max_f64,
        a_cos_nodes=a_cos_nodes,
        a_use_aniso=a_use_aniso,
        warn_on_auto_A=warn_on_auto_A,
    )

    bubbles = ntuple(3) do index
        flavor = _FLAVORS[index]
        _bubble_for_flavor(
            channel_norm,
            k0,
            k_norm,
            gamma,
            q.m[flavor],
            q.μ[flavor],
            t.T,
            t.Φ,
            t.Φbar,
            t.ξ,
            q.A[flavor],
            num_s[flavor],
            with_width,
            flavor,
        )
    end

    return (
        u=bubbles[1],
        d=bubbles[2],
        s=bubbles[3],
        channel=channel_norm,
        k0_inv_fm=k0,
        k_norm_inv_fm=k_norm,
        gamma_inv_fm=gamma,
        with_width=with_width,
        num_s_quark=num_s,
        a_auto_built=auto_A,
        quark_params=q,
        thermo_params=t,
    )
end

"""
    neutral_rpa_from_quark_params(kernel, k0_inv_fm, k_norm_inv_fm,
                                  quark_params, thermo_params; kwargs...) -> NamedTuple

Evaluate the diagnostic flavor-bubble adapter and feed it through the pure
neutral RPA matrix backend.  The result includes the bubble record, the
polarization matrix, the RPA inverse matrix, the propagator, and its
determinant.  No mean-field or meson-thermodynamic feedback is performed.
"""
function neutral_rpa_from_quark_params(kernel::FullKMTInteraction,
                                       k0_inv_fm::Real,
                                       k_norm_inv_fm::Real,
                                       quark_params,
                                       thermo_params;
                                       channel::Symbol=:P,
                                       kwargs...)
    bubbles = neutral_flavor_bubbles(
        channel,
        k0_inv_fm,
        k_norm_inv_fm,
        quark_params,
        thermo_params;
        kwargs...,
    )
    polarization = neutral_polarization_matrix(bubbles)
    inverse_matrix = neutral_rpa_inverse_matrix(kernel, polarization; channel=bubbles.channel)
    propagator = neutral_rpa_propagator(kernel, polarization; channel=bubbles.channel)
    determinant = neutral_rpa_determinant(kernel, polarization; channel=bubbles.channel)
    return (
        kernel=kernel,
        bubbles=bubbles,
        polarization=polarization,
        inverse_matrix=inverse_matrix,
        propagator=propagator,
        determinant=determinant,
        channel=bubbles.channel,
    )
end

end # module MesonRPAAdapter
