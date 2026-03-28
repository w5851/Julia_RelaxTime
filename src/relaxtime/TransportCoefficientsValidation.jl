module TransportCoefficientsValidation

using Main.ValidationUtils: require_finite, require_positive_finite,
    require_nonnegative_finite, require_positive_integer, validate_grid_weight_pair

const _SPECIES_ALL = (:u, :d, :s, :ubar, :dbar, :sbar)
const _FLAVORS = (:u, :d, :s)
const _CONSERVED_CHARGES = (:B, :Q, :S)

export _validate_tau_namedtuple,
    _validate_quark_thermo_inputs,
    _validate_bulk_coeffs_isentropic,
    _validate_transport_inputs,
    _validate_transport_request_contract,
    _validate_conserved_charge_symbol,
    _validate_diffusion_background

@inline function _tau_for_species(species::Symbol, tau::NamedTuple)::Float64
    hasproperty(tau, species) || error("tau is missing :$species")
    return getproperty(tau, species)
end

@inline function _validate_tau_namedtuple(tau::NamedTuple)
    for sp in _SPECIES_ALL
        τ = _tau_for_species(sp, tau)
        require_nonnegative_finite("tau.:$sp", τ)
    end
    return nothing
end

@inline function _validate_quark_thermo_inputs(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    provider=nothing,
)
    require_positive_finite("thermo_params.T", thermo_params.T)
    require_finite("thermo_params.Φ", thermo_params.Φ)
    require_finite("thermo_params.Φbar", thermo_params.Φbar)

    check_mass = provider === nothing || !hasproperty(provider, :mass_for_species)
    check_mu = provider === nothing || !hasproperty(provider, :mu_for_species)
    for fl in _FLAVORS
        m = getproperty(quark_params.m, fl)
        μ = getproperty(quark_params.μ, fl)
        if check_mass
            require_nonnegative_finite("quark_params.m.$fl", m)
        end
        if check_mu
            require_finite("quark_params.μ.$fl", μ)
        end
    end
    return nothing
end

@inline function _validate_bulk_coeffs_isentropic(coeffs::NamedTuple)
    require_finite("bulk_coeffs_isentropic.v_n_sq", coeffs.v_n_sq)
    require_finite("bulk_coeffs_isentropic.dμB_dT_sigma", coeffs.dμB_dT_sigma)

    length(coeffs.masses) == 3 || throw(ArgumentError("bulk_coeffs_isentropic.masses must have length 3"))
    length(coeffs.dM_dT) == 3 || throw(ArgumentError("bulk_coeffs_isentropic.dM_dT must have length 3"))
    length(coeffs.dM_dμB) == 3 || throw(ArgumentError("bulk_coeffs_isentropic.dM_dμB must have length 3"))

    all(isfinite, coeffs.masses) || throw(ArgumentError("bulk_coeffs_isentropic.masses must be finite"))
    if !all(m -> m >= 0.0, coeffs.masses)
        @warn "bulk_coeffs_isentropic.masses contains negative value(s) — may indicate unphysical gap solution; proceeding with abs(m)" masses=coeffs.masses
    end
    all(isfinite, coeffs.dM_dT) || throw(ArgumentError("bulk_coeffs_isentropic.dM_dT must be finite"))
    all(isfinite, coeffs.dM_dμB) || throw(ArgumentError("bulk_coeffs_isentropic.dM_dμB must be finite"))
    return nothing
end

@inline function _validate_transport_inputs(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    tau::NamedTuple,
    config;
    provider=nothing,
    charges::Union{Nothing,NamedTuple}=nothing,
    bulk_coeffs_isentropic::Union{Nothing,NamedTuple}=nothing,
)
    _validate_quark_thermo_inputs(quark_params, thermo_params; provider=provider)
    _validate_tau_namedtuple(tau)

    require_positive_integer("integration p_nodes", config.p_nodes)
    require_positive_integer("integration cos_nodes", config.cos_nodes)
    require_positive_finite("integration p_max", config.p_max)
    validate_grid_weight_pair("TransportIntegrationConfig", "p_grid", config.p_grid, "p_w", config.p_w)
    validate_grid_weight_pair("TransportIntegrationConfig", "cos_grid", config.cos_grid, "cos_w", config.cos_w)

    if charges !== nothing
        require_finite("charges.u", charges.u)
        require_finite("charges.d", charges.d)
        require_finite("charges.s", charges.s)
    end

    if bulk_coeffs_isentropic !== nothing
        _validate_bulk_coeffs_isentropic(bulk_coeffs_isentropic)
    end

    return nothing
end

@inline function _validate_conserved_charge_symbol(charge::Symbol)::Symbol
    charge in _CONSERVED_CHARGES || error("Unknown conserved charge: $charge. Allowed charges: $(_CONSERVED_CHARGES)")
    return charge
end

function _validate_diffusion_background(densities::NamedTuple, pressure::Real, energy::Real)
    u = densities.u - densities.ubar
    d = densities.d - densities.dbar
    s = densities.s - densities.sbar
    charge_densities = (
        B=(u + d + s) / 3,
        Q=(2 * u - d - s) / 3,
        S=-s,
    )
    h = Float64(pressure) + Float64(energy)
    h > 0.0 || error("enthalpy density must be > 0")
    return charge_densities, h
end

@inline function _validate_transport_request_contract(
    densities::Union{Nothing,NamedTuple},
    pressure::Union{Nothing,Real},
    energy::Union{Nothing,Real},
)
    has_diffusion_background = densities !== nothing || pressure !== nothing || energy !== nothing
    if has_diffusion_background
        densities !== nothing || error("transport_coefficients requires densities when pressure/energy is provided")
        pressure !== nothing || error("transport_coefficients requires pressure when densities/energy is provided")
        energy !== nothing || error("transport_coefficients requires energy when densities/pressure is provided")
    end
    return has_diffusion_background
end

end # module TransportCoefficientsValidation
