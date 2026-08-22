@inline function pnjl_magnetic_capabilities(; eB_fm2::Real)
    MagneticThermodynamics.validate_magnetic_eB(eB_fm2)
    return ModelCapabilities(
        supports_solve_gap=true,
        supports_model_thermo=true,
        # Magnetic fields expose net density through the dedicated API, not
        # the generic independent q/anti-q contract.
        supports_number_densities=false,
    )
end
