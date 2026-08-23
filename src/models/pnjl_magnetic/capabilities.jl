@inline pnjl_magnetic_capabilities(; eB_fm2::Real=0.0) = ModelCapabilities(
    supports_solve_gap=true,
    supports_model_thermo=true,
    # Nonzero magnetic fields expose net density through the dedicated
    # magnetic API, not the generic independent q/anti-q contract.
    supports_number_densities=abs(Float64(eB_fm2)) <= 1e-14,
)
