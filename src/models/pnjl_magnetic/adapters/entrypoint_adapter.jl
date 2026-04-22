@inline pnjl_magnetic_entrypoint_adapter() = (;
    solve_gap=pnjl_magnetic_api_solve_gap,
    model_thermo=pnjl_magnetic_api_model_thermo,
    number_densities=pnjl_magnetic_api_number_densities,
)
