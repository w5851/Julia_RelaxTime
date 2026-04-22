@inline njl2_entrypoint_adapter() = (;
    solve_gap=njl2_api_solve_gap,
    model_thermo=njl2_api_model_thermo,
    number_densities=njl2_api_number_densities,
)
