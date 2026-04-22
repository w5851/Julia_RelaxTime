@inline njl_entrypoint_adapter() = (;
    solve_gap=njl_api_solve_gap,
    model_thermo=njl_api_model_thermo,
    number_densities=njl_api_number_densities,
)
