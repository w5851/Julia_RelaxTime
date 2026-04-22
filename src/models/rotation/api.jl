@inline rotation_api_solve_gap(model, T, mu_vec; kwargs...) = solve_gap(model, T, mu_vec; kwargs...)
@inline rotation_api_model_thermo(model, x_state, mu_vec, T; kwargs...) = model_thermo(model, x_state, mu_vec, T; kwargs...)
@inline rotation_api_number_densities(model, x_state, T, mu_vec; kwargs...) = number_densities(model, x_state, T, mu_vec; kwargs...)
