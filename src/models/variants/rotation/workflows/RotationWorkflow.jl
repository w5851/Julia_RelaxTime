module RotationWorkflow

using ..RotationThermo: pressure_derivative_omega
import ...Models

@inline _create_model(args...; kwargs...) = getfield(Models, :create_model)(args...; kwargs...)
@inline _solve_gap(args...; kwargs...) = getfield(Models, :solve_gap)(args...; kwargs...)
@inline _meanfield_state(x) = getfield(Models, :MeanFieldState)(x)
@inline _omega_components(args...; kwargs...) = getfield(Models, :omega_components)(args...; kwargs...)
@inline _model_rho(args...; kwargs...) = getfield(Models, :model_rho)(args...; kwargs...)
@inline _model_thermo(args...; kwargs...) = getfield(Models, :model_thermo)(args...; kwargs...)

export solve_rotation_point

function solve_rotation_point(T_fm::Real, mu_fm::Real; omega::Real=0.0, model_kind::Symbol=:Rotation)
    model = _create_model(model_kind)
    x_state = _solve_gap(model, T_fm, mu_fm; omega=omega)

    st = x_state isa getfield(Models, :MeanFieldState) ? x_state : _meanfield_state(x_state)
    phi = st.phi[1]
    dP_domega = pressure_derivative_omega(phi, T_fm, mu_fm, omega, model.params)

    comp = _omega_components(model, x_state, T_fm, mu_fm; omega=omega)
    rho = _model_rho(model, x_state, mu_fm, T_fm; omega=omega)
    pressure, _, entropy, energy = _model_thermo(model, x_state, mu_fm, T_fm; omega=omega)

    return (
        model_kind=model_kind,
        T_fm=float(T_fm),
        mu_fm=float(mu_fm),
        omega=float(omega),
        x_state=x_state,
        omega_potential=float(comp.omega),
        pressure=float(pressure),
        rho=float(sum(rho) / 3),
        entropy=float(entropy),
        energy=float(energy),
        dP_domega=float(dP_domega),
        converged=true,
    )
end

end # module
