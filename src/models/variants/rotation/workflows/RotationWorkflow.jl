module RotationWorkflow

using ..RotationThermo: pressure_derivative_omega

export solve_rotation_point

function solve_rotation_point(T_fm::Real, mu_fm::Real; omega::Real=0.0, model_kind::Symbol=:Rotation)
    model = Main.Models.create_model(model_kind)
    x_state = Main.Models.solve_gap(model, T_fm, mu_fm; omega=omega)

    st = x_state isa Main.Models.MeanFieldState ? x_state : Main.Models.MeanFieldState(x_state)
    phi = st.phi[1]
    dP_domega = pressure_derivative_omega(phi, T_fm, mu_fm, omega, model.params)

    comp = Main.Models.omega_components(model, x_state, T_fm, mu_fm; omega=omega)
    rho = Main.Models.model_rho(model, x_state, mu_fm, T_fm; omega=omega)
    pressure, _, entropy, energy = Main.Models.model_thermo(model, x_state, mu_fm, T_fm; omega=omega)

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
