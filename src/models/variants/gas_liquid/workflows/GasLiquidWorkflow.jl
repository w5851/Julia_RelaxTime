module GasLiquidWorkflow

using StaticArrays

export solve_gas_liquid_point

"""最小单点 workflow：给定 (T, mu) 返回可观测量。"""
function solve_gas_liquid_point(T_fm::Real, mu_fm::Real; model_kind::Symbol=:GasLiquid)
    model = Main.Models.create_model(model_kind)
    x_state = Main.Models.solve_gap(model, T_fm, mu_fm)
    comp = Main.Models.omega_components(model, x_state, T_fm, mu_fm)
    rho = Main.Models.model_rho(model, x_state, mu_fm, T_fm)
    pressure, _, entropy, energy = Main.Models.model_thermo(model, x_state, mu_fm, T_fm)

    return (
        model_kind=model_kind,
        T_fm=float(T_fm),
        mu_fm=float(mu_fm),
        x_state=x_state,
        omega=float(comp.omega),
        pressure=float(pressure),
        rho=float(sum(rho) / 3),
        entropy=float(entropy),
        energy=float(energy),
        converged=true,
    )
end

end # module
