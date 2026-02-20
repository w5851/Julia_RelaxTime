using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL

function run_magnetic_point(; T_MeV::Float64=150.0, mu_MeV::Float64=0.0, eB_MeV2::Float64=2.0e4)
    T_fm = T_MeV / Constants_PNJL.ħc_MeV_fm
    mu_fm = mu_MeV / Constants_PNJL.ħc_MeV_fm
    eB_fm2 = eB_MeV2 / (Constants_PNJL.ħc_MeV_fm^2)

    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
    mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)

    conf = default_magnetic_config(eB_fm2=eB_fm2)
    comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
    rho = calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)

    println("=== PNJL Magnetic Single Point ===")
    println("T = $(T_MeV) MeV, mu = $(mu_MeV) MeV, eB = $(eB_MeV2) MeV^2")
    println("n_max = $(comp.n_max)")
    println("G(B) = $(comp.G_B)")
    println("omega = $(comp.omega)")
    println("pressure = $(-comp.omega)")
    println("rho_u,d,s = $(Tuple(rho))")
end

run_magnetic_point()
