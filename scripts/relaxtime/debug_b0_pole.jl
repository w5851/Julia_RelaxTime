#!/usr/bin/env julia

# Temporary diagnostic for tilde_B0_k_zero pole handling.

push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src", "integration"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src", "relaxtime"))

using Constants_PNJL
using OneLoopIntegrals
using Printf

function main()
    s = 29.68421052631579
    k0 = sqrt(s)
    λs = (+k0, -k0)
    m = 500 / Constants_PNJL.ħc_MeV_fm
    mprime = 300 / Constants_PNJL.ħc_MeV_fm
    T = 150 / Constants_PNJL.ħc_MeV_fm
    println("s=", s, " k0=", k0, " m=", m, " m'=", mprime)

    for λ in λs
        Emin = m
        Emax = OneLoopIntegrals.energy_cutoff(m)
        denom_term = (λ^2 + m^2 - mprime^2) / 2
        E0 = -denom_term / λ
        println("\nλ=", λ, " Emin=", Emin, " Emax=", Emax, " E0=", E0)
        integrand(E) = OneLoopIntegrals.real_integrand_k_zero(:plus, λ, m, denom_term, 0.0, T, 0.5, 0.5, E)
        # sample integrand around E0
        for δ in (-1e-4, -1e-6, 0.0, 1e-6, 1e-4)
            E = E0 + δ
            if E > Emin && E < Emax
                val = integrand(E)
                println(@sprintf("E=%.10f δ=%+.1e -> integrand=%.5e", E, δ, val))
            end
        end
        # try quadgk directly
        try
            real_part, _ = OneLoopIntegrals.quadgk(integrand, Emin, Emax; rtol=1e-6, atol=0.0)
            println("quadgk real_part=", real_part)
        catch err
            println("quadgk error: ", err)
        end
        # run tilde_B0_k_zero
        try
            println("tilde_B0_k_zero: ", OneLoopIntegrals.tilde_B0_k_zero(:plus, λ, m, mprime, 0.0, T, 0.5, 0.5))
        catch err
            println("tilde_B0_k_zero error: ", err)
        end
    end
end

main()
