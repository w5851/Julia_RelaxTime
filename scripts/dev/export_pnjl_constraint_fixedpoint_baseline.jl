#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_constraint_fixedpoints_v1.csv")
const HBARC_MEV_FM = 197.327

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
using .Models

const PNJL = Models.pnjl_module()

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_pnjl_constraint_fixedpoint_baseline.jl [--output <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return output
end

function _fixed_mu_points()
    return [
        (T_MeV=90.0, mu_MeV=0.0),
        (T_MeV=110.0, mu_MeV=20.0),
        (T_MeV=130.0, mu_MeV=40.0),
    ]
end

function _fixed_rho_points()
    return [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=130.0, rho_target=1.0),
    ]
end

function _to_fm_inv(x_mev::Real)
    return Float64(x_mev) / HBARC_MEV_FM
end

function _write_header(io)
    println(io, "kind,T_MeV,mu_MeV,rho_target,pnjl_converged,models_converged,pnjl_pressure,models_pressure,pnjl_rho_norm,models_rho_norm,pnjl_entropy,models_entropy,pnjl_energy,models_energy,pnjl_residual_norm,models_residual_norm")
end

function _write_row(io, kind::String, T_mev::Float64, mu_mev::Float64, rho_target::Float64, pnjl_res, models_res)
    @printf(
        io,
        "%s,%.6f,%.6f,%.6f,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
        kind,
        T_mev,
        mu_mev,
        rho_target,
        pnjl_res.converged ? 1 : 0,
        models_res.converged ? 1 : 0,
        pnjl_res.pressure,
        models_res.pressure,
        pnjl_res.rho_norm,
        models_res.rho_norm,
        pnjl_res.entropy,
        models_res.entropy,
        pnjl_res.energy,
        models_res.energy,
        pnjl_res.residual_norm,
        models_res.residual_norm,
    )
end

function main(args::Vector{String})
    output = parse_args(args)
    mkpath(dirname(output))

    model = Models.create_model(:PNJL)

    open(output, "w") do io
        _write_header(io)

        for pt in _fixed_mu_points()
            T_fm = _to_fm_inv(pt.T_MeV)
            mu_fm = _to_fm_inv(pt.mu_MeV)

            pnjl_res = PNJL.solve(PNJL.FixedMu(), T_fm, mu_fm; p_num=8, t_num=4, iterations=80)
            models_res = Models.solve_constraint(model, Models.FixedMu(), T_fm;
                μ_fm=mu_fm,
                seed_guess=PNJL.HADRON_SEED_5,
                p_num=8,
                t_num=4,
            )

            _write_row(io, "fixed_mu", pt.T_MeV, pt.mu_MeV, NaN, pnjl_res, models_res)
        end

        for pt in _fixed_rho_points()
            T_fm = _to_fm_inv(pt.T_MeV)

            pnjl_res = PNJL.solve(PNJL.FixedRho(pt.rho_target), T_fm; p_num=8, t_num=4, iterations=80)
            models_res = Models.solve_constraint(model, Models.FixedRho(pt.rho_target), T_fm;
                seed_guess=PNJL.HADRON_SEED_5,
                p_num=8,
                t_num=4,
            )

            _write_row(io, "fixed_rho", pt.T_MeV, NaN, pt.rho_target, pnjl_res, models_res)
        end
    end

    println("baseline exported to: " * output)
end

main(ARGS)
