#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "MesonMassWorkflow.jl"))

using Main.Constants_PNJL: ħc_MeV_fm
using Main.MesonMassWorkflow: solve_gap_and_meson_point

struct Opts
    out_csv::String
    xi::Float64
    muB_MeV::Float64
    T_min_MeV::Float64
    T_max_MeV::Float64
    T_step_MeV::Float64
    p_num::Int
    t_num::Int
    max_iter::Int
end

function _usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime/export_workflow_selected_path.jl [options]")
    println("Options:")
    println("  --out <path>         Output CSV")
    println("  --xi <value>         Xi (default -0.3)")
    println("  --muB <MeV>          muB in MeV (default 0)")
    println("  --tmin <MeV>         T min (default 240)")
    println("  --tmax <MeV>         T max (default 260)")
    println("  --tstep <MeV>        T step (default 5)")
    println("  --p-num <int>        Momentum nodes (default 8)")
    println("  --t-num <int>        Angle nodes (default 4)")
    println("  --max-iter <int>     Iteration cap (default 25)")
end

function _parse(args::Vector{String})
    out_csv = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase", "residual_surface", "workflow_selected_path_xi-0p3.csv")
    xi = -0.3
    muB_MeV = 0.0
    tmin = 240.0
    tmax = 260.0
    tstep = 5.0
    p_num = 8
    t_num = 4
    max_iter = 25

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end
        if arg == "--out"
            out_csv = require_value()
        elseif arg == "--xi"
            xi = parse(Float64, require_value())
        elseif arg == "--muB"
            muB_MeV = parse(Float64, require_value())
        elseif arg == "--tmin"
            tmin = parse(Float64, require_value())
        elseif arg == "--tmax"
            tmax = parse(Float64, require_value())
        elseif arg == "--tstep"
            tstep = parse(Float64, require_value())
        elseif arg == "--p-num"
            p_num = parse(Int, require_value())
        elseif arg == "--t-num"
            t_num = parse(Int, require_value())
        elseif arg == "--max-iter"
            max_iter = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            _usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    return Opts(out_csv, xi, muB_MeV, tmin, tmax, tstep, p_num, t_num, max_iter)
end

function main()
    opts = _parse(ARGS)
    mkpath(dirname(opts.out_csv))

    mu_fm = (opts.muB_MeV / ħc_MeV_fm) / 3.0
    meson_seed_state = nothing
    mixed_seed_tracking_state = nothing

    open(opts.out_csv, "w") do io
        println(io, "T_MeV,xi,muB_MeV,pre_joint_mass,pre_joint_gamma,pre_joint_residual,post_joint_mass,post_joint_gamma,post_joint_residual,governance_selection_reason,joint_pair_selection_reason,selected_method")
        T = opts.T_min_MeV
        while T <= opts.T_max_MeV + 1e-9
            T_fm = T / ħc_MeV_fm
            res = solve_gap_and_meson_point(
                T_fm,
                mu_fm;
                xi=opts.xi,
                mesons=(:eta, :eta_prime),
                p_num=opts.p_num,
                t_num=opts.t_num,
                solver_kwargs=(iterations=opts.max_iter,),
                mass_kwargs=(iterations=opts.max_iter,),
                meson_seed_state=meson_seed_state,
                mixed_seed_tracking_state=mixed_seed_tracking_state,
            )

            eta_p = res.meson_results[:eta_prime]
            diag = eta_p.root_diagnostics
            println(io, string(
                T, ",", opts.xi, ",", opts.muB_MeV, ",",
                getproperty(diag, :pre_joint_mass), ",",
                getproperty(diag, :pre_joint_gamma), ",",
                getproperty(diag, :pre_joint_residual), ",",
                getproperty(diag, :post_joint_mass), ",",
                getproperty(diag, :post_joint_gamma), ",",
                getproperty(diag, :post_joint_residual), ",",
                getproperty(diag, :governance_selection_reason), ",",
                getproperty(diag, :joint_pair_selection_reason), ",",
                getproperty(diag, :selected_method),
            ))

            meson_seed_state = res.meson_seed_state
            mixed_seed_tracking_state = res.mixed_seed_tracking
            T += opts.T_step_MeV
        end
    end

    println("Wrote workflow selected path CSV: ", opts.out_csv)
end

main()
