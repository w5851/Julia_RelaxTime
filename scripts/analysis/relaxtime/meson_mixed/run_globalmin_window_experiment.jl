#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "MesonMassWorkflow.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .MesonMassWorkflow: solve_gap_and_meson_point, build_equilibrium_params
using Main.MesonMass: solve_meson_mass

struct Opts
    outdir::String
    libs::Vector{Symbol}
    T_min_MeV::Float64
    T_max_MeV::Float64
    T_step_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    max_iter::Int
    global_steps::Int
    jump_threshold::Float64
end

function _print_usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime/meson_mixed/run_globalmin_window_experiment.jl [options]")
    println("Options:")
    println("  --outdir <path>                Output directory")
    println("  --libs <csv>                   Library list: nlopt,metaheuristics,evolutionary,cmaes,optim_samin")
    println("  --tmin <MeV>                   T min (default 220)")
    println("  --tmax <MeV>                   T max (default 280)")
    println("  --tstep <MeV>                  T step (default 2)")
    println("  --muB <MeV>                    muB (default 0)")
    println("  --xi <value>                   xi (default -0.3)")
    println("  --p-num <int>                  momentum nodes (default 12)")
    println("  --t-num <int>                  angle nodes (default 6)")
    println("  --max-iter <int>               solver iterations (default 40)")
    println("  --global-steps <int>           global optimizer steps/budget (default 50)")
    println("  --jump-threshold <value>       jump threshold (default 0.25)")
    println("  -h, --help                     Show help")
end

function _parse_libs(raw::AbstractString)
    vals = Symbol[]
    for seg in split(raw, ',')
        s = Symbol(lowercase(strip(seg)))
        s in (:nlopt, :metaheuristics, :evolutionary, :cmaes, :optim_samin) || throw(ArgumentError("unsupported lib: $(seg)"))
        push!(vals, s)
    end
    isempty(vals) && throw(ArgumentError("--libs cannot be empty"))
    return unique(vals)
end

function _parse_args(args::Vector{String})
    outdir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase", "mixed_window_globalmin")
    libs = Symbol[:nlopt, :metaheuristics, :evolutionary]
    T_min_MeV = 220.0
    T_max_MeV = 280.0
    T_step_MeV = 2.0
    muB_MeV = 0.0
    xi = -0.3
    p_num = 12
    t_num = 6
    max_iter = 40
    global_steps = 50
    jump_threshold = 0.25

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end

        if arg == "--outdir"
            outdir = require_value()
        elseif arg == "--libs"
            libs = _parse_libs(require_value())
        elseif arg == "--tmin"
            T_min_MeV = parse(Float64, require_value())
        elseif arg == "--tmax"
            T_max_MeV = parse(Float64, require_value())
        elseif arg == "--tstep"
            T_step_MeV = parse(Float64, require_value())
        elseif arg == "--muB"
            muB_MeV = parse(Float64, require_value())
        elseif arg == "--xi"
            xi = parse(Float64, require_value())
        elseif arg == "--p-num"
            p_num = parse(Int, require_value())
        elseif arg == "--t-num"
            t_num = parse(Int, require_value())
        elseif arg == "--max-iter"
            max_iter = parse(Int, require_value())
        elseif arg == "--global-steps"
            global_steps = parse(Int, require_value())
        elseif arg == "--jump-threshold"
            jump_threshold = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            _print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    T_step_MeV > 0 || throw(ArgumentError("T_step_MeV must be positive"))
    T_max_MeV >= T_min_MeV || throw(ArgumentError("T_max_MeV must be >= T_min_MeV"))

    return Opts(
        String(outdir),
        libs,
        T_min_MeV,
        T_max_MeV,
        T_step_MeV,
        muB_MeV,
        xi,
        p_num,
        t_num,
        max_iter,
        global_steps,
        jump_threshold,
    )
end

@inline function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    return string(x)
end

const _NLOPT_AVAILABLE = let ok = false
    try
        @eval import NLopt
        ok = true
    catch
        ok = false
    end
    ok
end

const _META_AVAILABLE = let ok = false
    try
        @eval import Metaheuristics
        ok = true
    catch
        ok = false
    end
    ok
end

const _EVO_AVAILABLE = let ok = false
    try
        @eval import Evolutionary
        ok = true
    catch
        ok = false
    end
    ok
end

const _CMAES_AVAILABLE = let ok = false
    try
        @eval import CMAEvolutionStrategy
        ok = true
    catch
        ok = false
    end
    ok
end

const _OPTIM_AVAILABLE = let ok = false
    try
        @eval import Optim
        ok = true
    catch
        ok = false
    end
    ok
end

function _default_seed()
    return Float64[2.8, 2.3]
end

@inline function _clamp_seed(x::Vector{Float64}, lower::Vector{Float64}, upper::Vector{Float64})
    return Float64[
        clamp(Float64(x[1]), lower[1], upper[1]),
        clamp(Float64(x[2]), lower[2], upper[2]),
    ]
end

function _optimize_seed_nlopt(objective::Function, lower::Vector{Float64}, upper::Vector{Float64}, steps::Int)
    _NLOPT_AVAILABLE || return nothing
    dim = 2
    opt = NLopt.Opt(:GN_DIRECT_L, dim)
    NLopt.lower_bounds!(opt, lower)
    NLopt.upper_bounds!(opt, upper)
    NLopt.min_objective!(opt, (x, _grad) -> objective(Float64[x[1], x[2]]))
    NLopt.maxeval!(opt, max(steps, 10))
    minf, minx, _ret = NLopt.optimize(opt, _default_seed())
    isfinite(minf) || return nothing
    return _clamp_seed(Float64[minx[1], minx[2]], lower, upper)
end

function _optimize_seed_metaheuristics(objective::Function, lower::Vector{Float64}, upper::Vector{Float64}, steps::Int)
    _META_AVAILABLE || return nothing
    bounds = Metaheuristics.boxconstraints(lower, upper)
    algo = Metaheuristics.DE(options=Metaheuristics.Options(iterations=max(steps, 10)))
    result = Metaheuristics.optimize(objective, bounds, algo)
    minx = Metaheuristics.minimizer(result)
    return _clamp_seed(Float64[minx[1], minx[2]], lower, upper)
end

function _optimize_seed_evolutionary(objective::Function, lower::Vector{Float64}, upper::Vector{Float64}, steps::Int)
    _EVO_AVAILABLE || return nothing
    wrapped = function (x)
        xc = _clamp_seed(Float64[x[1], x[2]], lower, upper)
        return objective(xc)
    end
    opts = Evolutionary.Options(iterations=max(steps, 10), show_trace=false)
    result = Evolutionary.optimize(wrapped, _default_seed(), Evolutionary.CMAES(), opts)
    minx = Evolutionary.minimizer(result)
    return _clamp_seed(Float64[minx[1], minx[2]], lower, upper)
end

function _optimize_seed_cmaes(objective::Function, lower::Vector{Float64}, upper::Vector{Float64}, steps::Int)
    _CMAES_AVAILABLE || return nothing
    wrapped = function (x)
        xc = _clamp_seed(Float64[x[1], x[2]], lower, upper)
        return objective(xc)
    end
    x0 = _default_seed()
    sigma0 = 0.2
    algo = CMAEvolutionStrategy.minimize(
        wrapped,
        x0,
        sigma0;
        maxfevals=max(steps * 20, 80),
        verbosity=0,
    )
    return _clamp_seed(Float64[algo.xbest[1], algo.xbest[2]], lower, upper)
end

function _optimize_seed_optim_samin(objective::Function, lower::Vector{Float64}, upper::Vector{Float64}, steps::Int)
    _OPTIM_AVAILABLE || return nothing
    x0 = _default_seed()
    res = Optim.optimize(
        objective,
        lower,
        upper,
        x0,
        Optim.SAMIN();
        iterations=max(steps * 20, 200),
    )
    minx = Optim.minimizer(res)
    return _clamp_seed(Float64[minx[1], minx[2]], lower, upper)
end

function _optimize_seed(lib::Symbol, objective::Function, lower::Vector{Float64}, upper::Vector{Float64}, steps::Int)
    try
        if lib === :nlopt
            return _optimize_seed_nlopt(objective, lower, upper, steps)
        elseif lib === :metaheuristics
            return _optimize_seed_metaheuristics(objective, lower, upper, steps)
        elseif lib === :evolutionary
            return _optimize_seed_evolutionary(objective, lower, upper, steps)
        elseif lib === :cmaes
            return _optimize_seed_cmaes(objective, lower, upper, steps)
        elseif lib === :optim_samin
            return _optimize_seed_optim_samin(objective, lower, upper, steps)
        end
    catch
        return nothing
    end
    return nothing
end

function _build_jump_summary(rows::Vector{Dict{String,Any}}, threshold::Float64)
    grouped = Dict{String,Vector{Dict{String,Any}}}()
    for r in rows
        k = String(r["lib"])
        get!(grouped, k, Dict{String,Any}[])
        push!(grouped[k], r)
    end

    out = Dict{String,Any}[]
    observables = (
        "M_eta",
        "M_eta_prime",
        "Gamma_eta",
        "Gamma_eta_prime",
    )
    for lib in sort(collect(keys(grouped)))
        rr = grouped[lib]
        sort!(rr; by=r -> Float64(r["T_MeV"]))
        s = Dict{String,Any}("lib" => lib)
        for obs in observables
            cnt = 0
            maxjump = 0.0
            for i in 2:length(rr)
                a = Float64(rr[i - 1][obs])
                b = Float64(rr[i][obs])
                dj = abs(b - a)
                if dj >= threshold
                    cnt += 1
                end
                if dj > maxjump
                    maxjump = dj
                end
            end
            s["jump_count_" * obs] = cnt
            s["max_jump_" * obs] = maxjump
        end
        push!(out, s)
    end
    return out
end

function _write_main_csv(path::String, rows::Vector{Dict{String,Any}})
    cols = [
        "lib", "T_MeV", "muB_MeV", "xi",
        "global_seed_mass", "global_seed_gamma", "global_obj",
        "M_eta", "M_eta_prime", "Gamma_eta", "Gamma_eta_prime",
        "residual_eta", "residual_eta_prime",
        "selected_method_eta", "selected_method_eta_prime",
        "status", "error_message",
    ]
    open(path, "w") do io
        println(io, "# script=scripts/analysis/relaxtime/meson_mixed/run_globalmin_window_experiment.jl")
        println(io, "# window.xi=-0.3")
        println(io, "# window.T_MeV=[220,280]")
        println(io, join(cols, ','))
        for r in rows
            println(io, join([_fmt(get(r, c, "")) for c in cols], ','))
        end
    end
end

function _write_jump_csv(path::String, summary::Vector{Dict{String,Any}})
    cols = [
        "lib",
        "jump_count_M_eta", "max_jump_M_eta",
        "jump_count_M_eta_prime", "max_jump_M_eta_prime",
        "jump_count_Gamma_eta", "max_jump_Gamma_eta",
        "jump_count_Gamma_eta_prime", "max_jump_Gamma_eta_prime",
    ]
    open(path, "w") do io
        println(io, join(cols, ','))
        for r in summary
            println(io, join([_fmt(get(r, c, "")) for c in cols], ','))
        end
    end
end

function _run_plot(csv_path::String, outdir::String)
    script = joinpath(PROJECT_ROOT, "scripts", "plot_scan_csv.py")
    mass_dir = joinpath(outdir, "fig_mass")
    gamma_dir = joinpath(outdir, "fig_gamma")
    mkpath(mass_dir)
    mkpath(gamma_dir)

    run(`python $script --mode lines --csv $csv_path --x T_MeV --ys M_eta,M_eta_prime --group lib --out-dir $mass_dir`)
    run(`python $script --mode lines --csv $csv_path --x T_MeV --ys Gamma_eta,Gamma_eta_prime --group lib --out-dir $gamma_dir`)
end

function _eta_prime_objective(qp, tp, max_iter::Int)
    qp_nt = (
        m=(u=Float64(qp.m.u), d=Float64(qp.m.d), s=Float64(qp.m.s)),
        μ=(u=Float64(qp.μ.u), d=Float64(qp.μ.d), s=Float64(qp.μ.s)),
    )
    tp_nt = (
        T=Float64(tp.T),
        Φ=Float64(tp.Φ),
        Φbar=Float64(tp.Φbar),
        ξ=Float64(tp.ξ),
    )
    return function (x::Vector{Float64})
        m = Float64(x[1])
        g = max(Float64(x[2]), 0.0)
        res = try
            solve_meson_mass(
                :eta_prime,
                qp_nt,
                tp_nt;
                k_norm=0.0,
                initial_mass=m,
                initial_gamma=g,
                method=:trust_region,
                iterations=max_iter,
            )
        catch
            nothing
        end
        res === nothing && return Inf
        v = Float64(res.residual_norm)
        return isfinite(v) ? v : Inf
    end
end

function _run_single_point(T_MeV::Float64, muB_MeV::Float64, xi::Float64, seed::Vector{Float64}, opts::Opts)
    T_fm = T_MeV / ħc_MeV_fm
    mu_fm = (muB_MeV / ħc_MeV_fm) / 3.0
    base = Main.EquilibriumFacade.solve_equilibrium_backend(
        T_fm,
        mu_fm;
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        solver_kwargs=(; iterations=opts.max_iter),
    )
    params = build_equilibrium_params(base, T_fm, mu_fm; xi=xi)
    qp = params.quark_params
    tp = params.thermo_params

    lower = Float64[2.2, 0.0]
    upper = Float64[3.3, 3.0]
    objective = _eta_prime_objective(qp, tp, opts.max_iter)
    best_seed = _clamp_seed(seed, lower, upper)
    best_obj = objective(best_seed)

    for lib in opts.libs
        candidate = _optimize_seed(lib, objective, lower, upper, opts.global_steps)
        candidate === nothing && continue
        obj = objective(candidate)
        if isfinite(obj) && obj < best_obj
            best_obj = obj
            best_seed = candidate
        end
    end

    seed_track = Dict{Symbol,Vector{Float64}}(:eta_plus => best_seed)
    return solve_gap_and_meson_point(
        T_fm,
        mu_fm;
        xi=xi,
        mesons=(:eta, :eta_prime),
        meson_seed_state=Dict{Symbol,Vector{Float64}}(:eta_prime => best_seed),
        mixed_seed_tracking_state=seed_track,
        mixed_branch_align=:strict_sign_binding,
        p_num=opts.p_num,
        t_num=opts.t_num,
        solver_kwargs=(; iterations=opts.max_iter),
        mass_kwargs=(; iterations=opts.max_iter),
        force_global_fallback=true,
    ), best_seed, best_obj
end

function main()
    opts = _parse_args(ARGS)
    mkpath(opts.outdir)

    rows = Dict{String,Any}[]
    for lib in opts.libs
        seed = _default_seed()
        T = opts.T_min_MeV
        while T <= opts.T_max_MeV + 1e-9
            row = Dict{String,Any}(
                "lib" => String(lib),
                "T_MeV" => T,
                "muB_MeV" => opts.muB_MeV,
                "xi" => opts.xi,
                "global_seed_mass" => NaN,
                "global_seed_gamma" => NaN,
                "global_obj" => NaN,
                "M_eta" => NaN,
                "M_eta_prime" => NaN,
                "Gamma_eta" => NaN,
                "Gamma_eta_prime" => NaN,
                "residual_eta" => NaN,
                "residual_eta_prime" => NaN,
                "selected_method_eta" => "",
                "selected_method_eta_prime" => "",
                "status" => "ok",
                "error_message" => "",
            )
            try
                run_opts = Opts(opts.outdir, [lib], opts.T_min_MeV, opts.T_max_MeV, opts.T_step_MeV, opts.muB_MeV, opts.xi, opts.p_num, opts.t_num, opts.max_iter, opts.global_steps, opts.jump_threshold)
                res, best_seed, best_obj = _run_single_point(T, opts.muB_MeV, opts.xi, seed, run_opts)
                eta = res.meson_results[:eta]
                etap = res.meson_results[:eta_prime]
                row["global_seed_mass"] = best_seed[1]
                row["global_seed_gamma"] = best_seed[2]
                row["global_obj"] = best_obj
                row["M_eta"] = Float64(eta.mass)
                row["M_eta_prime"] = Float64(etap.mass)
                row["Gamma_eta"] = Float64(eta.gamma)
                row["Gamma_eta_prime"] = Float64(etap.gamma)
                row["residual_eta"] = Float64(eta.residual)
                row["residual_eta_prime"] = Float64(etap.residual)
                row["selected_method_eta"] = String(eta.root_diagnostics.selected_method)
                row["selected_method_eta_prime"] = String(etap.root_diagnostics.selected_method)
                seed = Float64[Float64(etap.mass), max(Float64(etap.gamma), 0.0)]
            catch e
                row["status"] = "error"
                row["error_message"] = sprint(showerror, e)
            end
            push!(rows, row)
            T += opts.T_step_MeV
        end
    end

    csv_path = joinpath(opts.outdir, "window_eta_etap.csv")
    jump_path = joinpath(opts.outdir, "window_jump_summary.csv")
    _write_main_csv(csv_path, rows)
    summary = _build_jump_summary(rows, opts.jump_threshold)
    _write_jump_csv(jump_path, summary)
    _run_plot(csv_path, opts.outdir)

    println("Wrote window data: ", csv_path)
    println("Wrote jump summary: ", jump_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
