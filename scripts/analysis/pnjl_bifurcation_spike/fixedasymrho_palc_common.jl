module FixedAsymRhoPALCSpike

using Dates
using LinearAlgebra
using Printf

try
    @eval import BifurcationKit
    @eval import BifurcationKit: @optic
catch err
    error("BifurcationKit is not available. Run with --project=scripts/analysis/pnjl_bifurcation_spike and instantiate that environment. Original error: $(err)")
end

const U = Main.FixedMuPALCSpike
const HBARC_MEV_FM = U.HBARC_MEV_FM
const DEFAULT_OUTPUT_REL = joinpath("data", "outputs", "results", "analysis", "palc_fixedasymrho_spike")
const DEFAULT_T_MEV = 120.0
const DEFAULT_RHO_MIN = 1.05
const DEFAULT_RHO_MAX = 2.15
const DEFAULT_RHO_STEP = 0.05
const DEFAULT_ASYM_UD_RATIO = 0.876
const DEFAULT_ASYM_S_TARGET = 0.0
const DEFAULT_XI = 0.0
const DEFAULT_P_NUM = 8
const DEFAULT_T_NUM = 4
const DEFAULT_DS_RHO = 0.02
const DEFAULT_DSMAX_RHO = 0.08
const DEFAULT_MAX_STEPS = 120
const DEFAULT_BRANCH_JUMP_TOL = 0.5
const DEFAULT_PRESSURE_GAP_TOL = 1e-3
const PALC_NEWTON_TOL = 1e-7
const PALC_NEWTON_MAX_ITERATIONS = 10

struct CliConfig
    T_MeV::Float64
    rho_min::Float64
    rho_max::Float64
    rho_step::Float64
    asym_ud_ratio_target::Float64
    asym_s_target::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    output_dir::Union{Nothing, String}
    run_reference::Bool
    model_kind::Symbol
    ds_rho::Float64
    dsmax_rho::Float64
    max_steps::Int
    branch_jump_tol::Float64
    trho_reverse_rho::Bool
    pressure_gap_tol::Float64
end

Base.@kwdef struct EffectiveConfig
    repo_root::String
    output_dir::String
    run_id::String
    T_MeV::Float64
    rho_min::Float64
    rho_max::Float64
    rho_step::Float64
    rho_start::Float64
    asym_ud_ratio_target::Float64
    asym_s_target::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    run_reference::Bool
    model_kind::Symbol
    ds_rho::Float64
    dsmax_rho::Float64
    max_steps::Int
    branch_jump_tol::Float64
    trho_reverse_rho::Bool
    pressure_gap_tol::Float64
end

function _usage()
    return """
    Usage:
      julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedasymrho_palc.jl [options]

    Options:
      --T-MeV=<value>                    Fixed temperature in MeV. Default: 120.
      --rho-min=<value>                  Lower rho/rho0 continuation bound. Default: 1.05.
      --rho-max=<value>                  Upper rho/rho0 continuation bound. Default: 2.15.
      --rho-step=<value>                 Reference Trho scan step. Default: 0.05.
      --asym-ud-ratio-target=<value>     Target rho_u/rho_d. Default: 0.876.
      --asym-s-target=<value>            Target rho_s in fm^-3. Default: 0.
      --xi=<value>                       Anisotropy parameter. Default: 0.
      --p-num=<int>                      Momentum nodes. Default: 8.
      --t-num=<int>                      Theta nodes. Default: 4.
      --output-dir=<path>                Output directory. Default: data/outputs/results/analysis/palc_fixedasymrho_spike/<run_id>.
      --run-reference=<bool>             Run run_trho_scan reference.
      --trho-reverse-rho=<bool>          Run reference scan from high rho to low rho. Default: false.
      --ds-rho=<value>                   PALC initial pseudo-arclength scale in rho/rho0 units. Default: 0.02.
      --dsmax-rho=<value>                PALC max pseudo-arclength scale in rho/rho0 units. Default: 0.08.
      --max-steps=<int>                  PALC max continuation steps. Default: 120.
      --pressure-gap-tol=<value>         Pressure gap tolerance for branch comparison. Default: 1e-3.
    """
end

function _parse_cli(args::Vector{String}; default_run_reference::Bool)
    raw = Dict{String, String}()
    for arg in args
        arg in ("-h", "--help") && (println(_usage()); exit(0))
        startswith(arg, "--") || throw(ArgumentError("unexpected positional argument: $(arg)"))
        body = arg[3:end]
        parts = split(body, "="; limit=2)
        length(parts) == 2 || throw(ArgumentError("expected --key=value, got $(arg)"))
        raw[replace(lowercase(parts[1]), "-" => "_")] = parts[2]
    end

    get_float(key, default) = haskey(raw, key) ? parse(Float64, raw[key]) : default
    get_int(key, default) = haskey(raw, key) ? parse(Int, raw[key]) : default
    function get_bool(key, default)
        haskey(raw, key) || return default
        value = lowercase(strip(raw[key]))
        value in ("1", "true", "yes", "y") && return true
        value in ("0", "false", "no", "n") && return false
        throw(ArgumentError("$(key) must be true or false, got $(raw[key])"))
    end

    return CliConfig(
        get_float("t_mev", DEFAULT_T_MEV),
        get_float("rho_min", DEFAULT_RHO_MIN),
        get_float("rho_max", DEFAULT_RHO_MAX),
        get_float("rho_step", DEFAULT_RHO_STEP),
        get_float("asym_ud_ratio_target", DEFAULT_ASYM_UD_RATIO),
        get_float("asym_s_target", DEFAULT_ASYM_S_TARGET),
        get_float("xi", DEFAULT_XI),
        get_int("p_num", DEFAULT_P_NUM),
        get_int("t_num", DEFAULT_T_NUM),
        get(raw, "output_dir", nothing),
        get_bool("run_reference", default_run_reference),
        Symbol(get(raw, "model_kind", "PNJL")),
        get_float("ds_rho", DEFAULT_DS_RHO),
        get_float("dsmax_rho", DEFAULT_DSMAX_RHO),
        get_int("max_steps", DEFAULT_MAX_STEPS),
        get_float("branch_jump_tol", U.DEFAULT_BRANCH_JUMP_TOL),
        get_bool("trho_reverse_rho", false),
        get_float("pressure_gap_tol", DEFAULT_PRESSURE_GAP_TOL),
    )
end

function _ensure_models_loaded(repo_root::String)
    if !isdefined(Main, :Models)
        Base.include(Main, joinpath(repo_root, "src", "models", "Models.jl"))
    end
    return Main.Models
end

function _effective_config(repo_root::String, cli::CliConfig; run_id::String, output_dir::String)
    cli.rho_min < cli.rho_max || throw(ArgumentError("rho-min must be smaller than rho-max"))
    U._require_positive("rho-step", cli.rho_step)
    U._require_positive("p-num", cli.p_num)
    U._require_positive("t-num", cli.t_num)
    U._require_positive("ds-rho", cli.ds_rho)
    U._require_positive("dsmax-rho", cli.dsmax_rho)
    U._require_positive("max-steps", cli.max_steps)
    U._require_positive("pressure-gap-tol", cli.pressure_gap_tol)

    width = cli.rho_max - cli.rho_min
    rho_start = cli.rho_min + min(0.05, max(0.02, 0.05 * width))
    rho_start = min(max(rho_start, cli.rho_min), cli.rho_max)

    return EffectiveConfig(
        repo_root=repo_root,
        output_dir=output_dir,
        run_id=run_id,
        T_MeV=cli.T_MeV,
        rho_min=cli.rho_min,
        rho_max=cli.rho_max,
        rho_step=cli.rho_step,
        rho_start=rho_start,
        asym_ud_ratio_target=cli.asym_ud_ratio_target,
        asym_s_target=cli.asym_s_target,
        xi=cli.xi,
        p_num=cli.p_num,
        t_num=cli.t_num,
        run_reference=cli.run_reference,
        model_kind=cli.model_kind,
        ds_rho=cli.ds_rho,
        dsmax_rho=cli.dsmax_rho,
        max_steps=cli.max_steps,
        branch_jump_tol=cli.branch_jump_tol,
        trho_reverse_rho=cli.trho_reverse_rho,
        pressure_gap_tol=cli.pressure_gap_tol,
    )
end

function _mode(cfg::EffectiveConfig, rho_target::Real)
    return Main.Models.FixedAsymmetricRho(
        Float64(rho_target),
        cfg.asym_ud_ratio_target,
        cfg.asym_s_target,
    )
end

function _residual_fn(model, cfg::EffectiveConfig)
    T_fm = cfg.T_MeV / HBARC_MEV_FM
    thermal_nodes = Main.Models.cached_nodes(
        cfg.p_num,
        cfg.t_num;
        p_max_inv_fm=Main.Models.thermal_p_max_inv_fm(model),
    )
    params = Main.Models.GapParams(
        T_fm,
        thermal_nodes,
        cfg.xi;
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        model_kind=cfg.model_kind,
    )
    return function (x, p)
        mode = Main.Models.FixedAsymmetricRho(
            Float64(p.rho_target),
            Float64(p.asym_ud_ratio_target),
            Float64(p.asym_s_target),
        )
        F = Vector{Float64}(undef, 8)
        Main.Models.explicit_residual!(F, Float64.(x), Float64[p.T_fm], params, mode)
        return F
    end
end

function _safe_diagnostics(model, cfg::EffectiveConfig, x, rho_target::Real)
    T_fm = cfg.T_MeV / HBARC_MEV_FM
    x_state = Float64.(x[1:5])
    mu_vec = Float64.(x[6:8])
    residual_norm = try
        residual! = Main.Models.build_residual!(
            model,
            _mode(cfg, rho_target),
            T_fm;
            xi=cfg.xi,
            p_num=cfg.p_num,
            t_num=cfg.t_num,
        )
        F = zeros(Float64, 8)
        residual!(F, Float64.(x))
        norm(F)
    catch
        NaN
    end
    pressure = try
        Float64(Main.Models.model_pressure(model, x_state, mu_vec, T_fm; xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num))
    catch
        NaN
    end
    rho_vec = try
        Float64.(Main.Models.model_rho(model, x_state, mu_vec, T_fm; xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num))
    catch
        fill(NaN, 3)
    end
    masses = try
        st = Main.Models.meanfield_state(x_state)
        Float64.(Main.Models.calculate_mass_vec(model, st.phi)) .* HBARC_MEV_FM
    catch
        fill(NaN, 3)
    end
    rho_norm = all(isfinite, rho_vec) ? sum(rho_vec) / (3.0 * Main.Models.ρ0) : NaN
    ud_ratio = isfinite(rho_vec[1]) && isfinite(rho_vec[2]) && abs(rho_vec[2]) > 1e-12 ? rho_vec[1] / rho_vec[2] : NaN
    return residual_norm, pressure, rho_vec, rho_norm, ud_ratio, masses
end

function _run_palc(Models, cfg::EffectiveConfig)
    model = Models.create_model(cfg.model_kind)
    T_fm = cfg.T_MeV / HBARC_MEV_FM
    initial_elapsed = @elapsed initial = Models.solve(
        model,
        _mode(cfg, cfg.rho_start),
        T_fm;
        seed_strategy=Models.MultiSeed(),
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        mu0=[1.0, 1.0, 1.0],
        residual_norm_max=1e-3,
        iterations=300,
    )
    x0 = Float64.(initial.solution[1:8])

    residual_fn = _residual_fn(model, cfg)
    jacobian_fn = (x, p) -> U._finite_difference_jacobian(residual_fn, x, p)

    record_fn = function (x, rho_target; kwargs...)
        residual_norm, pressure, rho_vec, rho_norm, ud_ratio, masses =
            _safe_diagnostics(model, cfg, x, Float64(rho_target))
        return (
            phi_u=Float64(x[1]),
            phi_d=Float64(x[2]),
            phi_s=Float64(x[3]),
            Phi=Float64(x[4]),
            PhiBar=Float64(x[5]),
            mu_u_MeV=Float64(x[6]) * HBARC_MEV_FM,
            mu_d_MeV=Float64(x[7]) * HBARC_MEV_FM,
            mu_s_MeV=Float64(x[8]) * HBARC_MEV_FM,
            rho_target=Float64(rho_target),
            rho_norm=rho_norm,
            rho_u=rho_vec[1],
            rho_d=rho_vec[2],
            rho_s=rho_vec[3],
            ud_ratio=ud_ratio,
            residual_norm=residual_norm,
            pressure_fm4=pressure,
            M_u_MeV=masses[1],
            M_d_MeV=masses[2],
            M_s_MeV=masses[3],
        )
    end

    params = (
        T_fm=T_fm,
        rho_target=cfg.rho_start,
        asym_ud_ratio_target=cfg.asym_ud_ratio_target,
        asym_s_target=cfg.asym_s_target,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
    )
    prob = BifurcationKit.BifurcationProblem(
        residual_fn,
        x0,
        params,
        (@optic _.rho_target);
        J=jacobian_fn,
        record_from_solution=record_fn,
    )
    opts = BifurcationKit.ContinuationPar(
        p_min=cfg.rho_min,
        p_max=cfg.rho_max,
        ds=cfg.ds_rho,
        dsmax=cfg.dsmax_rho,
        max_steps=cfg.max_steps,
        detect_bifurcation=0,
        newton_options=BifurcationKit.NewtonPar(
            tol=PALC_NEWTON_TOL,
            max_iterations=PALC_NEWTON_MAX_ITERATIONS,
            verbose=false,
        ),
    )

    branch = Ref{Any}(nothing)
    continuation_elapsed = @elapsed begin
        branch[] = BifurcationKit.continuation(prob, BifurcationKit.PALC(), opts; bothside=true, verbosity=0)
    end
    return (
        result=branch[],
        initial_solution=initial,
        initial_solve_wall_time_s=initial_elapsed,
        continuation_wall_time_s=continuation_elapsed,
        wall_time_s=initial_elapsed + continuation_elapsed,
    )
end

function _branch_rows(br)
    rows = Vector{Vector}()
    for row in br.branch
        push!(rows, Any[
            getproperty(row, :step),
            getproperty(row, :param),
            getproperty(row, :rho_target),
            getproperty(row, :rho_norm),
            getproperty(row, :phi_u),
            getproperty(row, :phi_d),
            getproperty(row, :phi_s),
            getproperty(row, :Phi),
            getproperty(row, :PhiBar),
            getproperty(row, :mu_u_MeV),
            getproperty(row, :mu_d_MeV),
            getproperty(row, :mu_s_MeV),
            getproperty(row, :rho_u),
            getproperty(row, :rho_d),
            getproperty(row, :rho_s),
            getproperty(row, :ud_ratio),
            getproperty(row, :residual_norm),
            getproperty(row, :pressure_fm4),
            getproperty(row, :M_u_MeV),
            getproperty(row, :M_d_MeV),
            getproperty(row, :M_s_MeV),
            getproperty(row, :itnewton),
            getproperty(row, :itlinear),
            getproperty(row, :ds),
        ])
    end
    return rows
end

function _write_branch_csv(path::String, br)
    header = [
        "step", "param_rho_target", "rho_target", "rho_norm",
        "phi_u", "phi_d", "phi_s", "Phi", "PhiBar",
        "mu_u_MeV", "mu_d_MeV", "mu_s_MeV",
        "rho_u", "rho_d", "rho_s", "ud_ratio",
        "residual_norm", "pressure_fm4",
        "M_u_MeV", "M_d_MeV", "M_s_MeV",
        "itnewton", "itlinear", "ds",
    ]
    return U._write_csv(path, header, _branch_rows(br))
end

function _special_points(br)
    points = NamedTuple[]
    for sp in br.specialpoint
        param = try Float64(getproperty(sp, :param)) catch; NaN end
        push!(points, (
            source=:bifurcationkit,
            type=try getproperty(sp, :type) catch; :unknown end,
            idx=try Int(getproperty(sp, :idx)) catch; -1 end,
            step=try Int(getproperty(sp, :step)) catch; -1 end,
            rho_target=param,
            status=try getproperty(sp, :status) catch; :unknown end,
        ))
    end
    return points
end

function _local_param_extrema(br)
    branch = br.branch
    n = length(branch)
    points = NamedTuple[]
    n < 3 && return points
    vals = Float64[getproperty(row, :rho_target) for row in branch]
    for i in 2:(n - 1)
        d1 = vals[i] - vals[i - 1]
        d2 = vals[i + 1] - vals[i]
        if isfinite(d1) && isfinite(d2) && d1 * d2 < 0
            push!(points, (
                source=:local_extremum,
                type=d1 > 0 ? :param_max : :param_min,
                idx=i,
                step=try Int(getproperty(branch[i], :step)) catch; i end,
                rho_target=vals[i],
                status=:detected,
            ))
        end
    end
    return points
end

function _state_jump_metrics(rows::Vector{<:Any}; tol::Float64, cols=5:12)
    length(rows) < 2 && return (count=0, max_jump=0.0)
    max_jump = 0.0
    count = 0
    for i in 2:length(rows)
        a = Float64.(rows[i - 1][cols])
        b = Float64.(rows[i][cols])
        if length(a) == 8
            @views a[6:8] ./= HBARC_MEV_FM
            @views b[6:8] ./= HBARC_MEV_FM
        end
        jump = norm(b .- a)
        isfinite(jump) || continue
        max_jump = max(max_jump, jump)
        jump > tol && (count += 1)
    end
    return (count=count, max_jump=max_jump)
end

function _config_dict(cfg::EffectiveConfig)
    return Dict(
        "repo_root" => cfg.repo_root,
        "output_dir" => cfg.output_dir,
        "run_id" => cfg.run_id,
        "T_MeV" => cfg.T_MeV,
        "rho_min" => cfg.rho_min,
        "rho_max" => cfg.rho_max,
        "rho_step" => cfg.rho_step,
        "rho_start" => cfg.rho_start,
        "asym_ud_ratio_target" => cfg.asym_ud_ratio_target,
        "asym_s_target" => cfg.asym_s_target,
        "xi" => cfg.xi,
        "p_num" => cfg.p_num,
        "t_num" => cfg.t_num,
        "run_reference" => cfg.run_reference,
        "model_kind" => String(cfg.model_kind),
        "ds_rho" => cfg.ds_rho,
        "dsmax_rho" => cfg.dsmax_rho,
        "max_steps" => cfg.max_steps,
        "branch_jump_tol" => cfg.branch_jump_tol,
        "trho_reverse_rho" => cfg.trho_reverse_rho,
        "pressure_gap_tol" => cfg.pressure_gap_tol,
    )
end

function _palc_metrics(palc, cfg::EffectiveConfig)
    rows = _branch_rows(palc.result)
    finite_residuals = [Float64(row[17]) for row in rows if isfinite(Float64(row[17]))]
    jumps = _state_jump_metrics(rows; tol=cfg.branch_jump_tol, cols=5:12)
    rho_values = [Float64(row[3]) for row in rows if isfinite(Float64(row[3]))]
    special = vcat(_special_points(palc.result), _local_param_extrema(palc.result))
    fold_like = filter(p -> p.type in (:fold, :param_max, :param_min), special)
    return (
        status=:ok,
        branch_points=length(rows),
        finite_residual_count=length(finite_residuals),
        residual_norm_max=isempty(finite_residuals) ? NaN : maximum(finite_residuals),
        residual_norm_min=isempty(finite_residuals) ? NaN : minimum(finite_residuals),
        rho_min_seen=isempty(rho_values) ? NaN : minimum(rho_values),
        rho_max_seen=isempty(rho_values) ? NaN : maximum(rho_values),
        branch_jump_count=jumps.count,
        max_state_jump=jumps.max_jump,
        newton_iterations_total=sum(Int(row[22]) for row in rows),
        linear_iterations_total=sum(Int(row[23]) for row in rows),
        special_points=special,
        fold_like_count=length(fold_like),
        initial_solve_wall_time_s=palc.initial_solve_wall_time_s,
        continuation_wall_time_s=palc.continuation_wall_time_s,
        wall_time_s=palc.wall_time_s,
    )
end

function _read_trho_csv(path::String)
    rows = NamedTuple[]
    isfile(path) || return rows
    lines = readlines(path)
    isempty(lines) && return rows
    header = split(lines[1], ",")
    idx = Dict(strip(name) => i for (i, name) in pairs(header))
    for line in lines[2:end]
        isempty(strip(line)) && continue
        cols = split(line, ",")
        getcol(name, default="NaN") = get(cols, get(idx, name, 0), default)
        push!(rows, (
            T_MeV=parse(Float64, getcol("T_MeV")),
            rho=parse(Float64, getcol("rho")),
            xi=parse(Float64, getcol("xi")),
            mu_u_MeV=parse(Float64, getcol("mu_u_MeV")),
            mu_d_MeV=parse(Float64, getcol("mu_d_MeV")),
            mu_s_MeV=parse(Float64, getcol("mu_s_MeV")),
            phi_u=parse(Float64, getcol("phi_u")),
            phi_d=parse(Float64, getcol("phi_d")),
            phi_s=parse(Float64, getcol("phi_s")),
            Phi=parse(Float64, getcol("Phi1")),
            PhiBar=parse(Float64, getcol("Phi2")),
            pressure_fm4=parse(Float64, getcol("pressure_fm4")),
            residual_norm=parse(Float64, getcol("residual_norm")),
            iterations=parse(Int, getcol("iterations", "-1")),
            converged=lowercase(strip(getcol("converged", "false"))) == "true",
        ))
    end
    return rows
end

function _run_trho_reference(Models, cfg::EffectiveConfig)
    path = joinpath(cfg.output_dir, "trho_scan_reference.csv")
    rho_values = collect(cfg.rho_min:cfg.rho_step:cfg.rho_max)
    if isempty(rho_values) || last(rho_values) < cfg.rho_max
        push!(rho_values, cfg.rho_max)
    end
    elapsed = @elapsed stats = Models.run_trho_scan(
        T_values=[cfg.T_MeV],
        rho_values=rho_values,
        xi_values=[cfg.xi],
        output_path=path,
        overwrite=true,
        resume=false,
        reverse_rho=cfg.trho_reverse_rho,
        seed_policy=:hybrid_continuity,
        constraint_mode=:fixed_asymmetric_rho,
        asym_ud_ratio_target=cfg.asym_ud_ratio_target,
        asym_s_target=cfg.asym_s_target,
        solver_backend=:models,
        model_kind=cfg.model_kind,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        iterations=300,
        residual_norm_max=1e-3,
    )
    rows = _read_trho_csv(path)
    jumps = _state_jump_metrics([Any[
            row.T_MeV, row.rho, row.xi, row.rho,
            row.phi_u, row.phi_d, row.phi_s, row.Phi, row.PhiBar,
            row.mu_u_MeV, row.mu_d_MeV, row.mu_s_MeV,
        ] for row in rows]; tol=cfg.branch_jump_tol, cols=5:12)
    residuals = [row.residual_norm for row in rows if isfinite(row.residual_norm)]
    return (
        path=path,
        wall_time_s=elapsed,
        stats=stats,
        points=length(rows),
        failure_count=count(row -> !row.converged, rows),
        branch_jump_count=jumps.count,
        max_state_jump=jumps.max_jump,
        finite_residual_count=length(residuals),
        residual_norm_min=isempty(residuals) ? NaN : minimum(residuals),
        residual_norm_max=isempty(residuals) ? NaN : maximum(residuals),
        iterations_total=sum(max(row.iterations, 0) for row in rows),
    )
end

function _skipped_trho_reference()
    return (
        path=nothing,
        wall_time_s=0.0,
        stats=nothing,
        points=0,
        failure_count=0,
        branch_jump_count=0,
        max_state_jump=0.0,
        finite_residual_count=0,
        residual_norm_min=NaN,
        residual_norm_max=NaN,
        iterations_total=0,
    )
end

function _nearest_palc_row(palc_rows, rho::Float64)
    isempty(palc_rows) && return nothing
    best = nothing
    best_delta = Inf
    for row in palc_rows
        palc_rho = Float64(row[3])
        delta = abs(palc_rho - rho)
        if isfinite(delta) && delta < best_delta
            best = row
            best_delta = delta
        end
    end
    return best === nothing ? nothing : (row=best, rho_delta=best_delta)
end

function _branch_pressure_comparison(palc, trho_metrics, cfg::EffectiveConfig)
    palc_rows = _branch_rows(palc.result)
    trho_rows = trho_metrics.path === nothing ? NamedTuple[] : _read_trho_csv(trho_metrics.path)
    max_match_delta = max(cfg.rho_step, cfg.dsmax_rho)
    lower_pressure_examples = NamedTuple[]
    higher_pressure_examples = NamedTuple[]
    matched = 0
    lower_count = 0
    higher_count = 0
    max_lower_gap = 0.0
    max_higher_gap = 0.0

    for row in trho_rows
        isfinite(row.pressure_fm4) || continue
        nearest = _nearest_palc_row(palc_rows, Float64(row.rho))
        nearest === nothing && continue
        nearest.rho_delta <= max_match_delta || continue
        palc_row = nearest.row
        palc_pressure = Float64(palc_row[18])
        isfinite(palc_pressure) || continue

        matched += 1
        gap = palc_pressure - row.pressure_fm4
        mu_distance = norm((
            Float64(palc_row[10]) - row.mu_u_MeV,
            Float64(palc_row[11]) - row.mu_d_MeV,
            Float64(palc_row[12]) - row.mu_s_MeV,
        ))
        example = (
            rho=Float64(row.rho),
            nearest_palc_rho=Float64(palc_row[3]),
            rho_delta=Float64(nearest.rho_delta),
            trho_pressure_fm4=Float64(row.pressure_fm4),
            palc_pressure_fm4=palc_pressure,
            palc_minus_trho_pressure_fm4=gap,
            mu_distance_MeV=mu_distance,
            trho_mu_u_MeV=Float64(row.mu_u_MeV),
            trho_mu_d_MeV=Float64(row.mu_d_MeV),
            trho_mu_s_MeV=Float64(row.mu_s_MeV),
            palc_mu_u_MeV=Float64(palc_row[10]),
            palc_mu_d_MeV=Float64(palc_row[11]),
            palc_mu_s_MeV=Float64(palc_row[12]),
        )
        if gap > cfg.pressure_gap_tol
            lower_count += 1
            max_lower_gap = max(max_lower_gap, gap)
            length(lower_pressure_examples) < 5 && push!(lower_pressure_examples, example)
        elseif -gap > cfg.pressure_gap_tol
            higher_count += 1
            max_higher_gap = max(max_higher_gap, -gap)
            length(higher_pressure_examples) < 5 && push!(higher_pressure_examples, example)
        end
    end

    return (
        matched_points=matched,
        max_match_delta=max_match_delta,
        reference_lower_pressure_count=lower_count,
        reference_higher_pressure_count=higher_count,
        max_reference_lower_pressure_gap=max_lower_gap,
        max_reference_higher_pressure_gap=max_higher_gap,
        lower_pressure_examples=lower_pressure_examples,
        higher_pressure_examples=higher_pressure_examples,
    )
end

function _decision(palc_metrics, trho_metrics, branch_comparison)
    has_branch = palc_metrics.branch_points > 0 && palc_metrics.finite_residual_count > 0
    lower_pressure_reference = branch_comparison.reference_lower_pressure_count > 0
    better_continuity = has_branch && palc_metrics.branch_jump_count < trho_metrics.branch_jump_count
    recommended = has_branch && (lower_pressure_reference || better_continuity || palc_metrics.fold_like_count > 0)
    reason = if !has_branch
        "PALC produced no finite FixedAsymmetricRho branch rows"
    elseif lower_pressure_reference
        "reference continuation selected lower-pressure roots than the PALC pressure branch at matched rho values"
    elseif better_continuity
        "PALC produced a finite branch with fewer branch jumps than run_trho_scan"
    elseif palc_metrics.fold_like_count > 0
        "PALC produced a finite branch and fold-like diagnostics"
    else
        "PALC produced finite rows but no continuity or pressure-branch advantage over run_trho_scan in this window"
    end
    return (recommended=recommended, reason=reason)
end

function _write_report(path::String, cfg::EffectiveConfig, palc_metrics, trho_metrics, branch_comparison, decision)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# FixedAsymmetricRho PNJL PALC Spike Report")
        println(io)
        println(io, "## Configuration")
        println(io, "- T_MeV: $(cfg.T_MeV)")
        println(io, "- rho window: $(cfg.rho_min) to $(cfg.rho_max)")
        println(io, "- asym_ud_ratio_target: $(cfg.asym_ud_ratio_target)")
        println(io, "- asym_s_target: $(cfg.asym_s_target)")
        println(io, "- xi: $(cfg.xi)")
        println(io, "- p_num/t_num: $(cfg.p_num)/$(cfg.t_num)")
        println(io, "- run_reference: $(cfg.run_reference)")
        println(io, "- trho_reverse_rho: $(cfg.trho_reverse_rho)")
        println(io)
        println(io, "## PALC Branch")
        println(io, "- branch_points: $(palc_metrics.branch_points)")
        println(io, "- finite_residual_count: $(palc_metrics.finite_residual_count)")
        println(io, "- residual_norm_max: $(U._fmt(palc_metrics.residual_norm_max))")
        println(io, "- rho_seen: $(U._fmt(palc_metrics.rho_min_seen)) to $(U._fmt(palc_metrics.rho_max_seen))")
        println(io, "- branch_jump_count: $(palc_metrics.branch_jump_count)")
        println(io, "- max_state_jump: $(U._fmt(palc_metrics.max_state_jump))")
        println(io, "- newton_iterations_total: $(palc_metrics.newton_iterations_total)")
        println(io, "- wall_time_s: $(U._fmt(palc_metrics.wall_time_s))")
        println(io, "- fold_like_count: $(palc_metrics.fold_like_count)")
        println(io)
        println(io, "## Trho Reference")
        println(io, "- points: $(trho_metrics.points)")
        println(io, "- failure_count: $(trho_metrics.failure_count)")
        println(io, "- branch_jump_count: $(trho_metrics.branch_jump_count)")
        println(io, "- max_state_jump: $(U._fmt(trho_metrics.max_state_jump))")
        println(io, "- finite_residual_count: $(trho_metrics.finite_residual_count)")
        println(io, "- residual_norm_min/max: $(U._fmt(trho_metrics.residual_norm_min)) / $(U._fmt(trho_metrics.residual_norm_max))")
        println(io, "- iterations_total: $(trho_metrics.iterations_total)")
        println(io, "- wall_time_s: $(U._fmt(trho_metrics.wall_time_s))")
        println(io)
        println(io, "## Branch Pressure Comparison")
        println(io, "- matched_points: $(branch_comparison.matched_points)")
        println(io, "- reference_lower_pressure_count: $(branch_comparison.reference_lower_pressure_count)")
        println(io, "- max_reference_lower_pressure_gap: $(U._fmt(branch_comparison.max_reference_lower_pressure_gap))")
        if !isempty(branch_comparison.lower_pressure_examples)
            println(io, "- first_lower_pressure_example: $(branch_comparison.lower_pressure_examples[1])")
        end
        println(io)
        println(io, "## Decision")
        println(io, "- continue_recommended: $(decision.recommended)")
        println(io, "- reason: $(decision.reason)")
    end
    return path
end

function _write_outputs(cfg::EffectiveConfig, palc, trho_metrics)
    mkpath(cfg.output_dir)
    branch_path = joinpath(cfg.output_dir, "palc_fixedasymrho_branch.csv")
    summary_path = joinpath(cfg.output_dir, "comparison_summary.json")
    report_path = joinpath(cfg.output_dir, "comparison_report.md")

    _write_branch_csv(branch_path, palc.result)
    palc_metrics = _palc_metrics(palc, cfg)
    branch_comparison = _branch_pressure_comparison(palc, trho_metrics, cfg)
    decision = _decision(palc_metrics, trho_metrics, branch_comparison)
    summary = (
        config=_config_dict(cfg),
        artifacts=(
            palc_branch=branch_path,
            trho_scan_reference=trho_metrics.path,
            report=report_path,
        ),
        palc=palc_metrics,
        trho_reference=trho_metrics,
        branch_comparison=branch_comparison,
        decision=decision,
    )
    U._write_json(summary_path, summary)
    _write_report(report_path, cfg, palc_metrics, trho_metrics, branch_comparison, decision)
    return summary
end

function run(args::Vector{String}; repo_root::String, default_run_reference::Bool)
    cli = _parse_cli(args; default_run_reference=default_run_reference)
    Models = _ensure_models_loaded(repo_root)
    run_id = U._run_id()
    output_dir = cli.output_dir === nothing ?
        joinpath(repo_root, DEFAULT_OUTPUT_REL, run_id) :
        (isabspath(cli.output_dir) ? cli.output_dir : joinpath(repo_root, cli.output_dir))
    cfg = _effective_config(repo_root, cli; run_id=run_id, output_dir=output_dir)

    palc = _run_palc(Models, cfg)
    trho_metrics = cfg.run_reference ? _run_trho_reference(Models, cfg) : _skipped_trho_reference()
    summary = _write_outputs(cfg, palc, trho_metrics)

    println("FixedAsymmetricRho PALC spike output: $(cfg.output_dir)")
    println("continue_recommended: $(summary.decision.recommended)")
    println("reason: $(summary.decision.reason)")
    return summary
end

function main_run(args::Vector{String}; repo_root::String, default_run_reference::Bool)
    try
        run(args; repo_root=repo_root, default_run_reference=default_run_reference)
    catch err
        println(stderr, "FixedAsymmetricRho PALC spike failed.")
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        exit(1)
    end
    return nothing
end

end # module
