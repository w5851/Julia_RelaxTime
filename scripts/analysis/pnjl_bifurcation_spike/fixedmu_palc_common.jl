module FixedMuPALCSpike

using Dates
using LinearAlgebra
using Printf

try
    @eval import BifurcationKit
    @eval import BifurcationKit: @optic
catch err
    error("BifurcationKit is not available. Run with --project=scripts/analysis/pnjl_bifurcation_spike and instantiate that environment. Original error: $(err)")
end

try
    @eval import JSON
catch err
    error("JSON is not available from the stacked repository environment. Run the script from the repository root or keep the root project on LOAD_PATH. Original error: $(err)")
end

const HBARC_MEV_FM = 197.327
const DEFAULT_OUTPUT_REL = joinpath("data", "outputs", "results", "analysis", "palc_fixedmu_spike")
const DEFAULT_T_MEV = 120.0
const DEFAULT_MU_MIN_MEV = 240.0
const DEFAULT_MU_MAX_MEV = 380.0
const DEFAULT_MU_STEP_MEV = 5.0
const DEFAULT_XI = 0.0
const DEFAULT_P_NUM = 8
const DEFAULT_T_NUM = 4
const DEFAULT_DS_MEV = 2.0
const DEFAULT_DSMAX_MEV = 8.0
const DEFAULT_MAX_STEPS = 160
const DEFAULT_BRANCH_JUMP_TOL = 0.5
const PALC_NEWTON_TOL = 1e-7
const PALC_NEWTON_MAX_ITERATIONS = 10
const SPINODAL_HINT_PAD_MEV = 20.0
const FOLD_HINT_TARGET_MEV = 10.0
const PHASE_REFERENCE_T_GRID = collect(120.0:5.0:150.0)
const PHASE_REFERENCE_RHO_GRID = collect(0.1:0.1:3.0)
const PHASE_REFERENCE_ITERATIONS = 40

struct CliConfig
    T_MeV::Union{Nothing, Float64}
    mu_min_MeV::Union{Nothing, Float64}
    mu_max_MeV::Union{Nothing, Float64}
    mu_step_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    output_dir::Union{Nothing, String}
    run_reference::Bool
    model_kind::Symbol
    ds_MeV::Float64
    dsmax_MeV::Float64
    max_steps::Int
    branch_jump_tol::Float64
end

Base.@kwdef struct EffectiveConfig
    repo_root::String
    output_dir::String
    run_id::String
    T_MeV::Float64
    mu_min_MeV::Float64
    mu_max_MeV::Float64
    mu_step_MeV::Float64
    mu_start_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    run_reference::Bool
    model_kind::Symbol
    ds_MeV::Float64
    dsmax_MeV::Float64
    max_steps::Int
    branch_jump_tol::Float64
end

function _usage()
    return """
    Usage:
      julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedmu_palc.jl [options]

    Options:
      --T-MeV=<value>             Fixed temperature in MeV. Defaults to first phase spinodal T, then 120.
      --mu-min-MeV=<value>        Lower quark chemical potential bound in MeV.
      --mu-max-MeV=<value>        Upper quark chemical potential bound in MeV.
      --mu-step-MeV=<value>       Reference T-mu scan step in MeV. Default: 5.
      --xi=<value>                Anisotropy parameter. Default: 0.
      --p-num=<int>               Momentum nodes. Default: 8.
      --t-num=<int>               Theta nodes. Default: 4.
      --output-dir=<path>         Output directory. Default: data/outputs/results/analysis/palc_fixedmu_spike/<run_id>.
      --run-reference=<bool>      Run phase and run_tmu_scan references.
      --ds-MeV=<value>            PALC initial pseudo-arclength scale, reported in MeV units. Default: 2.
      --dsmax-MeV=<value>         PALC max pseudo-arclength scale, reported in MeV units. Default: 8.
      --max-steps=<int>           PALC max continuation steps. Default: 160.
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
        key = replace(lowercase(parts[1]), "-" => "_")
        raw[key] = parts[2]
    end

    get_float(key, default=nothing) = haskey(raw, key) ? parse(Float64, raw[key]) : default
    get_int(key, default) = haskey(raw, key) ? parse(Int, raw[key]) : default
    function get_bool(key, default)
        haskey(raw, key) || return default
        value = lowercase(strip(raw[key]))
        value in ("1", "true", "yes", "y") && return true
        value in ("0", "false", "no", "n") && return false
        throw(ArgumentError("$(key) must be true or false, got $(raw[key])"))
    end

    return CliConfig(
        get_float("t_mev", nothing),
        get_float("mu_min_mev", nothing),
        get_float("mu_max_mev", nothing),
        get_float("mu_step_mev", DEFAULT_MU_STEP_MEV),
        get_float("xi", DEFAULT_XI),
        get_int("p_num", DEFAULT_P_NUM),
        get_int("t_num", DEFAULT_T_NUM),
        get(raw, "output_dir", nothing),
        get_bool("run_reference", default_run_reference),
        Symbol(get(raw, "model_kind", "PNJL")),
        get_float("ds_mev", DEFAULT_DS_MEV),
        get_float("dsmax_mev", DEFAULT_DSMAX_MEV),
        get_int("max_steps", DEFAULT_MAX_STEPS),
        get_float("branch_jump_tol", DEFAULT_BRANCH_JUMP_TOL),
    )
end

function _require_positive(name::String, value::Real)
    value > 0 || throw(ArgumentError("$(name) must be positive, got $(value)"))
    return value
end

function _ensure_models_loaded(repo_root::String)
    if !isdefined(Main, :Models)
        Base.include(Main, joinpath(repo_root, "src", "models", "Models.jl"))
    end
    return Main.Models
end

function _run_id()
    return Dates.format(Dates.now(), dateformat"yyyymmdd_HHMMSS")
end

function _fmt(v)
    if v isa Real
        fv = Float64(v)
        isfinite(fv) || return "NaN"
        return @sprintf("%.12g", fv)
    end
    return string(v)
end

function _csv_quote(s)
    text = string(s)
    if occursin(",", text) || occursin("\"", text) || occursin("\n", text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function _write_csv(path::String, header::Vector{String}, rows::Vector{Vector})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join((_csv_quote(_fmt(v)) for v in row), ","))
        end
    end
    return path
end

function _json_safe(x)
    x === nothing && return nothing
    x isa Symbol && return String(x)
    if x isa Bool
        return x
    elseif x isa Integer
        return Int(x)
    elseif x isa Real
        y = Float64(x)
        return isfinite(y) ? y : nothing
    elseif x isa AbstractString
        return x
    elseif x isa AbstractDict
        return Dict(string(k) => _json_safe(v) for (k, v) in x)
    elseif x isa NamedTuple
        return Dict(string(k) => _json_safe(v) for (k, v) in pairs(x))
    elseif x isa AbstractVector || x isa Tuple
        return [_json_safe(v) for v in x]
    end
    return string(x)
end

function _json_safe(cfg::EffectiveConfig)
    return Dict(
        "repo_root" => cfg.repo_root,
        "output_dir" => cfg.output_dir,
        "run_id" => cfg.run_id,
        "T_MeV" => cfg.T_MeV,
        "mu_min_MeV" => cfg.mu_min_MeV,
        "mu_max_MeV" => cfg.mu_max_MeV,
        "mu_step_MeV" => cfg.mu_step_MeV,
        "mu_start_MeV" => cfg.mu_start_MeV,
        "xi" => cfg.xi,
        "p_num" => cfg.p_num,
        "t_num" => cfg.t_num,
        "run_reference" => cfg.run_reference,
        "model_kind" => String(cfg.model_kind),
        "ds_MeV" => cfg.ds_MeV,
        "dsmax_MeV" => cfg.dsmax_MeV,
        "max_steps" => cfg.max_steps,
        "branch_jump_tol" => cfg.branch_jump_tol,
    )
end

function _write_json(path::String, data)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, _json_safe(data), 4)
        println(io)
    end
    return path
end

function _finite_or_nan(value)
    try
        v = Float64(value)
        return isfinite(v) ? v : NaN
    catch
        return NaN
    end
end

function _safe_model_diagnostics(model, x, T_fm, mu_fm; xi::Float64, p_num::Int, t_num::Int)
    mu_vec = Main.Models.normalize_mu_vec(mu_fm)
    residual = try
        Vector(Main.Models.gap_residual(model, x, T_fm, mu_vec; xi=xi, p_num=p_num, t_num=t_num))
    catch
        fill(NaN, 5)
    end
    residual_norm = all(isfinite, residual) ? norm(Float64.(residual)) : NaN
    pressure = try
        Float64(Main.Models.model_pressure(model, x, mu_vec, T_fm; xi=xi, p_num=p_num, t_num=t_num))
    catch
        NaN
    end
    masses = try
        st = Main.Models.meanfield_state(x)
        Float64.(Main.Models.calculate_mass_vec(model, st.phi)) .* HBARC_MEV_FM
    catch
        fill(NaN, 3)
    end
    return residual_norm, pressure, masses
end

function _finite_difference_jacobian(residual_fn, x, p; rel_step::Float64=sqrt(eps(Float64)))
    x0 = Float64.(x)
    f0 = Float64.(residual_fn(x0, p))
    n_out = length(f0)
    n_in = length(x0)
    J = Matrix{Float64}(undef, n_out, n_in)
    xp = copy(x0)
    xm = copy(x0)
    for j in 1:n_in
        h = rel_step * max(1.0, abs(x0[j]))
        xp[j] = x0[j] + h
        xm[j] = x0[j] - h
        fp = Float64.(residual_fn(xp, p))
        fm = Float64.(residual_fn(xm, p))
        @views J[:, j] .= (fp .- fm) ./ (2h)
        xp[j] = x0[j]
        xm[j] = x0[j]
    end
    return J
end

function _run_phase_reference(Models, cfg::CliConfig, output_dir::String)
    phase_dir = joinpath(output_dir, "phase_reference")
    mkpath(phase_dir)
    elapsed = @elapsed result = Models.run_phase_pipeline(
        cfg.model_kind;
        mode=:research,
        T_grid=PHASE_REFERENCE_T_GRID,
        rho_grid=PHASE_REFERENCE_RHO_GRID,
        xi=cfg.xi,
        output_dir=phase_dir,
        profile=:palc_spike,
        solver_backend=:models,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        iterations=PHASE_REFERENCE_ITERATIONS,
        compute_crossover=false,
        cep_strategy=:interpolate,
        promote_reference=false,
    )

    selected = nothing
    spinodal_rows = sort(collect(result.spinodal); by=row -> row.T_MeV)
    for row in spinodal_rows
        mu_h = _finite_or_nan(row.mu_spinodal_hadron_MeV)
        mu_q = _finite_or_nan(row.mu_spinodal_quark_MeV)
        if isfinite(mu_h) && isfinite(mu_q)
            selected = (
                T_MeV=Float64(row.T_MeV),
                mu_spinodal_hadron_MeV=mu_h,
                mu_spinodal_quark_MeV=mu_q,
                mu_min_MeV=min(mu_h, mu_q) - SPINODAL_HINT_PAD_MEV,
                mu_max_MeV=max(mu_h, mu_q) + SPINODAL_HINT_PAD_MEV,
            )
            break
        end
    end

    return (
        ran=true,
        phase_output_dir=phase_dir,
        wall_time_s=elapsed,
        spinodal_count=length(result.spinodal),
        boundary_count=length(result.first_order_boundary),
        cep_found=result.cep.found,
        selected_from_phase=selected !== nothing,
        selected=selected,
        fallback_reason=selected === nothing ? "phase pipeline produced no finite two-sided spinodal row" : nothing,
    )
end

function _skipped_phase_reference(reason::String)
    return (
        ran=false,
        phase_output_dir=nothing,
        wall_time_s=0.0,
        spinodal_count=0,
        boundary_count=0,
        cep_found=false,
        selected_from_phase=false,
        selected=nothing,
        fallback_reason=reason,
    )
end

function _effective_config(repo_root::String, cli::CliConfig, phase_summary; run_id::String, output_dir::String)
    selected = getproperty(phase_summary, :selected)
    T_MeV = cli.T_MeV !== nothing ? cli.T_MeV :
        (selected !== nothing ? selected.T_MeV : DEFAULT_T_MEV)
    mu_min_MeV = cli.mu_min_MeV !== nothing ? cli.mu_min_MeV :
        (selected !== nothing ? selected.mu_min_MeV : DEFAULT_MU_MIN_MEV)
    mu_max_MeV = cli.mu_max_MeV !== nothing ? cli.mu_max_MeV :
        (selected !== nothing ? selected.mu_max_MeV : DEFAULT_MU_MAX_MEV)
    mu_min_MeV < mu_max_MeV || throw(ArgumentError("mu-min-MeV must be smaller than mu-max-MeV"))
    _require_positive("mu-step-MeV", cli.mu_step_MeV)
    _require_positive("p-num", cli.p_num)
    _require_positive("t-num", cli.t_num)
    _require_positive("ds-MeV", cli.ds_MeV)
    _require_positive("dsmax-MeV", cli.dsmax_MeV)
    _require_positive("max-steps", cli.max_steps)

    width = mu_max_MeV - mu_min_MeV
    start_offset = min(5.0, max(2.0, 0.05 * width))
    mu_start_MeV = mu_min_MeV + start_offset
    mu_start_MeV = min(max(mu_start_MeV, mu_min_MeV), mu_max_MeV)

    return EffectiveConfig(
        repo_root=repo_root,
        output_dir=output_dir,
        run_id=run_id,
        T_MeV=T_MeV,
        mu_min_MeV=mu_min_MeV,
        mu_max_MeV=mu_max_MeV,
        mu_step_MeV=cli.mu_step_MeV,
        mu_start_MeV=mu_start_MeV,
        xi=cli.xi,
        p_num=cli.p_num,
        t_num=cli.t_num,
        run_reference=cli.run_reference,
        model_kind=cli.model_kind,
        ds_MeV=cli.ds_MeV,
        dsmax_MeV=cli.dsmax_MeV,
        max_steps=cli.max_steps,
        branch_jump_tol=cli.branch_jump_tol,
    )
end

function _run_palc(Models, cfg::EffectiveConfig)
    model = Models.create_model(cfg.model_kind)
    T_fm = cfg.T_MeV / HBARC_MEV_FM
    mu_start_fm = cfg.mu_start_MeV / HBARC_MEV_FM
    mu_min_fm = cfg.mu_min_MeV / HBARC_MEV_FM
    mu_max_fm = cfg.mu_max_MeV / HBARC_MEV_FM
    ds_fm = cfg.ds_MeV / HBARC_MEV_FM
    dsmax_fm = cfg.dsmax_MeV / HBARC_MEV_FM

    initial_elapsed = @elapsed initial = Models.solve(
        model,
        Models.FixedMu(),
        T_fm,
        mu_start_fm;
        seed_strategy=Models.MultiSeed(),
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        residual_norm_max=1e-6,
    )
    x0 = Float64.(initial.solution)

    residual_fn = function (x, p)
        mu_vec = Main.Models.normalize_mu_vec(p.mu_fm)
        return Vector(Main.Models.gap_residual(model, x, p.T_fm, mu_vec; xi=p.xi, p_num=p.p_num, t_num=p.t_num))
    end

    jacobian_fn = (x, p) -> _finite_difference_jacobian(residual_fn, x, p)

    record_fn = function (x, mu_fm; kwargs...)
        residual_norm, pressure, masses = _safe_model_diagnostics(
            model,
            x,
            T_fm,
            Float64(mu_fm);
            xi=cfg.xi,
            p_num=cfg.p_num,
            t_num=cfg.t_num,
        )
        return (
            phi_u=Float64(x[1]),
            phi_d=Float64(x[2]),
            phi_s=Float64(x[3]),
            Phi=Float64(x[4]),
            PhiBar=Float64(x[5]),
            mu_MeV=Float64(mu_fm) * HBARC_MEV_FM,
            residual_norm=residual_norm,
            pressure_fm4=pressure,
            M_u_MeV=masses[1],
            M_d_MeV=masses[2],
            M_s_MeV=masses[3],
        )
    end

    params = (
        T_fm=T_fm,
        mu_fm=mu_start_fm,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
    )
    prob = BifurcationKit.BifurcationProblem(
        residual_fn,
        x0,
        params,
        (@optic _.mu_fm);
        J=jacobian_fn,
        record_from_solution=record_fn,
    )
    opts = BifurcationKit.ContinuationPar(
        p_min=mu_min_fm,
        p_max=mu_max_fm,
        ds=ds_fm,
        dsmax=dsmax_fm,
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
            getproperty(row, :mu_MeV),
            getproperty(row, :phi_u),
            getproperty(row, :phi_d),
            getproperty(row, :phi_s),
            getproperty(row, :Phi),
            getproperty(row, :PhiBar),
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

function _write_palc_branch_csv(path::String, br)
    header = [
        "step", "param_mu_fm", "mu_MeV",
        "phi_u", "phi_d", "phi_s", "Phi", "PhiBar",
        "residual_norm", "pressure_fm4",
        "M_u_MeV", "M_d_MeV", "M_s_MeV",
        "itnewton", "itlinear", "ds",
    ]
    return _write_csv(path, header, _branch_rows(br))
end

function _special_points(br)
    points = NamedTuple[]
    for sp in br.specialpoint
        param = try
            Float64(getproperty(sp, :param))
        catch
            NaN
        end
        push!(points, (
            source=:bifurcationkit,
            type=try getproperty(sp, :type) catch; :unknown end,
            idx=try Int(getproperty(sp, :idx)) catch; -1 end,
            step=try Int(getproperty(sp, :step)) catch; -1 end,
            mu_MeV=param * HBARC_MEV_FM,
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
    mus = Float64[getproperty(row, :mu_MeV) for row in branch]
    for i in 2:(n - 1)
        d1 = mus[i] - mus[i - 1]
        d2 = mus[i + 1] - mus[i]
        if isfinite(d1) && isfinite(d2) && d1 * d2 < 0
            push!(points, (
                source=:local_extremum,
                type=d1 > 0 ? :param_max : :param_min,
                idx=i,
                step=try Int(getproperty(branch[i], :step)) catch; i end,
                mu_MeV=mus[i],
                status=:detected,
            ))
        end
    end
    return points
end

function _state_jump_metrics(rows::Vector{<:Any}; tol::Float64)
    length(rows) < 2 && return (count=0, max_jump=0.0)
    max_jump = 0.0
    count = 0
    for i in 2:length(rows)
        a = Float64.(rows[i - 1][4:8])
        b = Float64.(rows[i][4:8])
        jump = norm(b .- a)
        isfinite(jump) || continue
        max_jump = max(max_jump, jump)
        jump > tol && (count += 1)
    end
    return (count=count, max_jump=max_jump)
end

function _palc_metrics(palc, cfg::EffectiveConfig, phase_summary)
    br = palc.result
    rows = _branch_rows(br)
    finite_residuals = [Float64(row[9]) for row in rows if isfinite(Float64(row[9]))]
    jumps = _state_jump_metrics(rows; tol=cfg.branch_jump_tol)
    special = vcat(_special_points(br), _local_param_extrema(br))
    fold_like = filter(p -> p.type in (:fold, :param_max, :param_min), special)

    hint_values = Float64[]
    selected = getproperty(phase_summary, :selected)
    if selected !== nothing
        push!(hint_values, selected.mu_spinodal_hadron_MeV)
        push!(hint_values, selected.mu_spinodal_quark_MeV)
    end

    nearest_delta = NaN
    nearest_fold_mu = NaN
    if !isempty(hint_values) && !isempty(fold_like)
        best = (delta=Inf, mu=NaN)
        for p in fold_like, hint in hint_values
            delta = abs(Float64(p.mu_MeV) - hint)
            if delta < best.delta
                best = (delta=delta, mu=Float64(p.mu_MeV))
            end
        end
        nearest_delta = best.delta
        nearest_fold_mu = best.mu
    end

    mu_values = [Float64(row[3]) for row in rows if isfinite(Float64(row[3]))]
    return (
        status=:ok,
        branch_points=length(rows),
        finite_residual_count=length(finite_residuals),
        residual_norm_max=isempty(finite_residuals) ? NaN : maximum(finite_residuals),
        residual_norm_min=isempty(finite_residuals) ? NaN : minimum(finite_residuals),
        mu_min_seen_MeV=isempty(mu_values) ? NaN : minimum(mu_values),
        mu_max_seen_MeV=isempty(mu_values) ? NaN : maximum(mu_values),
        branch_jump_count=jumps.count,
        max_state_jump=jumps.max_jump,
        newton_iterations_total=sum(Int(row[14]) for row in rows),
        linear_iterations_total=sum(Int(row[15]) for row in rows),
        special_points=special,
        fold_like_count=length(fold_like),
        nearest_fold_mu_MeV=nearest_fold_mu,
        nearest_phase_spinodal_delta_MeV=nearest_delta,
        initial_solve_wall_time_s=palc.initial_solve_wall_time_s,
        continuation_wall_time_s=palc.continuation_wall_time_s,
        wall_time_s=palc.wall_time_s,
    )
end

function _read_scan_csv(path::String)
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
            mu_MeV=parse(Float64, getcol("mu_MeV")),
            xi=parse(Float64, getcol("xi")),
            phi_u=parse(Float64, getcol("phi_u")),
            phi_d=parse(Float64, getcol("phi_d")),
            phi_s=parse(Float64, getcol("phi_s")),
            Phi=parse(Float64, getcol("Phi1")),
            PhiBar=parse(Float64, getcol("Phi2")),
            residual_norm=parse(Float64, getcol("residual_norm")),
            iterations=parse(Int, getcol("iterations", "-1")),
            converged=lowercase(strip(getcol("converged", "false"))) == "true",
        ))
    end
    return rows
end

function _run_tmu_reference(Models, cfg::EffectiveConfig)
    path = joinpath(cfg.output_dir, "tmu_scan_reference.csv")
    mu_values = collect(cfg.mu_min_MeV:cfg.mu_step_MeV:cfg.mu_max_MeV)
    if isempty(mu_values) || last(mu_values) < cfg.mu_max_MeV
        push!(mu_values, cfg.mu_max_MeV)
    end
    elapsed = @elapsed stats = Models.run_tmu_scan(
        T_values=[cfg.T_MeV],
        mu_values=mu_values,
        xi_values=[cfg.xi],
        output_path=path,
        overwrite=true,
        resume=false,
        use_phase_aware=true,
        bootstrap_multiseed=true,
        solver_backend=:models,
        model_kind=cfg.model_kind,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
    )
    rows = _read_scan_csv(path)
    model = Models.create_model(cfg.model_kind)
    diagnostic_residuals = Float64[]
    for row in rows
        state = [row.phi_u, row.phi_d, row.phi_s, row.Phi, row.PhiBar]
        residual_norm, _, _ = _safe_model_diagnostics(
            model,
            state,
            row.T_MeV / HBARC_MEV_FM,
            row.mu_MeV / HBARC_MEV_FM;
            xi=cfg.xi,
            p_num=cfg.p_num,
            t_num=cfg.t_num,
        )
        isfinite(residual_norm) && push!(diagnostic_residuals, residual_norm)
    end
    jumps = _state_jump_metrics([Any[
            row.T_MeV, row.mu_MeV, row.xi,
            row.phi_u, row.phi_d, row.phi_s, row.Phi, row.PhiBar,
        ] for row in rows]; tol=cfg.branch_jump_tol)
    return (
        path=path,
        wall_time_s=elapsed,
        stats=stats,
        points=length(rows),
        failure_count=count(row -> !row.converged, rows),
        branch_jump_count=jumps.count,
        max_state_jump=jumps.max_jump,
        finite_residual_count=length(diagnostic_residuals),
        residual_norm_min=isempty(diagnostic_residuals) ? NaN : minimum(diagnostic_residuals),
        residual_norm_max=isempty(diagnostic_residuals) ? NaN : maximum(diagnostic_residuals),
        iterations_total=sum(max(row.iterations, 0) for row in rows),
    )
end

function _skipped_tmu_reference()
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

function _continuation_decision(palc_metrics, phase_summary)
    has_finite_branch = palc_metrics.branch_points > 0 && palc_metrics.finite_residual_count > 0
    has_phase_hint = getproperty(phase_summary, :selected) !== nothing
    fold_close = isfinite(palc_metrics.nearest_phase_spinodal_delta_MeV) &&
        palc_metrics.nearest_phase_spinodal_delta_MeV <= FOLD_HINT_TARGET_MEV
    recommended = has_finite_branch && has_phase_hint && fold_close
    reason = if !has_finite_branch
        "PALC produced no finite branch rows"
    elseif has_phase_hint && !fold_close
        "nearest PALC fold-like point is farther than $(FOLD_HINT_TARGET_MEV) MeV from phase spinodal hint"
    elseif !has_phase_hint
        "PALC produced finite branch rows, but phase spinodal hint was unavailable; continuation decision remains incomplete"
    else
        "PALC produced finite branch rows and a fold-like point near the phase spinodal hint"
    end
    return (recommended=recommended, reason=reason)
end

function _write_report(path::String, cfg::EffectiveConfig, phase_summary, palc_metrics, tmu_metrics, decision)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# FixedMu PNJL PALC Spike Report")
        println(io)
        println(io, "## Configuration")
        println(io, "- T_MeV: $(cfg.T_MeV)")
        println(io, "- mu window MeV: $(cfg.mu_min_MeV) to $(cfg.mu_max_MeV)")
        println(io, "- xi: $(cfg.xi)")
        println(io, "- p_num/t_num: $(cfg.p_num)/$(cfg.t_num)")
        println(io, "- run_reference: $(cfg.run_reference)")
        println(io)
        println(io, "## Phase Reference")
        println(io, "- ran: $(phase_summary.ran)")
        println(io, "- spinodal_count: $(phase_summary.spinodal_count)")
        if phase_summary.selected !== nothing
            sel = phase_summary.selected
            println(io, "- selected T_MeV: $(sel.T_MeV)")
            println(io, "- spinodal hadron/quark mu_MeV: $(sel.mu_spinodal_hadron_MeV), $(sel.mu_spinodal_quark_MeV)")
        else
            println(io, "- selected: none ($(phase_summary.fallback_reason))")
        end
        println(io)
        println(io, "## PALC Branch")
        println(io, "- branch_points: $(palc_metrics.branch_points)")
        println(io, "- finite_residual_count: $(palc_metrics.finite_residual_count)")
        println(io, "- residual_norm_max: $(_fmt(palc_metrics.residual_norm_max))")
        println(io, "- mu_seen_MeV: $(_fmt(palc_metrics.mu_min_seen_MeV)) to $(_fmt(palc_metrics.mu_max_seen_MeV))")
        println(io, "- branch_jump_count: $(palc_metrics.branch_jump_count)")
        println(io, "- max_state_jump: $(_fmt(palc_metrics.max_state_jump))")
        println(io, "- newton_iterations_total: $(palc_metrics.newton_iterations_total)")
        println(io, "- wall_time_s: $(_fmt(palc_metrics.wall_time_s))")
        println(io, "- fold_like_count: $(palc_metrics.fold_like_count)")
        println(io, "- nearest_phase_spinodal_delta_MeV: $(_fmt(palc_metrics.nearest_phase_spinodal_delta_MeV))")
        println(io)
        println(io, "## T-mu Reference")
        println(io, "- points: $(tmu_metrics.points)")
        println(io, "- failure_count: $(tmu_metrics.failure_count)")
        println(io, "- branch_jump_count: $(tmu_metrics.branch_jump_count)")
        println(io, "- max_state_jump: $(_fmt(tmu_metrics.max_state_jump))")
        println(io, "- finite_residual_count: $(tmu_metrics.finite_residual_count)")
        println(io, "- residual_norm_min/max: $(_fmt(tmu_metrics.residual_norm_min)) / $(_fmt(tmu_metrics.residual_norm_max))")
        println(io, "- iterations_total: $(tmu_metrics.iterations_total)")
        println(io, "- wall_time_s: $(_fmt(tmu_metrics.wall_time_s))")
        println(io)
        println(io, "## Decision")
        println(io, "- continue_recommended: $(decision.recommended)")
        println(io, "- reason: $(decision.reason)")
    end
    return path
end

function _write_outputs(cfg::EffectiveConfig, phase_summary, palc, tmu_metrics)
    mkpath(cfg.output_dir)
    branch_path = joinpath(cfg.output_dir, "palc_fixedmu_branch.csv")
    phase_path = joinpath(cfg.output_dir, "phase_reference_summary.json")
    summary_path = joinpath(cfg.output_dir, "comparison_summary.json")
    report_path = joinpath(cfg.output_dir, "comparison_report.md")

    _write_palc_branch_csv(branch_path, palc.result)
    palc_metrics = _palc_metrics(palc, cfg, phase_summary)
    decision = _continuation_decision(palc_metrics, phase_summary)

    _write_json(phase_path, phase_summary)
    summary = (
        config=cfg,
        artifacts=(
            palc_branch=branch_path,
            phase_reference_summary=phase_path,
            tmu_scan_reference=tmu_metrics.path,
            report=report_path,
        ),
        phase_reference=phase_summary,
        palc=palc_metrics,
        tmu_reference=tmu_metrics,
        decision=decision,
    )
    _write_json(summary_path, summary)
    _write_report(report_path, cfg, phase_summary, palc_metrics, tmu_metrics, decision)
    return summary
end

function run(args::Vector{String}; repo_root::String, default_run_reference::Bool)
    cli = _parse_cli(args; default_run_reference=default_run_reference)
    Models = _ensure_models_loaded(repo_root)

    provisional_run_id = _run_id()
    provisional_output_dir = cli.output_dir === nothing ?
        joinpath(repo_root, DEFAULT_OUTPUT_REL, provisional_run_id) :
        (isabspath(cli.output_dir) ? cli.output_dir : joinpath(repo_root, cli.output_dir))
    mkpath(provisional_output_dir)

    phase_summary = cli.run_reference ?
        _run_phase_reference(Models, cli, provisional_output_dir) :
        _skipped_phase_reference("run_reference=false")
    cfg = _effective_config(repo_root, cli, phase_summary; run_id=provisional_run_id, output_dir=provisional_output_dir)

    palc = _run_palc(Models, cfg)
    tmu_metrics = cfg.run_reference ? _run_tmu_reference(Models, cfg) : _skipped_tmu_reference()
    summary = _write_outputs(cfg, phase_summary, palc, tmu_metrics)

    println("PALC spike output: $(cfg.output_dir)")
    println("continue_recommended: $(summary.decision.recommended)")
    println("reason: $(summary.decision.reason)")
    return summary
end

function main_run(args::Vector{String}; repo_root::String, default_run_reference::Bool)
    try
        run(args; repo_root=repo_root, default_run_reference=default_run_reference)
    catch err
        println(stderr, "FixedMu PALC spike failed.")
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        exit(1)
    end
    return nothing
end

end # module
