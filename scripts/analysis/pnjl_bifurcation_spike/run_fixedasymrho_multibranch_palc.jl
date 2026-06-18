const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !(REPO_ROOT in LOAD_PATH)
    insert!(LOAD_PATH, min(2, length(LOAD_PATH) + 1), REPO_ROOT)
end
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(REPO_ROOT, "src", "models", "Models.jl"))
end

include(joinpath(@__DIR__, "fixedmu_palc_common.jl"))
include(joinpath(@__DIR__, "fixedasymrho_palc_common.jl"))

module FixedAsymRhoMultiBranchPALCSpike

using Dates
using LinearAlgebra

try
    @eval import BifurcationKit
    @eval import BifurcationKit: @optic
catch err
    error("BifurcationKit is not available. Run with --project=scripts/analysis/pnjl_bifurcation_spike and instantiate that environment. Original error: $(err)")
end

const U = Main.FixedMuPALCSpike
const S = Main.FixedAsymRhoPALCSpike
const HBARC_MEV_FM = U.HBARC_MEV_FM
const DEFAULT_OUTPUT_REL = joinpath("data", "outputs", "results", "analysis", "palc_fixedasymrho_spike")

Base.@kwdef struct MultiBranchConfig
    repo_root::String
    output_dir::String
    run_id::String
    T_MeV::Float64 = 120.0
    rho_min::Float64 = 0.25
    rho_max::Float64 = 0.45
    rho_anchor::Float64 = 0.35
    rho_step::Float64 = 0.05
    production_source_rho::Float64 = 1.0
    production_source_step::Float64 = 0.05
    asym_ud_ratio_target::Float64 = 0.876
    asym_s_target::Float64 = 0.0
    xi::Float64 = 0.0
    p_num::Int = 8
    t_num::Int = 4
    ds_rho::Float64 = 0.01
    dsmax_rho::Float64 = 0.04
    max_steps::Int = 80
    branch_jump_tol::Float64 = S.DEFAULT_BRANCH_JUMP_TOL
    root_distance_tol::Float64 = 0.25
    pressure_gap_tol::Float64 = S.DEFAULT_PRESSURE_GAP_TOL
end

function _usage()
    return """
    Usage:
      julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/run_fixedasymrho_multibranch_palc.jl [options]

    Options:
      --T-MeV=<value>                    Fixed temperature in MeV. Default: 120.
      --rho-anchor=<value>               Root-discovery anchor rho/rho0. Default: 0.35.
      --rho-min=<value>                  Lower PALC rho/rho0 bound. Default: 0.25.
      --rho-max=<value>                  Upper PALC rho/rho0 bound. Default: 0.45.
      --rho-step=<value>                 Ground-state comparison rho/rho0 step. Default: 0.05.
      --production-source-rho=<value>    Starting rho/rho0 for production-like continuation seed. Default: 1.0.
      --production-source-step=<value>   Step size from production source to anchor. Default: 0.05.
      --asym-ud-ratio-target=<value>     Target rho_u/rho_d. Default: 0.876.
      --asym-s-target=<value>            Target rho_s in fm^-3. Default: 0.
      --xi=<value>                       Anisotropy parameter. Default: 0.
      --p-num=<int>                      Momentum nodes. Default: 8.
      --t-num=<int>                      Theta nodes. Default: 4.
      --ds-rho=<value>                   PALC initial pseudo-arclength scale. Default: 0.01.
      --dsmax-rho=<value>                PALC max pseudo-arclength scale. Default: 0.04.
      --max-steps=<int>                  PALC max continuation steps. Default: 80.
      --root-distance-tol=<value>        Natural-unit root dedup tolerance. Default: 0.25.
      --pressure-gap-tol=<value>         Pressure selection gap tolerance. Default: 1e-3.
      --output-dir=<path>                Output directory.
    """
end

function _parse_args(args::Vector{String}; repo_root::String)
    raw = Dict{String, String}()
    for arg in args
        arg in ("-h", "--help") && (println(_usage()); exit(0))
        startswith(arg, "--") || throw(ArgumentError("unexpected positional argument: $(arg)"))
        parts = split(arg[3:end], "="; limit=2)
        length(parts) == 2 || throw(ArgumentError("expected --key=value, got $(arg)"))
        raw[replace(lowercase(parts[1]), "-" => "_")] = parts[2]
    end

    get_float(key, default) = haskey(raw, key) ? parse(Float64, raw[key]) : default
    get_int(key, default) = haskey(raw, key) ? parse(Int, raw[key]) : default

    run_id = U._run_id()
    output_dir = get(raw, "output_dir", nothing)
    output_dir = output_dir === nothing ?
        joinpath(repo_root, DEFAULT_OUTPUT_REL, "multibranch_$(run_id)") :
        (isabspath(output_dir) ? output_dir : joinpath(repo_root, output_dir))

    cfg = MultiBranchConfig(
        repo_root=repo_root,
        output_dir=output_dir,
        run_id=run_id,
        T_MeV=get_float("t_mev", 120.0),
        rho_min=get_float("rho_min", 0.25),
        rho_max=get_float("rho_max", 0.45),
        rho_anchor=get_float("rho_anchor", 0.35),
        rho_step=get_float("rho_step", 0.05),
        production_source_rho=get_float("production_source_rho", 1.0),
        production_source_step=get_float("production_source_step", 0.05),
        asym_ud_ratio_target=get_float("asym_ud_ratio_target", 0.876),
        asym_s_target=get_float("asym_s_target", 0.0),
        xi=get_float("xi", 0.0),
        p_num=get_int("p_num", 8),
        t_num=get_int("t_num", 4),
        ds_rho=get_float("ds_rho", 0.01),
        dsmax_rho=get_float("dsmax_rho", 0.04),
        max_steps=get_int("max_steps", 80),
        branch_jump_tol=get_float("branch_jump_tol", S.DEFAULT_BRANCH_JUMP_TOL),
        root_distance_tol=get_float("root_distance_tol", 0.25),
        pressure_gap_tol=get_float("pressure_gap_tol", S.DEFAULT_PRESSURE_GAP_TOL),
    )
    _validate_config(cfg)
    return cfg
end

function _validate_config(cfg::MultiBranchConfig)
    cfg.rho_min < cfg.rho_max || throw(ArgumentError("rho-min must be smaller than rho-max"))
    cfg.rho_min <= cfg.rho_anchor <= cfg.rho_max || throw(ArgumentError("rho-anchor must be inside rho-min..rho-max"))
    U._require_positive("rho-step", cfg.rho_step)
    U._require_positive("production-source-step", cfg.production_source_step)
    U._require_positive("p-num", cfg.p_num)
    U._require_positive("t-num", cfg.t_num)
    U._require_positive("ds-rho", cfg.ds_rho)
    U._require_positive("dsmax-rho", cfg.dsmax_rho)
    U._require_positive("max-steps", cfg.max_steps)
    U._require_positive("root-distance-tol", cfg.root_distance_tol)
    U._require_positive("pressure-gap-tol", cfg.pressure_gap_tol)
    return nothing
end

function _mode(cfg::MultiBranchConfig, rho_target::Real)
    return Main.Models.FixedAsymmetricRho(
        Float64(rho_target),
        cfg.asym_ud_ratio_target,
        cfg.asym_s_target,
    )
end

function _effective_cfg(cfg::MultiBranchConfig)
    return S.EffectiveConfig(
        repo_root=cfg.repo_root,
        output_dir=cfg.output_dir,
        run_id=cfg.run_id,
        T_MeV=cfg.T_MeV,
        rho_min=cfg.rho_min,
        rho_max=cfg.rho_max,
        rho_step=cfg.rho_step,
        rho_start=cfg.rho_anchor,
        asym_ud_ratio_target=cfg.asym_ud_ratio_target,
        asym_s_target=cfg.asym_s_target,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        run_reference=false,
        model_kind=:PNJL,
        ds_rho=cfg.ds_rho,
        dsmax_rho=cfg.dsmax_rho,
        max_steps=cfg.max_steps,
        branch_jump_tol=cfg.branch_jump_tol,
        trho_reverse_rho=true,
        pressure_gap_tol=cfg.pressure_gap_tol,
    )
end

function _solve_kwargs(cfg::MultiBranchConfig)
    return (
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        iterations=300,
        residual_norm_max=1e-3,
    )
end

function _root_record(model, cfg::MultiBranchConfig, source::Symbol, result)
    ecfg = _effective_cfg(cfg)
    residual_norm, pressure, rho_vec, rho_norm, ud_ratio, masses =
        S._safe_diagnostics(model, ecfg, result.solution[1:8], cfg.rho_anchor)
    return (
        source=source,
        solution=Float64.(result.solution[1:8]),
        converged=Bool(result.converged),
        residual_norm=Float64(residual_norm),
        pressure_fm4=Float64(pressure),
        rho_norm=Float64(rho_norm),
        rho_u=Float64(rho_vec[1]),
        rho_d=Float64(rho_vec[2]),
        rho_s=Float64(rho_vec[3]),
        ud_ratio=Float64(ud_ratio),
        mu_u_MeV=Float64(result.mu_vec[1]) * HBARC_MEV_FM,
        mu_d_MeV=Float64(result.mu_vec[2]) * HBARC_MEV_FM,
        mu_s_MeV=Float64(result.mu_vec[3]) * HBARC_MEV_FM,
        M_u_MeV=Float64(masses[1]),
        M_d_MeV=Float64(masses[2]),
        M_s_MeV=Float64(masses[3]),
    )
end

function _state_distance(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    return norm(Float64.(a[1:8]) .- Float64.(b[1:8]))
end

function _push_distinct_root!(roots::Vector{NamedTuple}, candidate, cfg::MultiBranchConfig)
    candidate.converged || return false
    isfinite(candidate.residual_norm) || return false
    for root in roots
        if _state_distance(root.solution, candidate.solution) <= cfg.root_distance_tol
            return false
        end
    end
    push!(roots, candidate)
    return true
end

function _rho_path(source_rho::Float64, anchor_rho::Float64, step::Float64)
    if isapprox(source_rho, anchor_rho; atol=1e-12)
        return [anchor_rho]
    end
    direction = source_rho > anchor_rho ? -1.0 : 1.0
    vals = Float64[]
    rho = source_rho
    while direction > 0 ? rho < anchor_rho - 1e-9 : rho > anchor_rho + 1e-9
        push!(vals, rho)
        rho += direction * step
    end
    push!(vals, anchor_rho)
    return vals
end

function _production_like_anchor_root(model, cfg::MultiBranchConfig, T_fm::Float64)
    seed = nothing
    result = nothing
    for rho in _rho_path(cfg.production_source_rho, cfg.rho_anchor, cfg.production_source_step)
        kwargs = _solve_kwargs(cfg)
        result = if seed === nothing
            Main.Models.solve(
                model,
                _mode(cfg, rho),
                T_fm;
                seed_strategy=Main.Models.MultiSeed(),
                kwargs...,
            )
        else
            Main.Models.solve(
                model,
                _mode(cfg, rho),
                T_fm;
                seed_guess=seed,
                kwargs...,
            )
        end
        seed = Float64.(result.solution[1:8])
    end
    return result
end

function _discover_anchor_roots(model, cfg::MultiBranchConfig)
    T_fm = cfg.T_MeV / HBARC_MEV_FM
    mode = _mode(cfg, cfg.rho_anchor)
    roots = NamedTuple[]
    attempts = NamedTuple[]

    multiseeds = Main.Models.get_all_seeds(Main.Models.MultiSeed(), [T_fm], mode)
    push!(attempts, (source=:multiseed_catalog, seed_count=length(multiseeds)))
    pressure_root = Main.Models.solve_multi(
        model,
        mode,
        T_fm;
        seeds=multiseeds,
        semantic_mode=:ground_state,
        evaluate_all_attempts=true,
        _solve_kwargs(cfg)...,
    )
    _push_distinct_root!(roots, _root_record(model, cfg, :solve_multi_pressure, pressure_root), cfg)

    production_root = _production_like_anchor_root(model, cfg, T_fm)
    _push_distinct_root!(roots, _root_record(model, cfg, :production_like_continuation, production_root), cfg)

    sort!(roots; by=r -> -r.pressure_fm4)
    return roots, attempts
end

function _run_palc_branch(model, cfg::MultiBranchConfig, branch_id::String, root)
    T_fm = cfg.T_MeV / HBARC_MEV_FM
    ecfg = _effective_cfg(cfg)
    residual_fn = S._residual_fn(model, ecfg)
    jacobian_fn = (x, p) -> U._finite_difference_jacobian(residual_fn, x, p)
    record_fn = function (x, rho_target; kwargs...)
        residual_norm, pressure, rho_vec, rho_norm, ud_ratio, masses =
            S._safe_diagnostics(model, ecfg, x, Float64(rho_target))
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
        rho_target=cfg.rho_anchor,
        asym_ud_ratio_target=cfg.asym_ud_ratio_target,
        asym_s_target=cfg.asym_s_target,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
    )
    prob = BifurcationKit.BifurcationProblem(
        residual_fn,
        root.solution,
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
            tol=S.PALC_NEWTON_TOL,
            max_iterations=S.PALC_NEWTON_MAX_ITERATIONS,
            verbose=false,
        ),
    )

    branch_ref = Ref{Any}(nothing)
    elapsed = @elapsed begin
        branch_ref[] = BifurcationKit.continuation(prob, BifurcationKit.PALC(), opts; bothside=true, verbosity=0)
    end
    return (branch_id=branch_id, root=root, result=branch_ref[], wall_time_s=elapsed)
end

function _branch_rows(branches)
    rows = Vector{Vector}()
    for branch in branches
        for row in branch.result.branch
            push!(rows, Any[
                branch.branch_id,
                branch.root.source,
                branch.root.pressure_fm4,
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
    end
    return rows
end

function _write_multibranch_csv(path::String, branches)
    header = [
        "branch_id", "anchor_source", "anchor_pressure_fm4",
        "step", "param_rho_target", "rho_target", "rho_norm",
        "phi_u", "phi_d", "phi_s", "Phi", "PhiBar",
        "mu_u_MeV", "mu_d_MeV", "mu_s_MeV",
        "rho_u", "rho_d", "rho_s", "ud_ratio",
        "residual_norm", "pressure_fm4",
        "M_u_MeV", "M_d_MeV", "M_s_MeV",
        "itnewton", "itlinear", "ds",
    ]
    return U._write_csv(path, header, _branch_rows(branches))
end

function _write_anchor_csv(path::String, roots)
    header = [
        "branch_id", "source", "pressure_fm4", "residual_norm",
        "rho_norm", "rho_u", "rho_d", "rho_s", "ud_ratio",
        "mu_u_MeV", "mu_d_MeV", "mu_s_MeV",
        "phi_u", "phi_d", "phi_s", "Phi", "PhiBar",
        "M_u_MeV", "M_d_MeV", "M_s_MeV",
    ]
    rows = Vector{Vector}()
    for (idx, root) in enumerate(roots)
        push!(rows, Any[
            "branch_$(idx)",
            root.source,
            root.pressure_fm4,
            root.residual_norm,
            root.rho_norm,
            root.rho_u,
            root.rho_d,
            root.rho_s,
            root.ud_ratio,
            root.mu_u_MeV,
            root.mu_d_MeV,
            root.mu_s_MeV,
            root.solution[1],
            root.solution[2],
            root.solution[3],
            root.solution[4],
            root.solution[5],
            root.M_u_MeV,
            root.M_d_MeV,
            root.M_s_MeV,
        ])
    end
    return U._write_csv(path, header, rows)
end

function _nearest_branch_row(branch, rho::Float64)
    best = nothing
    best_delta = Inf
    for row in branch.result.branch
        delta = abs(Float64(getproperty(row, :rho_target)) - rho)
        if isfinite(delta) && delta < best_delta
            best = row
            best_delta = delta
        end
    end
    return best === nothing ? nothing : (row=best, rho_delta=best_delta)
end

function _rho_grid(cfg::MultiBranchConfig)
    vals = Float64[]
    rho = cfg.rho_min
    while rho <= cfg.rho_max + 1e-9
        push!(vals, rho)
        rho += cfg.rho_step
    end
    if isempty(vals) || abs(last(vals) - cfg.rho_max) > 1e-9
        push!(vals, cfg.rho_max)
    end
    return vals
end

function _groundstate_rows(branches, cfg::MultiBranchConfig)
    rows = Vector{Vector}()
    for rho in _rho_grid(cfg)
        candidates = NamedTuple[]
        for branch in branches
            nearest = _nearest_branch_row(branch, rho)
            nearest === nothing && continue
            nearest.rho_delta <= max(cfg.dsmax_rho, cfg.rho_step) || continue
            row = nearest.row
            pressure = Float64(getproperty(row, :pressure_fm4))
            residual = Float64(getproperty(row, :residual_norm))
            isfinite(pressure) && isfinite(residual) || continue
            push!(candidates, (
                branch_id=String(branch.branch_id),
                rho=Float64(getproperty(row, :rho_target)),
                rho_delta=Float64(nearest.rho_delta),
                pressure_fm4=pressure,
                residual_norm=residual,
                mu_u_MeV=Float64(getproperty(row, :mu_u_MeV)),
                mu_d_MeV=Float64(getproperty(row, :mu_d_MeV)),
                mu_s_MeV=Float64(getproperty(row, :mu_s_MeV)),
            ))
        end
        isempty(candidates) && continue
        sort!(candidates; by=c -> -c.pressure_fm4)
        best = first(candidates)
        runner_up = length(candidates) >= 2 ? candidates[2] : nothing
        gap = runner_up === nothing ? NaN : best.pressure_fm4 - runner_up.pressure_fm4
        push!(rows, Any[
            rho,
            best.branch_id,
            best.rho,
            best.rho_delta,
            best.pressure_fm4,
            runner_up === nothing ? "" : runner_up.branch_id,
            runner_up === nothing ? "" : runner_up.pressure_fm4,
            gap,
            best.residual_norm,
            best.mu_u_MeV,
            best.mu_d_MeV,
            best.mu_s_MeV,
            length(candidates),
        ])
    end
    return rows
end

function _write_groundstate_csv(path::String, branches, cfg::MultiBranchConfig)
    header = [
        "rho_target", "selected_branch_id", "nearest_branch_rho",
        "rho_delta", "selected_pressure_fm4",
        "runner_up_branch_id", "runner_up_pressure_fm4", "pressure_gap_fm4",
        "selected_residual_norm",
        "selected_mu_u_MeV", "selected_mu_d_MeV", "selected_mu_s_MeV",
        "candidate_branch_count",
    ]
    return U._write_csv(path, header, _groundstate_rows(branches, cfg))
end

function _branch_summary(branch, cfg::MultiBranchConfig)
    rows = _branch_rows((branch,))
    residuals = [Float64(row[20]) for row in rows if isfinite(Float64(row[20]))]
    pressures = [Float64(row[21]) for row in rows if isfinite(Float64(row[21]))]
    rhos = [Float64(row[6]) for row in rows if isfinite(Float64(row[6]))]
    jumps = S._state_jump_metrics(rows; tol=cfg.branch_jump_tol, cols=8:15)
    return (
        branch_id=branch.branch_id,
        anchor_source=branch.root.source,
        anchor_pressure_fm4=branch.root.pressure_fm4,
        points=length(rows),
        rho_min_seen=isempty(rhos) ? NaN : minimum(rhos),
        rho_max_seen=isempty(rhos) ? NaN : maximum(rhos),
        pressure_min=isempty(pressures) ? NaN : minimum(pressures),
        pressure_max=isempty(pressures) ? NaN : maximum(pressures),
        residual_norm_max=isempty(residuals) ? NaN : maximum(residuals),
        branch_jump_count=jumps.count,
        max_state_jump=jumps.max_jump,
        wall_time_s=branch.wall_time_s,
    )
end

function _write_report(path::String, cfg::MultiBranchConfig, roots, branch_summaries, groundstate_rows)
    mkpath(dirname(path))
    anchor_pressure_gap = length(roots) >= 2 ? roots[1].pressure_fm4 - roots[2].pressure_fm4 : NaN
    experimental_candidate = length(roots) >= 2 && !isempty(groundstate_rows)
    open(path, "w") do io
        println(io, "# FixedAsymmetricRho Multi-Branch PALC Spike")
        println(io)
        println(io, "## Configuration")
        println(io, "- T_MeV: $(cfg.T_MeV)")
        println(io, "- rho window: $(cfg.rho_min) to $(cfg.rho_max)")
        println(io, "- rho_anchor: $(cfg.rho_anchor)")
        println(io, "- production_source_rho: $(cfg.production_source_rho)")
        println(io, "- asym_ud_ratio_target: $(cfg.asym_ud_ratio_target)")
        println(io, "- asym_s_target: $(cfg.asym_s_target)")
        println(io, "- p_num/t_num: $(cfg.p_num)/$(cfg.t_num)")
        println(io)
        println(io, "## Anchor Roots")
        println(io, "- distinct_root_count: $(length(roots))")
        println(io, "- anchor_pressure_gap_fm4: $(U._fmt(anchor_pressure_gap))")
        for (idx, root) in enumerate(roots)
            println(io, "- branch_$(idx): source=$(root.source), pressure_fm4=$(U._fmt(root.pressure_fm4)), mu_MeV=($(U._fmt(root.mu_u_MeV)), $(U._fmt(root.mu_d_MeV)), $(U._fmt(root.mu_s_MeV)))")
        end
        println(io)
        println(io, "## PALC Branches")
        for summary in branch_summaries
            println(io, "- $(summary.branch_id): points=$(summary.points), rho_seen=$(U._fmt(summary.rho_min_seen))..$(U._fmt(summary.rho_max_seen)), pressure=$(U._fmt(summary.pressure_min))..$(U._fmt(summary.pressure_max)), residual_max=$(U._fmt(summary.residual_norm_max)), jumps=$(summary.branch_jump_count), wall_time_s=$(U._fmt(summary.wall_time_s))")
        end
        println(io)
        println(io, "## Ground-State Selection")
        println(io, "- sampled_rho_count: $(length(groundstate_rows))")
        if !isempty(groundstate_rows)
            first_row = first(groundstate_rows)
            println(io, "- first_selection: rho=$(first_row[1]), branch=$(first_row[2]), pressure_fm4=$(U._fmt(first_row[5])), candidates=$(first_row[13])")
        end
        println(io)
        println(io, "## Decision")
        println(io, "- experimental_backend_candidate: $(experimental_candidate)")
        println(io, "- reason: multi-branch PALC is useful only if anchor discovery finds distinct roots and branch-level pressure selection remains explicit.")
    end
    return path
end

function _config_dict(cfg::MultiBranchConfig)
    return Dict(
        "repo_root" => cfg.repo_root,
        "output_dir" => cfg.output_dir,
        "run_id" => cfg.run_id,
        "T_MeV" => cfg.T_MeV,
        "rho_min" => cfg.rho_min,
        "rho_max" => cfg.rho_max,
        "rho_anchor" => cfg.rho_anchor,
        "rho_step" => cfg.rho_step,
        "production_source_rho" => cfg.production_source_rho,
        "production_source_step" => cfg.production_source_step,
        "asym_ud_ratio_target" => cfg.asym_ud_ratio_target,
        "asym_s_target" => cfg.asym_s_target,
        "xi" => cfg.xi,
        "p_num" => cfg.p_num,
        "t_num" => cfg.t_num,
        "ds_rho" => cfg.ds_rho,
        "dsmax_rho" => cfg.dsmax_rho,
        "max_steps" => cfg.max_steps,
        "branch_jump_tol" => cfg.branch_jump_tol,
        "root_distance_tol" => cfg.root_distance_tol,
        "pressure_gap_tol" => cfg.pressure_gap_tol,
    )
end

function run(args::Vector{String}; repo_root::String)
    cfg = _parse_args(args; repo_root=repo_root)
    model = Main.Models.create_model(:PNJL)
    mkpath(cfg.output_dir)

    roots, attempts = _discover_anchor_roots(model, cfg)
    isempty(roots) && error("anchor discovery produced no converged distinct roots")
    branches = NamedTuple[]
    for (idx, root) in enumerate(roots)
        branch_id = "branch_$(idx)"
        push!(branches, _run_palc_branch(model, cfg, branch_id, root))
    end

    anchor_path = joinpath(cfg.output_dir, "palc_fixedasymrho_anchor_roots.csv")
    branch_path = joinpath(cfg.output_dir, "palc_fixedasymrho_multibranch.csv")
    groundstate_path = joinpath(cfg.output_dir, "palc_fixedasymrho_groundstate_selection.csv")
    summary_path = joinpath(cfg.output_dir, "multibranch_summary.json")
    report_path = joinpath(cfg.output_dir, "multibranch_report.md")

    _write_anchor_csv(anchor_path, roots)
    _write_multibranch_csv(branch_path, branches)
    _write_groundstate_csv(groundstate_path, branches, cfg)
    branch_summaries = [_branch_summary(branch, cfg) for branch in branches]
    ground_rows = _groundstate_rows(branches, cfg)
    _write_report(report_path, cfg, roots, branch_summaries, ground_rows)

    summary = (
        config=_config_dict(cfg),
        artifacts=(
            anchor_roots=anchor_path,
            multibranch=branch_path,
            groundstate_selection=groundstate_path,
            report=report_path,
        ),
        anchor_attempts=attempts,
        anchor_roots=[(
            branch_id="branch_$(idx)",
            source=root.source,
            pressure_fm4=root.pressure_fm4,
            residual_norm=root.residual_norm,
            mu_u_MeV=root.mu_u_MeV,
            mu_d_MeV=root.mu_d_MeV,
            mu_s_MeV=root.mu_s_MeV,
        ) for (idx, root) in enumerate(roots)],
        branch_summaries=branch_summaries,
        groundstate_sample_count=length(ground_rows),
        experimental_backend_candidate=(length(roots) >= 2 && !isempty(ground_rows)),
    )
    U._write_json(summary_path, summary)

    println("FixedAsymmetricRho multi-branch PALC output: $(cfg.output_dir)")
    println("distinct_root_count: $(length(roots))")
    println("experimental_backend_candidate: $(summary.experimental_backend_candidate)")
    return summary
end

function main_run(args::Vector{String}; repo_root::String)
    try
        run(args; repo_root=repo_root)
    catch err
        println(stderr, "FixedAsymmetricRho multi-branch PALC spike failed.")
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        exit(1)
    end
    return nothing
end

end # module

using .FixedAsymRhoMultiBranchPALCSpike

FixedAsymRhoMultiBranchPALCSpike.main_run(ARGS; repo_root=REPO_ROOT)
