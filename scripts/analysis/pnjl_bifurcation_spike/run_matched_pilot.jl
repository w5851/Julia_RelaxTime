const MATCHED_PILOT_LOAD_STARTED = time()
include(joinpath(@__DIR__, "run_fixedasymrho_multibranch_palc.jl"))
include(joinpath(@__DIR__, "matched_pilot_support.jl"))

module MatchedFixedAsymRhoPilot

using LinearAlgebra
import JSON
import NLsolve
import BifurcationKit
import BifurcationKit: @optic

const Models = Main.Models
const M = Main.FixedAsymRhoMultiBranchPALCSpike
const S = Main.FixedAsymRhoPALCSpike
const U = Main.FixedMuPALCSpike
const Support = Main.MatchedPALCPilotSupport
include(joinpath(@__DIR__, "..", "..", "relaxtime", "combined_meson_density_equilibrium.jl"))

const RESIDUAL_TOL = 1e-6
const ROOT_TOL = 1e-9
const ROOT_DISTANCE = 1e-4
const METHODS = ("palc", "seed", "multi", "production", "dense", "legacy_seed", "seed_reverse")

Base.@kwdef mutable struct PilotWork
    solver::Models.SolverWorkTelemetry = Models.SolverWorkTelemetry()
    bk::Support.NewtonWork = Support.NewtonWork()
    direct_solves::Int = 0
    direct_residual_evaluations::Int = 0
    direct_nlsolve_f_calls::Int = 0
    direct_nlsolve_g_calls::Int = 0
    direct_iterations::Int = 0
    direct_nonconverged::Int = 0
    direct_exceptions::Int = 0
    bk_residual_evaluations::Int = 0
    bk_jacobian_calls::Int = 0
    bk_exceptions::Int = 0
    diagnostic_residual_calls::Int = 0
    thermo_evaluations::Int = 0
    anchor_requests::Int = 0
    anchor_count::Int = 0
    branch_count::Int = 0
    continuation_points::Int = 0
    correction_requests::Int = 0
end

function config(case, output; profile=get(ENV, "PALC_PILOT_PROFILE", "local"))
    case in ("t120_rho035", "t130_rho080") || throw(ArgumentError("unknown case: $case"))
    profile in ("local", "long", "higher", "transition") || throw(ArgumentError("unknown profile: $profile"))
    if profile == "transition"
        case == "t120_rho035" || throw(ArgumentError("transition profile uses T=120 only"))
        return M.MultiBranchConfig(repo_root=Main.REPO_ROOT, output_dir=dirname(output), run_id=case,
            T_MeV=120.0, rho_min=0.6, rho_max=3.0, rho_anchor=1.8, rho_step=0.1,
            production_source_rho=3.2, production_source_step=0.1,
            p_num=24, t_num=8, ds_rho=0.02, dsmax_rho=0.1, max_steps=100,
            root_distance_tol=ROOT_DISTANCE, run_phase3_review=false)
    end
    first_case = case == "t120_rho035"
    anchor = first_case ? 0.35 : 0.80
    width = profile == "local" ? 0.05 : 0.10
    return M.MultiBranchConfig(repo_root=Main.REPO_ROOT, output_dir=dirname(output), run_id=case,
        T_MeV=first_case ? 120.0 : 130.0,
        rho_min=anchor-width, rho_max=anchor+width, rho_anchor=anchor,
        rho_step=profile == "local" ? 0.025 : 0.01,
        p_num=profile == "higher" ? 24 : 8, t_num=profile == "higher" ? 8 : 4,
        root_distance_tol=ROOT_DISTANCE, run_phase3_review=false)
end

solve_kwargs(cfg, work; policy=:single_pass) = (; xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num,
    iterations=300, residual_norm_max=RESIDUAL_TOL, work_telemetry=work.solver,
    fixedasymrho_seed_policy=policy)
function rho_grid(cfg; dense=false)
    intervals = round(Int, (cfg.rho_max-cfg.rho_min)/cfg.rho_step)
    return reverse(round.(collect(range(cfg.rho_min, cfg.rho_max;
        length=(dense ? 2 : 1)*intervals+1)); digits=10))
end
mode_at(cfg, rho) = Models.FixedAsymmetricRho(rho, cfg.asym_ud_ratio_target, cfg.asym_s_target)

function failed_point(rho; reason="solver_failed")
    return (; rho, solution=fill(NaN, 8), pressure=NaN, rho_norm=NaN,
        residual_norm=Inf, accepted=false, reason)
end

function checked_point(model, cfg, rho, solution, solver_ok, work)
    x = Float64.(solution)
    T = cfg.T_MeV / U.HBARC_MEV_FM
    F = zeros(8)
    residual! = Models.build_residual!(model, mode_at(cfg, rho), T;
        xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num)
    work.diagnostic_residual_calls += 1
    residual!(F, x)
    residual_norm = norm(F, Inf)
    work.thermo_evaluations += 1
    thermo = Models.compute_thermo_from_solution(model, x, T;
        xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num, rho0_scale=Models.ρ0)
    physical = Models._thermo_quantities_finite(thermo) &&
        Models._is_mass_positive(thermo.masses) && Models._is_phi_in_range(thermo.x_state) &&
        Models._is_mass_hierarchy_physical(thermo.masses)
    accepted = solver_ok && physical && isfinite(residual_norm) && residual_norm <= RESIDUAL_TOL &&
        abs(thermo.rho_norm - rho) <= RESIDUAL_TOL
    return (; rho, solution=x, pressure=Float64(thermo.pressure),
        mu_inv_fm=Float64.(thermo.mu_vec), density_inv_fm3=Float64.(thermo.rho_vec),
        helmholtz_inv_fm4=Float64(-thermo.pressure + dot(thermo.mu_vec, thermo.rho_vec)),
        masses_inv_fm=Float64.(thermo.masses),
        rho_norm=Float64(thermo.rho_norm), residual_norm, accepted,
        reason=accepted ? "accepted" : !Models._is_mass_hierarchy_physical(thermo.masses) ?
            "mass_hierarchy_invalid" : "common_residual_or_physicality_gate")
end

function direct_root(model, cfg, rho, seed, work)
    work.direct_solves += 1
    T = cfg.T_MeV / U.HBARC_MEV_FM
    raw! = Models.build_residual!(model, mode_at(cfg, rho), T;
        xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num)
    counted! = function (F, x)
        work.direct_residual_evaluations += 1
        raw!(F, x)
    end
    try
        result = NLsolve.nlsolve(counted!, copy(seed); autodiff=:forward,
            method=:trust_region, xtol=ROOT_TOL, ftol=ROOT_TOL, iterations=300)
        work.direct_nlsolve_f_calls += result.f_calls
        work.direct_nlsolve_g_calls += result.g_calls
        work.direct_iterations += result.iterations
        point = checked_point(model, cfg, rho, result.zero, result.f_converged, work)
        work.direct_nonconverged += !point.accepted
        return point
    catch err
        err isa InterruptException && rethrow()
        work.direct_exceptions += 1
        work.direct_nonconverged += 1
        return failed_point(rho; reason=sprint(showerror, err))
    end
end

function add_distinct!(roots, point)
    point.accepted || return
    any(r -> norm(r.solution - point.solution) <= ROOT_DISTANCE, roots) || push!(roots, point)
end

function discover_anchors(model, cfg, work)
    T, rho = cfg.T_MeV / U.HBARC_MEV_FM, cfg.rho_anchor
    mode = mode_at(cfg, rho)
    roots = NamedTuple[]
    catalog = Models.get_all_seeds(Models.MultiSeed(), [T], mode)
    work.anchor_requests += 1
    result = Models.solve_multi(model, mode, T; seeds=catalog,
        semantic_mode=:ground_state, evaluate_all_attempts=true, solve_kwargs(cfg, work)...)
    add_distinct!(roots, checked_point(model, cfg, rho, result.solution, result.converged, work))

    # Single-seed solves retain candidate roots which an inner residual-ranked
    # production request may discard before the outer pressure selector sees them.
    for seed in catalog
        work.anchor_requests += 1
        add_distinct!(roots, direct_root(model, cfg, rho, seed, work))
    end
    seed = nothing
    for source_rho in M._rho_path(cfg.production_source_rho, rho, cfg.production_source_step)
        work.anchor_requests += 1
        result = seed === nothing ?
            Models.solve(model, mode_at(cfg, source_rho), T; seed_strategy=Models.MultiSeed(),
                evaluate_all_attempts=true, solve_kwargs(cfg, work)...) :
            Models.solve(model, mode_at(cfg, source_rho), T; seed_guess=seed,
                evaluate_all_attempts=true, solve_kwargs(cfg, work)...)
        seed = copy(result.solution)
    end
    add_distinct!(roots, checked_point(model, cfg, rho, result.solution, result.converged, work))
    sort!(roots; by=r -> -r.pressure)
    work.anchor_count = length(roots)
    return roots
end

function palc_branch(model, cfg, root, work)
    raw = S._residual_fn(model, M._effective_cfg(cfg))
    residual = function (x, p)
        work.bk_residual_evaluations += 1
        raw(x, p)
    end
    jacobian = function (x, p)
        work.bk_jacobian_calls += 1
        U._finite_difference_jacobian(residual, x, p)
    end
    record = (x, rho; kwargs...) -> (; solution=copy(x), rho=Float64(rho))
    params = (; T_fm=cfg.T_MeV / U.HBARC_MEV_FM, rho_target=cfg.rho_anchor,
        asym_ud_ratio_target=cfg.asym_ud_ratio_target, asym_s_target=cfg.asym_s_target)
    problem = BifurcationKit.BifurcationProblem(residual, copy(root.solution), params,
        (@optic _.rho_target); J=jacobian, record_from_solution=record)
    options = BifurcationKit.ContinuationPar(p_min=cfg.rho_min, p_max=cfg.rho_max,
        ds=cfg.ds_rho, dsmax=cfg.dsmax_rho, max_steps=cfg.max_steps, detect_bifurcation=0,
        newton_options=BifurcationKit.NewtonPar(tol=ROOT_TOL, max_iterations=10))
    callback = (state; kwargs...) -> Support.observe_newton!(work.bk, state)
    try
        branch = BifurcationKit.continuation(problem, BifurcationKit.PALC(), options;
            bothside=true, verbosity=0, callback_newton=callback)
        points = [(rho=Float64(p.param), solution=Float64.(p.solution)) for p in branch.branch]
        work.branch_count += 1
        work.continuation_points += length(points)
        return points
    finally
        Support.finish_attempt!(work.bk)
    end
end

function run_palc(model, cfg, work)
    roots = NamedTuple[]
    anchor_seconds = @elapsed roots = discover_anchors(model, cfg, work)
    branches = Vector{NamedTuple}[]
    errors = String[]
    continuation_seconds = @elapsed for root in roots
        try
            push!(branches, palc_branch(model, cfg, root, work))
        catch err
            err isa InterruptException && rethrow()
            work.bk_exceptions += 1
            push!(errors, sprint(showerror, err))
        end
    end
    points = NamedTuple[]
    candidates = NamedTuple[]
    correction_seconds = @elapsed for rho in rho_grid(cfg)
        local_roots = NamedTuple[]
        for (branch_id, branch) in enumerate(branches), seed in Support.target_seeds(branch, rho)
            work.correction_requests += 1
            point = direct_root(model, cfg, rho, seed, work)
            push!(candidates, (; branch_id, point...))
            add_distinct!(local_roots, point)
        end
        if isempty(local_roots)
            push!(points, failed_point(rho; reason="no_accepted_corrected_branch"))
        else
            push!(points, local_roots[argmax(p.pressure for p in local_roots)])
        end
    end
    return (; points, anchors=roots, branches, candidates, errors,
        stages=(; anchor_seconds, continuation_seconds, correction_seconds))
end

function run_baseline(method, model, cfg, work)
    T = cfg.T_MeV / U.HBARC_MEV_FM
    grid = rho_grid(cfg; dense=method == "dense")
    method == "seed_reverse" && reverse!(grid)
    points = NamedTuple[]
    if method in ("seed", "dense", "legacy_seed", "seed_reverse")
        path = Models.FixedAsymmetricRhoPath(T, grid; ud_ratio_target=cfg.asym_ud_ratio_target,
            s_target=cfg.asym_s_target, xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num)
        result = Models.solve_path(model, path; continuation_strategy=Models.SeedContinuation(),
            branch_policy=Models.GroundStateBranchPolicy(cfg.pressure_gap_tol, RESIDUAL_TOL),
            iterations=300, residual_norm_max=RESIDUAL_TOL, work_telemetry=work.solver,
            fixedasymrho_seed_policy=method == "legacy_seed" ? :nested : :single_pass)
        for p in only(result.branches).points
            push!(points, checked_point(model, cfg, p.param_value, p.solution, p.converged, work))
        end
    else
        previous = nothing
        for rho in grid
            mode = mode_at(cfg, rho)
            result = if method == "production"
                seeds = _trho_asymmetric_seed_candidates(T, mode, previous)
                _solve_trho_asymmetric_equilibrium(model, mode, T; seeds, solve_kwargs(cfg, work)...)
            else
                Models.solve_multi(model, mode, T; semantic_mode=:ground_state,
                    evaluate_all_attempts=true, solve_kwargs(cfg, work)...)
            end
            push!(points, checked_point(model, cfg, rho, result.solution, result.converged, work))
            previous = copy(result.solution)
        end
    end
    return (; points, errors=String[])
end

function work_snapshot(work)
    counters = Dict{String,Any}(String(k) => getfield(work, k) for k in fieldnames(PilotWork) if k ∉ (:solver, :bk))
    counters["solver"] = Models.solver_work_snapshot(work.solver)
    counters["bk_newton"] = (; attempts=work.bk.attempts, iterations=work.bk.iterations,
        nonconverged=work.bk.nonconverged)
    return counters
end

function run(method, case, output)
    method in METHODS || throw(ArgumentError("unknown method: $method"))
    cfg = config(case, output)
    model = Models.create_model(:PNJL)
    repetitions = Any[]
    repetition_count = parse(Int, get(ENV, "PALC_PILOT_REPETITIONS", "6"))
    repetition_count >= 3 || throw(ArgumentError("PALC_PILOT_REPETITIONS must be at least 3"))
    # Reset all seeds and anchor discovery every time; only Julia compilation
    # and model quadrature caches survive between identical repetitions.
    for repetition in 1:repetition_count
        work = PilotWork()
        GC.gc()
        timed = @timed method == "palc" ? run_palc(model, cfg, work) : run_baseline(method, model, cfg, work)
        push!(repetitions, (; repetition, phase=repetition == 1 ? "first_call" : "warm",
            wall_seconds=timed.time, gc_seconds=timed.gctime, allocated_bytes=timed.bytes,
            compile_seconds=hasproperty(timed, :compile_time) ? timed.compile_time : nothing,
            work=work_snapshot(work), result=timed.value))
        U._write_json(output, (; schema="palc_matched_pilot_v2", method, case,
            source_sha=readchomp(`git -C $(cfg.repo_root) rev-parse HEAD`),
            julia_version=string(VERSION), bifurcationkit_version=string(Base.pkgversion(BifurcationKit)),
            threads=Threads.nthreads(), blas_threads=BLAS.get_num_threads(),
            load_seconds=Main.MATCHED_PILOT_LOAD_SECONDS,
            expected_repetitions=repetition_count,
            config=(; T_MeV=cfg.T_MeV, targets=rho_grid(cfg), xi=cfg.xi, p_num=cfg.p_num,
                profile=get(ENV, "PALC_PILOT_PROFILE", "local"),
                t_num=cfg.t_num, ud_ratio=cfg.asym_ud_ratio_target, s_target=cfg.asym_s_target,
                residual_tol=RESIDUAL_TOL, root_tol=ROOT_TOL, root_distance_tol=ROOT_DISTANCE,
                ds=cfg.ds_rho, dsmax=cfg.dsmax_rho, max_steps=cfg.max_steps,
                anchor_source_rho=cfg.production_source_rho, anchor_source_step=cfg.production_source_step),
            repetitions))
        println("$(case) $(method) $(repetition): $(round(timed.time; digits=3)) s; accepted=$(count(p -> p.accepted, timed.value.points))/$(length(timed.value.points))")
        flush(stdout)
    end
end

end

const MATCHED_PILOT_LOAD_SECONDS = time() - MATCHED_PILOT_LOAD_STARTED
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 || throw(ArgumentError("usage: run_matched_pilot.jl METHOD CASE OUTPUT_JSON"))
    MatchedFixedAsymRhoPilot.run(ARGS...)
end
