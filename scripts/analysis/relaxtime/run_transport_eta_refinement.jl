#!/usr/bin/env julia
# Quote-remote only. SOURCE_ROOT INPUT_DIR OUTPUT_DIR STAGE SELECTION PROFILE
# Stages: locate, transport. PROFILE is base/fine/tighter for transport.
using StaticArrays, NLsolve
include(joinpath(@__DIR__, "run_transport_mott_probe.jl"))
include(joinpath(SOURCE_ROOT, "scripts", "relaxtime", "transport_sampling_scan_config.jl"))
include(joinpath(SOURCE_ROOT, "scripts", "relaxtime", "gap_transport_scan_phase_equilibrium.jl"))
const SP = Models.TransportSamplingProfiles
const PHASE = Main.GapTransportScanPhaseEquilibrium
const INPUT_CONFIG = JSON3.read(read(joinpath(INPUT_DIR, "config.json"), String))
const SOURCE_ROWS = collect(CSV.File(joinpath(INPUT_DIR, "source_rows.csv")))
const SPECIES_FINE = (:u, :d, :s, :ubar, :dbar, :sbar)

function original_seed(model, row)
    target = SVector{3}(Float64(row.m_u), Float64(row.m_d), Float64(row.m_s))
    residual = x -> Vector(Models.calculate_mass_vec(model, SVector{3}(x))-target)
    solved = nlsolve(residual, [-1.0,-1.0,-1.5]; xtol=1e-12, ftol=1e-12)
    maximum(abs, residual(solved.zero)) < 1e-10 || error("original condensate reconstruction failed")
    return SVector{5}(solved.zero..., Float64(row.Phi), Float64(row.Phibar))
end

function numerical_level(profile, name)
    name == "base" && return (name=name, ap=32, ac=16, p=96, c=24, limit=25.0)
    name == "A64" && return (name=name, ap=64, ac=32, p=96, c=24, limit=25.0)
    name == "background_ref" && return (name=name, ap=64, ac=32, p=128, c=32, limit=30.0)
    error("unknown background/A level")
end

function background(target, xi, profile, level; seed=nothing)
    model = SP.profile_model(profile; thermal_p_max_inv_fm=level.limit)
    opts = (model=model, p_num=level.p, t_num=level.c, residual_norm_max=1e-10, max_iter=80)
    T = Float64(target.T_MeV)/HC
    mu = Float64(target.muB_MeV)/(3HC)
    eq = PHASE.solve_with_multiseed_governance(T, mu, Float64(xi), opts; continuation_seed=seed)
    eq === nothing && error("stable multiseed background failed at xi=$xi")
    eq.residual_norm <= 1e-10 || error("insufficient background residual")
    row = (point_id=String(target.curve_id), window=String(target.window),
        T_MeV=Float64(target.T_MeV), muB_MeV=Float64(target.muB_MeV),
        xi=Float64(xi), T_fm=T, muq_fm=mu,
        m_u=eq.masses[1], m_d=eq.masses[2], m_s=eq.masses[3],
        Phi=eq.x_state[4], Phibar=eq.x_state[5],
        phi_u=eq.x_state[1], phi_d=eq.x_state[2], phi_s=eq.x_state[3],
        background_residual=eq.residual_norm, omega=eq.omega)
    state = frozen_state(row, profile, level.ap, level.ac)
    return (row=row, equilibrium=eq, state=state)
end

include(joinpath(@__DIR__, "transport_eta_threshold_locator.jl"))

function compute_transport()
    curve_id, variant = ARGS[5], ARGS[6]
    profile = SP.load_profile(joinpath(INPUT_DIR, variant*".toml"))
    runtime = SP.build_runtime(profile)
    plan_path = joinpath(get(ENV,"ETA_LOCATE_DIR",joinpath(dirname(OUTPUT_DIR),"locate_v1")), "transport_plan.csv")
    selected = filter(r->String(r.curve_id)==curve_id, collect(CSV.File(plan_path)))
    wanted = get(ENV,"ETA_POINT_IDS","all")
    evaluate_poles = get(ENV,"ETA_EVAL_POLES","false") == "true"
    if wanted != "all"
        ids=Set(split(wanted,','))
        selected=filter(r->String(r.point_id) in ids,selected)
    end
    isempty(selected) && error("no transport points for selected curve")
    target = only(t for t in INPUT_CONFIG.targets if String(t.curve_id)==curve_id)
    original = only(r for r in SOURCE_ROWS if String(r.point_id)==String(target.point_left))
    seed = original_seed(runtime.model,original)
    level = numerical_level(profile,"base")
    observables, rates, poles, candidates = NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[]
    for point in selected
        began = time()
        bg = background(target,Float64(point.xi),profile,level;seed=seed)
        seed=bg.equilibrium.x_state
        on_channel = function (r)
            push!(rates,merge((point_id=String(point.point_id),curve_id=curve_id,window=String(point.window),
                T_MeV=Float64(point.T_MeV),muB_MeV=Float64(point.muB_MeV),xi=Float64(point.xi)),r))
            CSV.write(joinpath(OUTPUT_DIR,"rates.csv"),rates)
        end
        println("POINT_START ",point.point_id," xi=",point.xi," delta=",point.delta_xi);flush(stdout)
        result=SP.solve_profile_transport(runtime,bg.equilibrium,bg.row.T_fm,bg.row.muq_fm;
            xi=bg.row.xi,compute_bulk=true,on_channel=on_channel)
        row=merge(bg.row,(point_id=String(point.point_id),curve_id=curve_id,
            xi_mott=Float64(point.xi_mott),delta_xi=Float64(point.delta_xi),
            point_kind=String(point.kind),original_point_id=ismissing(point.original_point_id) ? "" : String(point.original_point_id),
            entropy=result.thermo_background.entropy,s_fm3inv=result.thermo_background.entropy,eta=result.transport.eta,
            zeta=result.transport.zeta,sigma=result.transport.sigma,
            eta_over_s=result.transport.eta/result.thermo_background.entropy,
            zeta_over_s=result.transport.zeta/result.thermo_background.entropy,
            sigma_over_T=result.transport.sigma/bg.row.T_fm,
            seconds=time()-began),
            NamedTuple{Tuple(Symbol("tau_",s) for s in SPECIES_FINE)}(Tuple(result.tau)),
            NamedTuple{Tuple(Symbol("tauinv_",s) for s in SPECIES_FINE)}(Tuple(result.tau_inv)),
            NamedTuple{Tuple(Symbol("n_",s) for s in SPECIES_FINE)}(Tuple(result.densities)))
        all(v->isfinite(v)&&v>0,(row.eta_over_s,row.zeta_over_s,row.sigma_over_T,Tuple(result.tau)...)) ||
            error("nonfinite or nonpositive transport output")
        push!(observables,row)
        CSV.write(joinpath(OUTPUT_DIR,"observables.csv"),observables)
        if evaluate_poles
            cs,sm,_=probe_roots(row,bg.state,:eta,"point_multistart",nothing)
            append!(candidates,cs);push!(poles,sm)
            CSV.write(joinpath(OUTPUT_DIR,"poles.csv"),poles)
            CSV.write(joinpath(OUTPUT_DIR,"pole_candidates.csv"),candidates)
        end
        println("POINT_DONE ",point.point_id," tau_ubar=",row.tau_ubar," seconds=",row.seconds);flush(stdout)
        GC.gc()
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",curve_id=curve_id,
        variant=variant,points=length(observables),channels=length(rates),point_selection=wanted,
        per_point_poles=evaluate_poles,location_summary_sha256=digest(joinpath(dirname(plan_path),"summary.json")),
        background_solver=(residual_norm_max=1e-10,max_iter=80,stable_multiseed=true),
        calculation_sha=String(INPUT_CONFIG.calculation_sha),effective_profile=SP.effective_config(profile),
        input_config_sha256=digest(joinpath(INPUT_DIR,"config.json")),transport_plan_sha256=digest(plan_path),
        analysis_script_sha256=digest(@__FILE__),helper_sha256=digest(joinpath(@__DIR__,"run_transport_mott_probe.jl")),
        output_sha256=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

length(ARGS)==6 || error("SOURCE_ROOT INPUT_DIR OUTPUT_DIR STAGE SELECTION PROFILE required")
if ARGS[4]=="locate"
    locate_all()
elseif ARGS[4]=="transport"
    compute_transport()
else
    error("unknown refinement stage")
end
