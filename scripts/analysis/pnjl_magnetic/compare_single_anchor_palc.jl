#!/usr/bin/env julia

"""
    compare_single_anchor_palc.jl

Small, isolated validation probe for the magnetic PNJL entrypoint.  It is not
a production scan and does not modify the magnetic solver.  A single accepted
root is continued with BifurcationKit/PALC at fixed T and eB; roots obtained
from PALC target crossings are corrected by the magnetic residual and compared
with the independent multi-seed `solve_magnetic_gap` catalogue.

Run in the root Julia environment.  `current_mfir` uses the production
regularization; `paper_smooth` uses the fifth chapter's N=10 Landau vacuum.
Both use the upstream IMC coefficient a=0.0108805.  The printed thesis-table
variant a=0.108805 remains explicit as `table_smooth`, rather than silently
changing the production model.  Results are local validation evidence.
"""

using Dates
using JSON
using LinearAlgebra: norm
using Printf
using StaticArrays
using ForwardDiff
import BifurcationKit
import BifurcationKit: @optic

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
using .Models

const HBARC = 197.3269804

Base.@kwdef struct ProbeConfig
    T_MeV::Float64 = 5.0
    muB_min_MeV::Float64 = 500.0
    muB_max_MeV::Float64 = 1700.0
    target_muB_MeV::Vector{Float64} = collect(900.0:50.0:1650.0)
    eB_GeV2::Vector{Float64} = [0.05, 0.155, 0.40]
    p_num::Int = 32
    t_num::Int = 4
    pz_max::Float64 = 40.0
    n_max::Int = 32
    ds::Float64 = 0.004
    dsmax::Float64 = 0.018
    max_steps::Int = 900
    newton_tol::Float64 = 2e-8
    newton_max_iterations::Int = 16
    residual_norm_max::Float64 = 2e-6
    root_merge_tol::Float64 = 2e-5
    anchor_side::Symbol = :high
    jacobian_method::Symbol = :forward
    linesearch::Bool = false
    compare_catalog::Bool = true
end

@inline _muq_fm(muB_MeV) = Float64(muB_MeV) / (3 * HBARC)
@inline _T_fm(T_MeV) = Float64(T_MeV) / HBARC
@inline _eB_fm2(eB_GeV2) = Float64(eB_GeV2) * (1000 / HBARC)^2

# Failure evidence must remain valid JSON and must not truncate a previous
# complete case. A missing crossing distance is null, never a physical zero.
_json_value(x::AbstractFloat) = isfinite(x) ? x : nothing
_json_value(x::AbstractArray) = map(_json_value, x)
_json_value(x::NamedTuple) = Dict(String(k) => _json_value(v) for (k,v) in pairs(x))
_json_value(x::AbstractDict) = Dict(k => _json_value(v) for (k,v) in pairs(x))
_json_value(x) = x

function _write_result(path, result)
    temporary=path*".pending"
    open(temporary,"w") do io
        JSON.print(io,_json_value(result),2)
        println(io)
    end
    mv(temporary,path;force=true)
end

function _convention(name::Symbol)
    name in (:paper_smooth, :table_smooth) && return (
        route=:landau_legacy,
        imc=Models.MagneticThermodynamics.MagneticIMCParams(
            a=name === :table_smooth ? 0.108805 : 0.0108805,
            b=-1.0133e-4, c=0.02228, d=1.84558e-4,
            Λ_QCD_MeV=300.0,
        ),
    )
    name === :current_mfir && return (
        route=:mfir,
        imc=Models.MagneticThermodynamics.MagneticIMCParams(
            a=0.0108805, b=-1.0133e-4, c=0.02228, d=1.84558e-4,
            Λ_QCD_MeV=300.0,
        ),
    )
    throw(ArgumentError("unknown convention $(name)"))
end

function _model(eB_GeV2::Real, convention::Symbol)
    c = _convention(convention)
    return Models.PNJLMagneticModel(
        eB_fm2=_eB_fm2(eB_GeV2),
        route=c.route,
        imc=c.imc,
    )
end

function _solve_one(model, T_fm, muB_MeV, seed, cfg::ProbeConfig)
    muq = _muq_fm(muB_MeV)
    return Models.solve_magnetic_gap(
        model, T_fm, Models.normalize_mu_vec(muq);
        p_num=cfg.p_num, t_num=cfg.t_num, pz_max=cfg.pz_max,
        n_max=cfg.n_max, initial_guess=seed, include_default_seeds=false,
        method=:trust_region, fallback_method=:newton, iterations=260,
        residual_norm_max=cfg.residual_norm_max, root_merge_tol=cfg.root_merge_tol,
    )
end

function _anchor(model, cfg::ProbeConfig, anchor_muB_MeV::Float64)
    # The low-density side in the reference curves is the hadronic root.  A
    # single seed is used here; no catalog is consulted during anchor solve.
    seed = cfg.anchor_side === :high ? Models.QUARK_SEED_5 : Models.HADRON_SEED_5
    result = _solve_one(model, _T_fm(cfg.T_MeV), anchor_muB_MeV, seed, cfg)
    result.converged || error("single-anchor solve failed at muB=$(anchor_muB_MeV)")
    return result.candidates[result.selected_index].x_state
end

function _state_record(model, x, muB_MeV, cfg)
    state = SVector{5,Float64}(x)
    mu = Models.normalize_mu_vec(_muq_fm(muB_MeV))
    controls = (;p_num=cfg.p_num, t_num=cfg.t_num, pz_max=cfg.pz_max, n_max=cfg.n_max)
    comp = Models.calculate_magnetic_omega_components(state,mu,_T_fm(cfg.T_MeV),model.magnetic;controls...)
    densities = Models.calculate_magnetic_rho(state,mu,_T_fm(cfg.T_MeV),model.magnetic;controls...)
    residual = Models.magnetic_gap_residual(model,state,_T_fm(cfg.T_MeV),mu;controls...)
    return (; x_state=collect(state), omega=comp.omega, residual_norm=norm(residual),
        masses_MeV=collect(comp.masses .* HBARC), density_inv_fm3=collect(densities),
        rhoB_over_rho0=sum(densities)/(3*model.base.params.rho0_fm3),
        physical=all(>(0),comp.masses) && all(v->0<=v<=1,state[4:5]))
end

function _residual(model, cfg::ProbeConfig)
    T_fm = _T_fm(cfg.T_MeV)
    p_num, t_num, pz_max, n_max = cfg.p_num, cfg.t_num, cfg.pz_max, cfg.n_max
    return function (x, p)
        mu = Models.normalize_mu_vec(p.muq_fm)
        Vector(Models.magnetic_gap_residual(
            model, x, T_fm, mu;
            p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max,
        ))
    end
end

function _jacobian(residual, x, p)
    x0 = Float64.(x)
    f0 = residual(x0, p)
    J = Matrix{Float64}(undef, length(f0), length(x0))
    for j in eachindex(x0)
        h = cbrt(eps(Float64)) * max(1.0, abs(x0[j]))
        xp, xm = copy(x0), copy(x0)
        xp[j] += h
        xm[j] -= h
        @views J[:, j] .= (residual(xp, p) .- residual(xm, p)) ./ (2h)
    end
    return J
end

function _palc(model, cfg::ProbeConfig, anchor, anchor_muB_MeV::Float64)
    residual = _residual(model, cfg)
    p0 = (; muq_fm=_muq_fm(anchor_muB_MeV))
    problem = BifurcationKit.BifurcationProblem(
        residual, collect(Float64.(anchor)), p0, (@optic _.muq_fm);
        J=(x, p) -> cfg.jacobian_method===:forward ?
            ForwardDiff.jacobian(y->residual(y,p),Float64.(x)) : _jacobian(residual,x,p),
        record_from_solution=(x, p; kwargs...) -> (; solution=Float64.(x)),
    )
    opts = BifurcationKit.ContinuationPar(
        p_min=_muq_fm(cfg.muB_min_MeV),
        p_max=_muq_fm(cfg.muB_max_MeV),
        ds=cfg.ds, dsmax=cfg.dsmax, max_steps=cfg.max_steps,
        detect_bifurcation=0,
        newton_options=BifurcationKit.NewtonPar(
            tol=cfg.newton_tol, max_iterations=cfg.newton_max_iterations,
            verbose=false, linesearch=cfg.linesearch,
        ),
    )
    return BifurcationKit.continuation(problem, BifurcationKit.PALC(), opts;
        bothside=true, verbosity=0)
end

function _crossing_seeds(branch, target_muq)
    points = branch.branch
    seeds = Vector{Vector{Float64}}()
    for point in points
        Float64(point.param) == target_muq && push!(seeds, Float64.(point.solution))
    end
    for i in 2:length(points)
        p0, p1 = Float64(points[i-1].param), Float64(points[i].param)
        p0 == p1 && continue
        min(p0, p1) < target_muq < max(p0, p1) || continue
        w = (target_muq - p0) / (p1 - p0)
        x0, x1 = Float64.(points[i-1].solution), Float64.(points[i].solution)
        push!(seeds, (1 - w) .* x0 .+ w .* x1)
    end
    unique_seeds = Vector{Vector{Float64}}()
    for seed in seeds
        any(norm(seed - other) <= 1e-5 for other in unique_seeds) || push!(unique_seeds, seed)
    end
    return unique_seeds
end

function _correct_crossings(model, branch, muB_MeV, cfg)
    corrected = NamedTuple[]
    for seed in _crossing_seeds(branch, _muq_fm(muB_MeV))
        result = try
            _solve_one(model, _T_fm(cfg.T_MeV), muB_MeV, seed, cfg)
        catch
            nothing
        end
        result === nothing && continue
        result.converged || continue
        for candidate in result.candidates
            any(norm(candidate.x_state - item.x_state) <= cfg.root_merge_tol for item in corrected) ||
                push!(corrected, _state_record(model,candidate.x_state,muB_MeV,cfg))
        end
    end
    return corrected
end

function _catalog(model, muB_MeV, cfg)
    result = Models.solve_magnetic_gap(
        model, _T_fm(cfg.T_MeV), Models.normalize_mu_vec(_muq_fm(muB_MeV));
        p_num=cfg.p_num, t_num=cfg.t_num, pz_max=cfg.pz_max, n_max=cfg.n_max,
        include_default_seeds=true, method=:trust_region, fallback_method=:newton,
        iterations=260, residual_norm_max=cfg.residual_norm_max,
        root_merge_tol=cfg.root_merge_tol,
    )
    return [candidate for candidate in result.candidates if candidate.converged]
end

function _matches(palc_roots, catalog_roots, tol)
    isempty(catalog_roots) && return 0
    return count(catalog_roots) do c
        any(norm(c.x_state - p.x_state) <= tol for p in palc_roots)
    end
end

function _branch_summary(branch, cfg)
    params = Float64[point.param * 3 * HBARC for point in branch.branch]
    isempty(params) && return (; points=0, muB_min_MeV=NaN, muB_max_MeV=NaN, folds=0)
    folds = 0
    for i in 2:(length(params)-1)
        (params[i] - params[i-1]) * (params[i+1] - params[i]) < 0 && (folds += 1)
    end
    return (; points=length(params), muB_min_MeV=minimum(params),
        muB_max_MeV=maximum(params), folds=folds)
end

function run_case(cfg::ProbeConfig, convention::Symbol, eB_GeV2::Float64)
    model = _model(eB_GeV2, convention)
    anchor_muB_MeV = cfg.anchor_side === :high ? cfg.muB_max_MeV :
        (eB_GeV2 >= 0.30 ? 500.0 : 850.0)
    started = time()
    anchor = _anchor(model, cfg, anchor_muB_MeV)
    branch = _palc(model, cfg, anchor, anchor_muB_MeV)
    palc_seconds = time()-started
    println("  continued ",length(branch.branch)," points in ",palc_seconds," s"); flush(stdout)
    summary = _branch_summary(branch, cfg)
    points = [merge((;index=i,muB_MeV=Float64(p.param)*3*HBARC),
        _state_record(model,p.solution,Float64(p.param)*3*HBARC,cfg))
        for (i,p) in enumerate(branch.branch)]
    targets = NamedTuple[]
    for muB in cfg.target_muB_MeV
        palc_roots = _correct_crossings(model, branch, muB, cfg)
        catalog_error = nothing
        catalog = try
            cfg.compare_catalog ? _catalog(model,muB,cfg) : Models.MagneticGapCandidate[]
        catch err
            catalog_error = sprint(showerror,err)
            Models.MagneticGapCandidate[]
        end
        catalog_records = [_state_record(model,c.x_state,muB,cfg) for c in catalog]
        distances = [isempty(palc_roots) ? Inf : minimum(norm(c.x_state-p.x_state) for p in palc_roots)
            for c in catalog]
        push!(targets, (; muB_MeV=muB, palc_root_count=length(palc_roots),
            catalog_root_count=length(catalog), matched_catalog_roots=
            _matches(palc_roots, catalog, 5e-4),
            catalog_error, catalog_min_state_distances=distances,
            palc_roots, catalog_roots=catalog_records,
            catalog_residual_max=isempty(catalog) ? NaN : maximum(c.residual_norm for c in catalog),
            palc_residual_max=isempty(palc_roots) ? NaN : maximum(p.residual_norm for p in palc_roots)))
        println("  target ",muB," roots PALC/catalog=",length(palc_roots),"/",length(catalog)); flush(stdout)
    end
    return (; convention=String(convention), eB_GeV2, T_MeV=cfg.T_MeV,
        route=String(model.magnetic.route), imc_a=model.magnetic.imc.a,
        n_max=cfg.n_max, p_num=cfg.p_num, pz_max=cfg.pz_max,
        anchor_muB_MeV, anchor_side=String(cfg.anchor_side), anchor_state=collect(anchor),
        jacobian_method=String(cfg.jacobian_method),
        ds=cfg.ds, dsmax=cfg.dsmax, max_steps=cfg.max_steps,
        newton_tol=cfg.newton_tol, newton_max_iterations=cfg.newton_max_iterations,
        linesearch=cfg.linesearch, compare_catalog=cfg.compare_catalog,
        palc_seconds,total_seconds=time()-started, branch=summary, points, targets)
end

function main(args=ARGS)
    options = Dict{String,String}()
    output = joinpath(PROJECT_ROOT,"tmp","magnetic_palc_validation.json")
    for arg in args
        if startswith(arg,"--")
            parts=split(arg[3:end],"=";limit=2)
            length(parts)==2 || throw(ArgumentError("expected --key=value"))
            options[String(parts[1])]=String(parts[2])
        else
            output=arg
        end
    end
    allowed=("T","p-num","pz-max","n-max","ds","dsmax","max-steps","anchor-side","convention","eB","targets","jacobian","linesearch","catalog","newton-tol","newton-iterations")
    all(key->key in allowed,keys(options)) || throw(ArgumentError("unknown probe option"))
    cfg = ProbeConfig(T_MeV=parse(Float64,get(options,"T","5")),
        p_num=parse(Int,get(options,"p-num","32")),pz_max=parse(Float64,get(options,"pz-max","40")),
        n_max=parse(Int,get(options,"n-max","32")),ds=parse(Float64,get(options,"ds","0.004")),
        dsmax=parse(Float64,get(options,"dsmax","0.018")),
        max_steps=parse(Int,get(options,"max-steps","900")),
        anchor_side=Symbol(get(options,"anchor-side","high")),
        jacobian_method=Symbol(get(options,"jacobian","forward")),
        linesearch=parse(Bool,get(options,"linesearch","false")),
        compare_catalog=parse(Bool,get(options,"catalog","true")),
        newton_tol=parse(Float64,get(options,"newton-tol","2e-8")),
        newton_max_iterations=parse(Int,get(options,"newton-iterations","16")),
        target_muB_MeV=parse.(Float64,split(get(options,"targets","900,1000,1050,1100,1150,1200,1500,1650"),',')))
    cfg.anchor_side in (:low,:high) || throw(ArgumentError("anchor-side must be low or high"))
    cfg.jacobian_method in (:forward,:finite) || throw(ArgumentError("jacobian must be forward or finite"))
    cases = Any[]
    conventions=Symbol.(split(get(options,"convention","current_mfir"),','))
    fields=parse.(Float64,split(get(options,"eB","0.05,0.155,0.4"),','))
    mkpath(dirname(output))
    for convention in conventions
        for eB in fields
            println("running convention=$(convention) eB=$(eB) GeV^2")
            flush(stdout)
            result = try
                run_case(cfg,convention,eB)
            catch err
                err isa InterruptException && rethrow()
                (;convention=String(convention),eB_GeV2=eB,T_MeV=cfg.T_MeV,
                    error=sprint(showerror,err))
            end
            push!(cases,result)
            _write_result(output,(;generated_at=string(Dates.now()),julia_version=string(VERSION),cases))
        end
    end
    result = (; generated_at=string(Dates.now()), julia_version=string(VERSION), cases)
    mkpath(dirname(output))
    _write_result(output,result)
    println("wrote $(output)")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
