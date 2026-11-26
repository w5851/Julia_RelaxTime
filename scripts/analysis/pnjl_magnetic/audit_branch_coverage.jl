#!/usr/bin/env julia

"""Match corrected production roots to independent MFIR roots.

Also examine the narrow extra p512 folds near muB=872.08 MeV and any
unmatched production catalogue root.  Every corrective solve uses just its
saved branch state as an initial guess.  No production entrypoint is changed.
"""

using NLsolve
include(joinpath(@__DIR__, "audit_palc_production_gap.jl"))

function select_case(data, field)
    return only(c for c in data["cases"] if c["eB_GeV2"] == field)
end

function state_distance(a, b)
    return norm(Float64.(a["x_state"]) - Float64.(b["x_state"]))
end

function root_matching(production, reference)
    forward = [minimum(state_distance(root, other) for other in reference)
        for root in production]
    reverse = [minimum(state_distance(root, other) for other in production)
        for root in reference]
    separation = length(production) < 2 ? nothing : minimum(
        state_distance(production[i], production[j])
        for i in 1:length(production)-1 for j in i+1:length(production))
    return (; production_count=length(production), reference_count=length(reference),
        production_to_reference=forward, reference_to_production=reverse,
        max_state_distance=max(maximum(forward), maximum(reverse)),
        minimum_production_root_separation=separation)
end

function interpolate_crossings(points, target)
    seeds = Vector{Float64}[]
    for point in points
        abs(point["muB_MeV"] - target) < 1e-10 &&
            push!(seeds, Float64.(point["x_state"]))
    end
    for i in 2:length(points)
        a, b = points[i-1], points[i]
        left, right = a["muB_MeV"], b["muB_MeV"]
        min(left, right) < target < max(left, right) || continue
        w = (target-left)/(right-left)
        push!(seeds, (1-w).*Float64.(a["x_state"]) .+ w.*Float64.(b["x_state"]))
    end
    return seeds
end

function correct_production(seed, c, target, order)
    eB = c["eB_GeV2"]*(1000/HBARC)^2
    model = Models.PNJLMagneticModel(eB_fm2=eB, p_num=order, n_max=32, pz_max=40.0)
    mu = Models.normalize_mu_vec(target/(3*HBARC))
    result = Models.solve_magnetic_gap(model, c["T_MeV"]/HBARC, mu;
        p_num=order, n_max=32, pz_max=40.0, t_num=4,
        initial_guess=seed, include_default_seeds=false,
        method=:trust_region, fallback_method=:newton, iterations=260,
        ftol=1e-11, xtol=1e-12, residual_norm_max=2e-7, root_merge_tol=2e-5)
    result.converged || return (; converged=false, x_state=nothing, audit=nothing)
    root = result.candidates[result.selected_index]
    record = Dict("x_state"=>collect(root.x_state))
    return (; converged=true, x_state=record["x_state"],
        audit=audit_state(c, target, record, order))
end

function correct_reference(seed, c, target)
    eB = c["eB_GeV2"]*(1000/HBARC)^2
    T, mu = c["T_MeV"]/HBARC, fill(target/(3*HBARC), 3)
    params, ctl = MagneticReference.Parameters(), MagneticReference.Controls(order=48)
    function residual(x)
        try
            return MagneticReference.evaluate(x, T, mu, eB, params, ctl).residual
        catch err
            err isa DomainError || rethrow()
            return fill(1e6, 5)
        end
    end
    residual!(r, x) = (r .= residual(x))
    function jacobian!(J, x)
        for j in eachindex(x)
            h = 2e-5*max(1.0, abs(x[j]))
            xp, xm = copy(x), copy(x)
            xp[j] += h
            xm[j] -= h
            J[:, j] .= (residual(xp) .- residual(xm)) ./ (2h)
        end
    end
    try
        solved = nlsolve(residual!, jacobian!, copy(seed);
            method=:trust_region, ftol=1e-11, xtol=1e-12, iterations=180)
        ev = MagneticReference.evaluate(solved.zero, T, mu, eB, params, ctl)
        return (; converged=norm(ev.residual)<1e-7,
            x_state=collect(solved.zero), residual_norm=norm(ev.residual),
            mass_jacobian_sigma_min=ev.mass_jacobian_sigma_min)
    catch err
        err isa InterruptException && rethrow()
        return (; converged=false, x_state=nothing, error=sprint(showerror, err))
    end
end

function count_unique(records)
    states = Vector{Float64}[]
    for record in records
        record.converged || continue
        any(norm(record.x_state - other)<2e-5 for other in states) ||
            push!(states, record.x_state)
    end
    return length(states)
end

function run_coverage_audit(args=ARGS)
    length(args) == 4 || throw(ArgumentError(
        "production branch JSON, independent branch JSON, corrected roots JSON, output JSON required"))
    production_path, reference_path, corrected_path, output = args
    branch = only(JSON.parsefile(production_path)["cases"])
    reference = select_case(JSON.parsefile(reference_path), branch["eB_GeV2"])
    corrected = select_case(JSON.parsefile(corrected_path), branch["eB_GeV2"])

    targets = Any[]
    for target in corrected["targets"]
        roots = filter(r->r["converged"], target["roots"])
        ref = only(t for t in reference["targets"] if t["muB_MeV"]==target["muB_MeV"])
        push!(targets, (; muB_MeV=target["muB_MeV"],
            root_matching(roots, ref["palc_roots"])...))
    end

    catalogue = Any[]
    for target in branch["targets"], (i, root) in enumerate(target["catalog_roots"])
        distance = minimum(state_distance(root, r) for r in target["palc_roots"])
        distance > 5e-4 || continue
        push!(catalogue, (; muB_MeV=target["muB_MeV"], root_index=i, distance,
            production_audit=audit_state(branch, target["muB_MeV"], root, branch["p_num"])))
    end

    # These targets straddle the two extra p512 turning points.  Crossing
    # counts and correction outcomes are diagnostic, not a global branch proof.
    narrow = Any[]
    for target in (branch["eB_GeV2"] == 0.155 ? (872.08,) : ())
        seeds = interpolate_crossings(branch["points"], target)
        ref_seeds = interpolate_crossings(reference["points"], target)
        println("narrow muB=",target," crossings production/reference=",length(seeds),"/",length(ref_seeds))
        flush(stdout)
        reference_roots = [correct_reference(seed, branch, target) for seed in ref_seeds]
        checks = Any[]
        for order in (512, 1024)
            roots = [correct_production(seed, branch, target, order) for seed in seeds]
            push!(checks, (; p_num=order, converged=count(r->r.converged, roots),
                unique_roots=count_unique(roots), roots))
        end
        corrected_reference = [correct_reference(seed, branch, target) for seed in seeds]
        push!(narrow, (; muB_MeV=target, production_crossings=length(seeds),
            reference_crossings=length(ref_seeds), reference_unique_roots=count_unique(reference_roots),
            reference_roots, reference_from_production_unique=count_unique(corrected_reference),
            reference_from_production=corrected_reference, checks))
    end

    mkpath(dirname(output))
    open(output,"w") do io
        JSON.print(io, (; julia_version=string(VERSION),
            production_source=production_path, reference_source=reference_path,
            corrected_source=corrected_path, eB_GeV2=branch["eB_GeV2"],
            targets, unmatched_catalogue=catalogue, narrow_window=narrow), 2)
    end
    println("wrote ", output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_coverage_audit()
end
