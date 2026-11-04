#!/usr/bin/env julia
#= Probe meson poles on frozen transport backgrounds; run on quote-remote.

Arguments: SOURCE_ROOT INPUT_DIR OUTPUT_DIR [POINT_IDS|all] [A_P_NODES A_COS_NODES]
This diagnostic deliberately consumes m, mu, Phi and Phibar directly. The meson
workflow equilibrium adapter would rebuild A with its separate default grid.
No gap solve, transport production, or production root-policy changes occur.
Mixed results retain the algebraic minus/plus eigenmode and every root candidate;
they are not automatically assigned an accepted physical mode identity.
=#
using CSV, JSON3, SHA, Dates, LinearAlgebra

length(ARGS) in (3, 4, 6) || error("expected SOURCE_ROOT INPUT_DIR OUTPUT_DIR [POINT_IDS|all] [A_P_NODES A_COS_NODES]")
const SOURCE_ROOT = abspath(ARGS[1])
const INPUT_DIR = abspath(ARGS[2])
const OUTPUT_DIR = abspath(ARGS[3])
isdir(OUTPUT_DIR) && !isempty(readdir(OUTPUT_DIR)) && error("refusing to overwrite existing probe output")
mkpath(OUTPUT_DIR)
include(joinpath(SOURCE_ROOT, "src", "models", "Models.jl"))
const MP = Main.RelaxTime.MesonMass
const MT = Main.RelaxTime.MottTransition
const AS = Main.RelaxTime.AverageScatteringRate
const EC = Main.RelaxTime.EffectiveCouplings
const RTM = Main.RelaxTime.RelaxationTime
const HC = Main.Constants_PNJL.ħc_MeV_fm
const MESONS = (:pi, :K, :sigma_pi, :sigma_K, :eta, :eta_prime, :sigma, :sigma_prime)
const MIXED = (:eta, :eta_prime, :sigma, :sigma_prime)
const ROOT_TOL = 1e-7
LinearAlgebra.BLAS.set_num_threads(1)

digest(path) = bytes2hex(sha256(read(path)))
function write_json(path, value)
    open(path, "w") do io
        JSON3.pretty(io, value)
        println(io)
    end
end

function frozen_state(row, profile, ap, ac)
    qp0 = (m=(u=Float64(row.m_u), d=Float64(row.m_d), s=Float64(row.m_s)),
           μ=(u=Float64(row.muq_fm), d=Float64(row.muq_fm), s=Float64(row.muq_fm)))
    tp = (T=Float64(row.T_fm), Φ=Float64(row.Phi), Φbar=Float64(row.Phibar), ξ=Float64(row.xi))
    all(x -> isfinite(x) && x > 0, qp0.m) || error("invalid frozen quark masses")
    abs(tp.T * HC - row.T_MeV) < 1e-8 || error("temperature unit mismatch")
    abs(qp0.μ.u * HC * 3 - row.muB_MeV) < 1e-8 || error("chemical-potential unit mismatch")
    p, w = AS._cached_semi_infinite_momentum_grid(ap, profile.a_integral.scale_inv_fm)
    av = Main.RelaxTime.AFieldBuilder.build_A_triplet(qp0, tp, p, w; cos_nodes=ac, use_aniso=true)
    qp = merge(qp0, (A=av,))
    gu, gs = EC.calculate_G_from_A(av.u, qp.m.u), EC.calculate_G_from_A(av.s, qp.m.s)
    kc = EC.calculate_effective_couplings(Main.Constants_PNJL.G_fm2, Main.Constants_PNJL.K_fm5, gu, gs)
    return (qp=qp, tp=tp, kc=kc)
end

point_fields(r) = (point_id=String(r.point_id), window=String(r.window),
    T_MeV=Float64(r.T_MeV), muB_MeV=Float64(r.muB_MeV), xi=Float64(r.xi))

function threshold_rows(row, state)
    out = NamedTuple[]
    qp, tp, kc = state.qp, state.tp, state.kc
    for meson in MESONS
        cuts = meson in MIXED ? ((:uu, 2qp.m.u), (:ss, 2qp.m.s)) :
            ((:pair, MT.mott_threshold_mass(meson, qp)),)
        for (cut, threshold) in cuts, offset in (-1e-4, -1e-6, 0.0, 1e-6, 1e-4)
            m = threshold + offset
            f = MP.meson_mass_equation(meson, m, 0.0, 0.0, qp, tp, kc)
            det = if meson in MIXED
                ch = meson in (:eta, :eta_prime) ? :P : :S
                e = MP._mixed_matrix_elements(ch, m, 0.0, 0.0, qp, tp, kc)
                (e.M00 * e.M88 - e.M08^2) / (ch == :P ? kc.det_K_plus : kc.det_K_minus)
            else
                ComplexF64(NaN, NaN)
            end
            push!(out, merge(point_fields(row), (meson=String(meson), cut=String(cut),
                threshold_MeV=threshold*HC, offset_inv_fm=offset,
                inverse_re=real(f), inverse_im=imag(f), inverse_abs=abs(f),
                determinant_re=real(det), determinant_im=imag(det),
                inverse_units=meson in MIXED ? "fm2" : "dimensionless")))
        end
    end
    out
end

function probe_roots(row, state, meson, direction, previous)
    qp, tp, kc = state.qp, state.tp, state.kc
    threshold = MT.mott_threshold_mass(meson, qp)
    m0 = MP.default_meson_mass_guess(meson, qp)
    masses = unique([0.6m0, 0.95m0, 1.02m0, 1.4m0, 0.98*2qp.m.s, 1.03*2qp.m.s, 4.8])
    seeds = Tuple{Float64,Float64}[(m,g) for m in masses for g in (0.0, 0.1, 0.6)]
    previous === nothing || pushfirst!(seeds, previous)
    candidates = NamedTuple[]
    for (seed_index, (m,g)) in enumerate(unique(seeds))
        raw_m, raw_g, residual, conv, message = NaN, NaN, Inf, false, ""
        try
            root = MP.solve_meson_mass(meson, qp, tp; K_coeffs=kc,
                initial_mass=m, initial_gamma=g, method=:trust_region,
                ftol=1e-10, xtol=1e-11, iterations=100)
            raw_m, raw_g, conv = root.mass, root.gamma, root.converged
            # Recheck the actual positive-mass pole, including the existing sign convention.
            mm = raw_m < 0 ? -raw_m : raw_m
            gg = raw_m < 0 ? abs(raw_g) : raw_g
            f = MP.meson_mass_equation(meson, mm, gg, 0.0, qp, tp, kc)
            residual = abs(f)
        catch err
            message = sprint(showerror, err)
        end
        mm = raw_m < 0 ? -raw_m : raw_m
        gg = raw_m < 0 ? abs(raw_g) : raw_g
        qualified = conv && isfinite(mm) && mm > 0 && isfinite(gg) &&
                    gg >= -1e-8 && residual <= ROOT_TOL
        push!(candidates, merge(point_fields(row), (meson=String(meson), direction=direction,
            seed_index=seed_index, seed_mass_inv_fm=m, seed_gamma_inv_fm=g,
            raw_mass_inv_fm=raw_m, raw_gamma_inv_fm=raw_g,
            mass_MeV=mm*HC, gamma_MeV=gg*HC, residual_abs=residual,
            converged=conv, qualified=qualified, error=message)))
    end
    good = filter(c -> c.qualified, candidates)
    unique_roots = NamedTuple[]
    for candidate in sort(good; by=c -> c.residual_abs)
        any(c -> hypot(c.mass_MeV-candidate.mass_MeV,c.gamma_MeV-candidate.gamma_MeV) < 1e-3,
            unique_roots) || push!(unique_roots, candidate)
    end
    if isempty(good)
        selected = candidates[argmin([c.residual_abs for c in candidates])]
    elseif previous === nothing
        selected = good[argmin([hypot(c.mass_MeV/HC-m0,c.gamma_MeV/HC) for c in good])]
    else
        selected = good[argmin([hypot(c.mass_MeV/HC-previous[1],c.gamma_MeV/HC-previous[2]) for c in good])]
    end
    branch = meson in (:eta,:sigma) ? "minus" : meson in (:eta_prime,:sigma_prime) ? "plus" : "simple"
    summary = merge(point_fields(row), (meson=String(meson), algebraic_branch=branch,
        direction=direction, mass_MeV=selected.mass_MeV, gamma_MeV=selected.gamma_MeV,
        qualified=selected.qualified, residual_abs=selected.residual_abs,
        unique_qualified_roots=length(unique_roots),
        threshold_MeV=threshold*HC, gap_MeV=selected.mass_MeV-threshold*HC,
        light_cut_MeV=2qp.m.u*HC, strange_cut_MeV=2qp.m.s*HC,
        gap_light_MeV=selected.mass_MeV-2qp.m.u*HC,
        gap_strange_MeV=selected.mass_MeV-2qp.m.s*HC))
    next_seed = selected.qualified ? (selected.mass_MeV/HC, selected.gamma_MeV/HC) : previous
    return candidates, summary, next_seed
end

function rate_contributions(rows, rates)
    out = NamedTuple[]
    for row in rows
        rr = filter(r -> String(r.point_id) == String(row.point_id), rates)
        length(rr) == length(RTM.REQUIRED_RATE_KEYS_FOR_TAU) || error("incomplete existing channel rates")
        values = Dict(Symbol(r.channel) => Float64(r.rate) for r in rr)
        densities = (; (q => Float64(getproperty(row, Symbol("n_", q))) for q in (:u,:d,:s,:ubar,:dbar,:sbar))...)
        complete = RTM.relaxation_rates(densities, values)
        for species in (:u,:s,:ubar,:sbar)
            expected = Float64(getproperty(row, Symbol("tauinv_",species)))
            abs(getproperty(complete,species)/expected-1) < 1e-10 || error("rate reconstruction mismatch")
        end
        for r in rr
            one = Dict(k => (k == Symbol(r.channel) ? Float64(r.rate) : 0.0) for k in keys(values))
            partial = RTM.relaxation_rates(densities, one)
            for species in (:u,:s,:ubar,:sbar)
                push!(out, merge(point_fields(row), (species=String(species), channel=String(r.channel),
                    rate_fm2=Float64(r.rate), inverse_tau_contribution_inv_fm=getproperty(partial,species),
                    inverse_tau_inv_fm=getproperty(complete,species),
                    fraction=getproperty(partial,species)/getproperty(complete,species))))
            end
        end
    end
    out
end

function main()
    profile = Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    ap = length(ARGS) == 6 ? parse(Int, ARGS[5]) : profile.a_integral.p_nodes
    ac = length(ARGS) == 6 ? parse(Int, ARGS[6]) : profile.a_integral.cos_nodes
    rows = collect(CSV.File(joinpath(INPUT_DIR,"source_rows.csv")))
    rates = collect(CSV.File(joinpath(INPUT_DIR,"source_rates.csv")))
    wanted = length(ARGS) >= 4 && ARGS[4] != "all" ? Set(split(ARGS[4],',')) : nothing
    selected = filter(r -> wanted === nothing || String(r.point_id) in wanted, rows)
    isempty(selected) && error("no selected points")
    states = Dict(String(r.point_id) => frozen_state(r,profile,ap,ac) for r in rows)
    thresholds = reduce(vcat,[threshold_rows(r,states[String(r.point_id)]) for r in rows])
    CSV.write(joinpath(OUTPUT_DIR,"threshold_inverses.csv"),thresholds)
    CSV.write(joinpath(OUTPUT_DIR,"channel_contributions.csv"),rate_contributions(rows,rates))
    CSV.write(joinpath(OUTPUT_DIR,"background_A.csv"),[merge(point_fields(r),
        (A_u=states[String(r.point_id)].qp.A.u,A_s=states[String(r.point_id)].qp.A.s,
        m_u_MeV=Float64(r.m_u)*HC,m_s_MeV=Float64(r.m_s)*HC)) for r in rows])
    candidates, summaries = NamedTuple[], NamedTuple[]
    for direction in ("forward","reverse"), window in unique(String(r.window) for r in selected)
        scan = sort(filter(r -> String(r.window)==window,selected); by=r -> Float64(r.xi), rev=direction=="reverse")
        previous = Dict{Symbol,Tuple{Float64,Float64}}()
        for row in scan
            started = time()
            for meson in MESONS
                cs, sm, next_seed = probe_roots(row, states[String(row.point_id)], meson,
                    direction, get(previous, meson, nothing))
                append!(candidates, cs)
                push!(summaries, sm)
                next_seed === nothing || (previous[meson]=next_seed)
            end
            CSV.write(joinpath(OUTPUT_DIR,"pole_candidates.csv"),candidates)
            CSV.write(joinpath(OUTPUT_DIR,"poles.csv"),summaries)
            println("completed ",direction," ",row.point_id," xi=",row.xi," seconds=",round(time()-started;digits=2))
            flush(stdout)
        end
    end
    files = Dict(name => digest(joinpath(OUTPUT_DIR,name)) for name in readdir(OUTPUT_DIR))
    write_json(joinpath(OUTPUT_DIR,"summary.json"), (
        status="completed_diagnostic", timestamp_utc=string(now(UTC)),
        julia_version=string(VERSION), julia_threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(), hbarc_MeV_fm=HC,
        source_root=SOURCE_ROOT, input_rows_sha256=digest(joinpath(INPUT_DIR,"source_rows.csv")),
        input_rates_sha256=digest(joinpath(INPUT_DIR,"source_rates.csv")),
        analysis_script_sha256=digest(@__FILE__),
        effective_profile=Models.TransportSamplingProfiles.effective_config(profile),
        effective_A=(p_nodes=ap,cos_nodes=ac,scale_inv_fm=profile.a_integral.scale_inv_fm),
        background_points=length(rows), root_points=length(selected),
        pole_rows=length(summaries), qualified_pole_rows=count(r->r.qualified,summaries),
        ambiguous_pole_rows=count(r->r.unique_qualified_roots>1,summaries),
        root_residual_limit=ROOT_TOL, root_method="existing real B0 plus width-coupled pole equation",
        mixed_identity="algebraic eigenbranches; multistart and bidirectional diagnostic only",
        rate_reconstruction_relative_tolerance=1e-10, output_sha256=files,
    ))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
