#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const CONSTANTS_SCRIPT = joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl")
const MODELS_SCRIPT = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")

isfile(CONSTANTS_SCRIPT) || error("missing constants script: $CONSTANTS_SCRIPT")
isfile(MODELS_SCRIPT) || error("missing models script: $MODELS_SCRIPT")

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, CONSTANTS_SCRIPT)
end
if !isdefined(Main, :Models)
    Base.include(Main, MODELS_SCRIPT)
end

using Main.MesonMass: solve_meson_mass

const WF = Main.Models.meson_workflow_module()
const HBARC = Main.Constants_PNJL.ħc_MeV_fm

function _read_bad_Ts(path::String)
    Ts = Float64[]
    i = 0
    for ln in eachline(path)
        i += 1
        i == 1 && continue
        s = strip(ln)
        isempty(s) && continue
        cols = split(s, ',')
        push!(Ts, parse(Float64, cols[1]))
    end
    return sort(unique(Ts))
end

function _as_namedtuple_params(quark, thermo)
    qp = (
        m=(u=Float64(quark.m.u), d=Float64(quark.m.d), s=Float64(quark.m.s)),
        μ=(u=Float64(quark.μ.u), d=Float64(quark.μ.d), s=Float64(quark.μ.s)),
    )
    tp = (
        T=Float64(thermo.T),
        Φ=Float64(thermo.Φ),
        Φbar=Float64(thermo.Φbar),
        ξ=Float64(thermo.ξ),
    )
    return qp, tp
end

function _run_eta_prime(qp, tp, m0::Float64, g0::Float64; method::Union{Nothing,Symbol}=nothing)
    kwargs = method === nothing ? (;) : (method=method,)
    return try
        solve_meson_mass(
            :eta_prime,
            qp,
            tp;
            k_norm=0.0,
            initial_mass=m0,
            initial_gamma=max(g0, 0.0),
            iterations=50,
            kwargs...,
        )
    catch
        nothing
    end
end

function _candidates(base_mass::Float64, base_gamma::Float64, prev_good::Union{Nothing,Vector{Float64}})
    cand = Tuple{Float64,Float64,String,Union{Nothing,Symbol}}[]

    push!(cand, (base_mass, base_gamma, "baseline", nothing))
    push!(cand, (base_mass, max(base_gamma, 0.0), "baseline_tr", :trust_region))
    push!(cand, (base_mass, 1.5 * max(base_gamma, 0.0), "baseline_gamma1p5", nothing))
    push!(cand, (base_mass, 2.0 * max(base_gamma, 0.0), "baseline_gamma2p0", :trust_region))

    if prev_good !== nothing
        m = prev_good[1]
        g = prev_good[2]
        push!(cand, (m, g, "prev_good", nothing))
        push!(cand, (m, g, "prev_good_tr", :trust_region))
        push!(cand, (m * 0.98, g * 1.2, "prev_good_scale_a", nothing))
        push!(cand, (m * 1.02, g * 0.8, "prev_good_scale_b", nothing))
        push!(cand, (m, g * 1.8, "prev_good_gamma1p8", :trust_region))
        push!(cand, (m + 10.0 / HBARC, g, "prev_good_m_plus10", :trust_region))
        push!(cand, (m - 10.0 / HBARC, g, "prev_good_m_minus10", :trust_region))
    end

    return cand
end

function main()
    bad_csv = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics", "fortran_conv_julia_nonconv_eta_prime_muB600.csv")
    isfile(bad_csv) || error("missing mismatch csv: $bad_csv")
    bad_Ts = Set(_read_bad_Ts(bad_csv))

    all_Ts = collect(100.0:1.0:300.0)
    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "experimental_eta_prime_strategy_boost_muB600.csv")
    out_txt = joinpath(out_dir, "experimental_eta_prime_strategy_boost_muB600_summary.txt")

    seed_state = nothing
    tracking_state = nothing
    prev_good = nothing

    rows = Vector{NamedTuple}()

    for T_MeV in all_Ts
        T_fm = T_MeV / HBARC
        muq_fm = (600.0 / 3.0) / HBARC

        base = WF.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=40,),
            meson_seed_state=seed_state,
            mixed_seed_tracking_state=tracking_state,
            mixed_branch_align=:identity_track_label_output,
        )

        seed_state = base.meson_seed_state
        tracking_state = base.mixed_seed_tracking

        base_res = base.meson_results[:eta_prime]
        qp, tp = _as_namedtuple_params(base.quark_params, base.thermo_params)

        best = nothing
        best_name = "none"
        best_resid = Inf
        best_conv = false

        for (m0, g0, name, method) in _candidates(Float64(base_res.mass), Float64(base_res.gamma), prev_good)
            rr = _run_eta_prime(qp, tp, m0, g0; method=method)
            rr === nothing && continue
            rr_resid = Float64(rr.residual_norm)
            if isfinite(rr_resid) && rr_resid < best_resid
                best = rr
                best_name = name
                best_resid = rr_resid
                best_conv = Bool(rr.converged)
            end
        end

        if best !== nothing && best_conv && best_resid <= 1e-6
            prev_good = Float64[Float64(best.mass), Float64(best.gamma)]
        end

        push!(rows, (
            T_MeV=T_MeV,
            in_bad_window=T_MeV in bad_Ts,
            baseline_residual=Float64(base_res.residual),
            baseline_converged=Bool(base_res.converged),
            boosted_best_seed=best_name,
            boosted_residual=best_resid,
            boosted_converged=best_conv,
            boosted_recovered_strict=best !== nothing && best_conv && best_resid <= 1e-6,
            boosted_recovered_loose=best !== nothing && best_conv && best_resid <= 1e-4,
            boosted_mass_MeV=best === nothing ? NaN : Float64(best.mass) * HBARC,
            boosted_gamma_MeV=best === nothing ? NaN : Float64(best.gamma) * HBARC,
        ))
    end

    open(out_csv, "w") do io
        write(io, "T_MeV,in_bad_window,baseline_residual,baseline_converged,boosted_best_seed,boosted_residual,boosted_converged,boosted_recovered_strict,boosted_recovered_loose,boosted_mass_MeV,boosted_gamma_MeV\n")
        for r in rows
            write(io, string(
                r.T_MeV, ',', r.in_bad_window, ',', r.baseline_residual, ',', r.baseline_converged, ',',
                r.boosted_best_seed, ',', r.boosted_residual, ',', r.boosted_converged, ',',
                r.boosted_recovered_strict, ',', r.boosted_recovered_loose, ',',
                r.boosted_mass_MeV, ',', r.boosted_gamma_MeV,
                '\n',
            ))
        end
    end

    bad_rows = filter(r -> r.in_bad_window, rows)
    nbad = length(bad_rows)
    rec_strict = count(r -> r.boosted_recovered_strict, bad_rows)
    rec_loose = count(r -> r.boosted_recovered_loose, bad_rows)
    base_good = count(r -> r.baseline_converged, bad_rows)

    open(out_txt, "w") do io
        write(io, "n_bad=$nbad\n")
        write(io, "baseline_converged_in_bad=$base_good\n")
        write(io, "boosted_recovered_strict=$rec_strict\n")
        write(io, "boosted_recovered_loose=$rec_loose\n")
        if nbad > 0
            write(io, "boosted_strict_rate=$(rec_strict / nbad)\n")
            write(io, "boosted_loose_rate=$(rec_loose / nbad)\n")
            bmean = sum(r.baseline_residual for r in bad_rows) / nbad
            nmean = sum(r.boosted_residual for r in bad_rows) / nbad
            write(io, "baseline_residual_mean_bad=$bmean\n")
            write(io, "boosted_residual_mean_bad=$nmean\n")
        end
    end

    println("Wrote experimental boost table: $out_csv")
    println("Wrote experimental boost summary: $out_txt")
end

main()
