#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FORTRAN_ROOT = raw"D:\Desktop\fortran_mott\PNJL-mott-mu_T\PNJL-mu-T"

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

function _read_numeric_rows(path::String)
    rows = Vector{Vector{Float64}}()
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        vals = Float64[]
        ok = true
        for c in split(s)
            v = tryparse(Float64, c)
            if v === nothing
                ok = false
                break
            end
            push!(vals, v)
        end
        ok || continue
        push!(rows, vals)
    end
    return rows
end

function _load_fortran_eta_prime_by_T()
    root_diag_path = joinpath(FORTRAN_ROOT, "quark_phase", "root_diag.dat")
    out = Dict{Float64,NamedTuple}()
    for r in _read_numeric_rows(root_diag_path)
        length(r) >= 14 || continue
        idx = Int(round(r[2]))
        idx == 4 || continue
        isapprox(r[14], 600.0; atol=1e-9) || continue
        out[r[1]] = (mass_MeV=r[3], gamma_MeV=r[4])
    end
    return out
end

function _load_bad_Ts(path::String)
    Ts = Float64[]
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        s = strip(ln)
        isempty(s) && continue
        cols = split(s, ',')
        push!(Ts, parse(Float64, cols[1]))
    end
    return sort(unique(Ts))
end

function _try_one_seed(quark_params, thermo_params, m0_fm::Float64, g0_fm::Float64; method::Union{Nothing,Symbol}=nothing)
    kwargs = method === nothing ? (;) : (method=method,)
    return try
        solve_meson_mass(
            :eta_prime,
            quark_params,
            thermo_params;
            k_norm=0.0,
            initial_mass=m0_fm,
            initial_gamma=max(g0_fm, 0.0),
            iterations=40,
            kwargs...,
        )
    catch
        nothing
    end
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

function main()
    bad_csv = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics", "fortran_conv_julia_nonconv_eta_prime_muB600.csv")
    isfile(bad_csv) || error("missing bad-case csv: $bad_csv")

    bad_Ts = Set(_load_bad_Ts(bad_csv))
    ft_etap = _load_fortran_eta_prime_by_T()

    Ts = sort(collect(keys(ft_etap)))
    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "eta_prime_anchor_recovery_experiment_muB600.csv")
    out_summary = joinpath(out_dir, "eta_prime_anchor_recovery_experiment_muB600_summary.txt")

    seed_state = nothing
    tracking_state = nothing

    rows = Vector{NamedTuple}()
    recovered_strict = 0
    recovered_loose = 0

    for T_MeV in Ts
        T_fm = T_MeV / HBARC
        muq_fm = (600.0 / 3.0) / HBARC

        prev_eta_plus = tracking_state === nothing ? nothing : get(tracking_state, :eta_plus, nothing)

        res = WF.solve_gap_and_meson_point(
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

        seed_state = res.meson_seed_state
        tracking_state = res.mixed_seed_tracking

        in_bad = T_MeV in bad_Ts
        in_bad || continue

        jl = res.meson_results[:eta_prime]
        qp, tp = _as_namedtuple_params(res.quark_params, res.thermo_params)
        ft = ft_etap[T_MeV]

        candidates = Tuple{Float64,Float64,String,Union{Nothing,Symbol}}[]

        if prev_eta_plus !== nothing
            m = prev_eta_plus[1]
            g = prev_eta_plus[2]
            push!(candidates, (m, g, "prev_plus", nothing))
            push!(candidates, (m * 0.98, g * 1.02, "prev_plus_scale_a", nothing))
            push!(candidates, (m * 1.02, g * 0.98, "prev_plus_scale_b", nothing))
            push!(candidates, (m, g * 0.5, "prev_plus_half_gamma", nothing))
            push!(candidates, (m, g * 1.5, "prev_plus_1p5_gamma", nothing))
            push!(candidates, ((m * HBARC + 10.0) / HBARC, g, "prev_plus_m_plus10", nothing))
            push!(candidates, ((m * HBARC - 10.0) / HBARC, g, "prev_plus_m_minus10", nothing))
        end

        push!(candidates, (ft.mass_MeV / HBARC, max(ft.gamma_MeV, 0.0) / HBARC, "fortran_eta_prime", nothing))
        push!(candidates, (ft.mass_MeV / HBARC, max(ft.gamma_MeV, 0.0) / HBARC, "fortran_eta_prime_tr", :trust_region))

        best = nothing
        best_resid = Inf
        best_name = "none"
        best_conv = false
        best_mass_MeV = NaN
        best_gamma_MeV = NaN

        for (m0, g0, name, method) in candidates
            rr = _try_one_seed(qp, tp, m0, g0; method=method)
            rr === nothing && continue
            resid = Float64(rr.residual_norm)
            if isfinite(resid) && resid < best_resid
                best = rr
                best_resid = resid
                best_name = name
                best_conv = Bool(rr.converged)
                best_mass_MeV = Float64(rr.mass) * HBARC
                best_gamma_MeV = Float64(rr.gamma) * HBARC
            end
        end

        strict_ok = best !== nothing && best_conv && best_resid <= 1e-6
        loose_ok = best !== nothing && best_conv && best_resid <= 1e-4
        strict_ok && (recovered_strict += 1)
        loose_ok && (recovered_loose += 1)

        push!(rows, (
            T_MeV=T_MeV,
            baseline_residual=Float64(jl.residual),
            baseline_converged=Bool(jl.converged),
            baseline_mass_MeV=Float64(jl.mass) * HBARC,
            baseline_gamma_MeV=Float64(jl.gamma) * HBARC,
            best_seed=best_name,
            best_residual=best_resid,
            best_converged=best_conv,
            best_mass_MeV=best_mass_MeV,
            best_gamma_MeV=best_gamma_MeV,
            recovered_strict=strict_ok,
            recovered_loose=loose_ok,
        ))
    end

    open(out_csv, "w") do io
        write(io, "T_MeV,baseline_residual,baseline_converged,baseline_mass_MeV,baseline_gamma_MeV,best_seed,best_residual,best_converged,best_mass_MeV,best_gamma_MeV,recovered_strict,recovered_loose\n")
        for r in rows
            write(io, string(
                r.T_MeV, ',', r.baseline_residual, ',', r.baseline_converged, ',',
                r.baseline_mass_MeV, ',', r.baseline_gamma_MeV, ',',
                r.best_seed, ',', r.best_residual, ',', r.best_converged, ',',
                r.best_mass_MeV, ',', r.best_gamma_MeV, ',',
                r.recovered_strict, ',', r.recovered_loose,
                '\n',
            ))
        end
    end

    n = length(rows)
    open(out_summary, "w") do io
        write(io, "n_bad=$n\n")
        write(io, "recovered_strict=$recovered_strict\n")
        write(io, "recovered_loose=$recovered_loose\n")
        if n > 0
            write(io, "strict_rate=$(recovered_strict / n)\n")
            write(io, "loose_rate=$(recovered_loose / n)\n")
        end
    end

    println("Wrote anchor recovery table: $out_csv")
    println("Wrote anchor recovery summary: $out_summary")
end

main()
