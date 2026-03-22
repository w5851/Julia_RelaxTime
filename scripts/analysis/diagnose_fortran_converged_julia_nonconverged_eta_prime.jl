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

const WF = Main.Models.meson_workflow_module()
const HBARC = Main.Constants_PNJL.ħc_MeV_fm

function _read_numeric_rows(path::String)
    rows = Vector{Vector{Float64}}()
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        cols = split(s)
        vals = Float64[]
        ok = true
        for c in cols
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

@inline function _seed_dist(m::Float64, g::Float64, seed::Union{Nothing,Vector{Float64}})
    seed === nothing && return NaN
    return hypot(m - seed[1], g - seed[2])
end

function main()
    root_diag_path = joinpath(FORTRAN_ROOT, "quark_phase", "root_diag.dat")
    isfile(root_diag_path) || error("missing fortran root diagnostics: $root_diag_path")

    frows = Dict{Float64,Dict{Symbol,NamedTuple}}()
    for r in _read_numeric_rows(root_diag_path)
        length(r) >= 14 || continue
        idx = Int(round(r[2]))
        idx in (3, 4) || continue
        isapprox(r[14], 600.0; atol=1e-9) || continue
        meson = idx == 3 ? :eta : :eta_prime
        d = get!(frows, r[1], Dict{Symbol,NamedTuple}())
        d[meson] = (
            mass_MeV=r[3],
            gamma_MeV=r[4],
            resid=r[5],
            iters=Int(round(r[6])),
            converged=Int(round(r[7])) == 1,
        )
    end

    Ts = sort(collect(keys(frows)))
    isempty(Ts) && error("no muB=600 eta/eta_prime rows found")

    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "fortran_conv_julia_nonconv_eta_prime_muB600.csv")
    summary_txt = joinpath(out_dir, "fortran_conv_julia_nonconv_eta_prime_muB600_summary.txt")

    seed_state = nothing
    tracking_state = nothing
    bad_rows = Vector{NamedTuple}()

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

        jl = res.meson_results[:eta_prime]
        ft = frows[T_MeV][:eta_prime]

        if ft.converged && !Bool(jl.converged)
            jm = Float64(jl.mass) * HBARC
            jg = Float64(jl.gamma) * HBARC
            push!(bad_rows, (
                T_MeV=T_MeV,
                fortran_resid=Float64(ft.resid),
                fortran_iters=Int(ft.iters),
                fortran_mass_MeV=Float64(ft.mass_MeV),
                fortran_gamma_MeV=Float64(ft.gamma_MeV),
                julia_residual=Float64(jl.residual),
                julia_root_quality=String(jl.root_quality),
                julia_mass_MeV=jm,
                julia_gamma_MeV=jg,
                branch_score_selected=Float64(jl.branch_score_selected),
                branch_score_plus=Float64(jl.branch_score_plus),
                branch_score_minus=Float64(jl.branch_score_minus),
                dist_to_prev_eta_plus=_seed_dist(jm, jg, prev_eta_plus),
            ))
        end
    end

    open(out_csv, "w") do io
        write(io, "T_MeV,fortran_resid,fortran_iters,fortran_mass_MeV,fortran_gamma_MeV,julia_residual,julia_root_quality,julia_mass_MeV,julia_gamma_MeV,branch_score_selected,branch_score_plus,branch_score_minus,dist_to_prev_eta_plus\n")
        for r in bad_rows
            write(io, string(
                r.T_MeV, ',', r.fortran_resid, ',', r.fortran_iters, ',',
                r.fortran_mass_MeV, ',', r.fortran_gamma_MeV, ',',
                r.julia_residual, ',', r.julia_root_quality, ',',
                r.julia_mass_MeV, ',', r.julia_gamma_MeV, ',',
                r.branch_score_selected, ',', r.branch_score_plus, ',', r.branch_score_minus, ',',
                r.dist_to_prev_eta_plus,
                '\n',
            ))
        end
    end

    n = length(bad_rows)
    if n == 0
        write(summary_txt, "No Fortran-converged/Julia-nonconverged eta_prime rows found.\n")
        println("Wrote diagnostics: $out_csv")
        println("Wrote summary: $summary_txt")
        return
    end

    Ts_bad = [r.T_MeV for r in bad_rows]
    jl_res = [r.julia_residual for r in bad_rows]
    ft_res = [r.fortran_resid for r in bad_rows]
    dseed = [r.dist_to_prev_eta_plus for r in bad_rows if isfinite(r.dist_to_prev_eta_plus)]
    branch_sel = [r.branch_score_selected for r in bad_rows]

    open(summary_txt, "w") do io
        write(io, "count=$n\n")
        write(io, "T_range=$(minimum(Ts_bad))..$(maximum(Ts_bad))\n")
        write(io, "fortran_resid_mean=$(sum(ft_res)/length(ft_res))\n")
        write(io, "fortran_resid_max=$(maximum(ft_res))\n")
        write(io, "julia_residual_mean=$(sum(jl_res)/length(jl_res))\n")
        write(io, "julia_residual_max=$(maximum(jl_res))\n")
        write(io, "branch_selected_mean=$(sum(branch_sel)/length(branch_sel))\n")
        write(io, "branch_selected_max=$(maximum(branch_sel))\n")
        if !isempty(dseed)
            write(io, "dist_to_prev_eta_plus_mean=$(sum(dseed)/length(dseed))\n")
            write(io, "dist_to_prev_eta_plus_max=$(maximum(dseed))\n")
        else
            write(io, "dist_to_prev_eta_plus_mean=NaN\n")
            write(io, "dist_to_prev_eta_plus_max=NaN\n")
        end
    end

    println("Wrote diagnostics: $out_csv")
    println("Wrote summary: $summary_txt")
end

main()
