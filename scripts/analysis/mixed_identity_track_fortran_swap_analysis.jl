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

@inline _mass_rel(a::Float64, b::Float64) = abs(a - b) / max(abs(b), 1e-12)
@inline _gamma_rel_floor(a::Float64, b::Float64) = abs(a - b) / max(abs(b), 1.0)
@inline _gamma_rel_raw(a::Float64, b::Float64) = abs(a - b) / max(abs(b), 1e-12)

function main()
    root_diag_path = joinpath(FORTRAN_ROOT, "quark_phase", "root_diag.dat")
    isfile(root_diag_path) || error("missing fortran root diagnostics: $root_diag_path")

    fortran_rows = Dict{Float64,Dict{Symbol,NamedTuple}}()
    for r in _read_numeric_rows(root_diag_path)
        length(r) >= 14 || continue
        idx = Int(round(r[2]))
        idx in (3, 4) || continue
        isapprox(r[14], 600.0; atol=1e-9) || continue
        T = r[1]
        meson = idx == 3 ? :eta : :eta_prime
        d = get!(fortran_rows, T, Dict{Symbol,NamedTuple}())
        d[meson] = (mass_MeV=r[3], gamma_MeV=r[4])
    end

    Ts = sort(collect(keys(fortran_rows)))
    isempty(Ts) && error("no muB=600 eta/eta_prime rows found in $root_diag_path")

    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "fortran_julia_identity_track_swap_compare_muB600.csv")

    seed_state = nothing
    tracking_state = nothing

    open(out_csv, "w") do io
        write(io, "T_MeV,slot,fortran_mapping,j_mass_MeV,j_gamma_MeV,f_mass_MeV,f_gamma_MeV,mass_rel,gamma_rel_floor1,gamma_rel_raw\n")
        for T_MeV in Ts
            T_fm = T_MeV / HBARC
            muq_fm = (600.0 / 3.0) / HBARC

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

            j_eta = (mass_MeV=res.meson_results[:eta].mass * HBARC, gamma_MeV=res.meson_results[:eta].gamma * HBARC)
            j_etap = (mass_MeV=res.meson_results[:eta_prime].mass * HBARC, gamma_MeV=res.meson_results[:eta_prime].gamma * HBARC)

            f_eta = fortran_rows[T_MeV][:eta]
            f_etap = fortran_rows[T_MeV][:eta_prime]

            a_cost = _mass_rel(j_eta.mass_MeV, f_eta.mass_MeV) + _mass_rel(j_etap.mass_MeV, f_etap.mass_MeV)
            b_cost = _mass_rel(j_eta.mass_MeV, f_etap.mass_MeV) + _mass_rel(j_etap.mass_MeV, f_eta.mass_MeV)

            if b_cost < a_cost
                mapping = "swap"
                f1 = f_etap
                f2 = f_eta
            else
                mapping = "label"
                f1 = f_eta
                f2 = f_etap
            end

            mrel1 = _mass_rel(j_eta.mass_MeV, f1.mass_MeV)
            grel1 = _gamma_rel_floor(j_eta.gamma_MeV, f1.gamma_MeV)
            grel1_raw = _gamma_rel_raw(j_eta.gamma_MeV, f1.gamma_MeV)
            write(io, string(T_MeV, ',', :eta, ',', mapping, ',', j_eta.mass_MeV, ',', j_eta.gamma_MeV, ',', f1.mass_MeV, ',', f1.gamma_MeV, ',', mrel1, ',', grel1, ',', grel1_raw, '\n'))

            mrel2 = _mass_rel(j_etap.mass_MeV, f2.mass_MeV)
            grel2 = _gamma_rel_floor(j_etap.gamma_MeV, f2.gamma_MeV)
            grel2_raw = _gamma_rel_raw(j_etap.gamma_MeV, f2.gamma_MeV)
            write(io, string(T_MeV, ',', :eta_prime, ',', mapping, ',', j_etap.mass_MeV, ',', j_etap.gamma_MeV, ',', f2.mass_MeV, ',', f2.gamma_MeV, ',', mrel2, ',', grel2, ',', grel2_raw, '\n'))
        end
    end

    println("Wrote identity-track swap-compare CSV: $out_csv")
end

main()
