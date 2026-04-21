if !isdefined(Main, :Constants_PNJL)
    include(joinpath(@__DIR__, "..", "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
end

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "..", "src", "models", "Models.jl"))
end

const WF = Models.meson_workflow_module()
const HBARC = Main.Constants_PNJL.ħc_MeV_fm

const _MESONS = (:pi, :K, :eta, :eta_prime, :sigma_pi, :sigma_K, :sigma, :sigma_prime)
const _RTOL_BY_MESON = Dict(
    :pi => 1e-2,
    :K => 1e-2,
    :eta => 4e-2,
    :eta_prime => 3e-1,
    :sigma_pi => 1e-2,
    :sigma_K => 1e-2,
    :sigma => 1e-2,
    :sigma_prime => 1e-2,
)

@inline _julia_symbol_for_legacy(meson::Symbol) =
    meson == :eta ? :eta_prime :
    meson == :eta_prime ? :eta :
    meson == :sigma ? :sigma_prime :
    meson == :sigma_prime ? :sigma : meson

function _read_csv_dict_rows(path::String)
    lines = readlines(path)
    isempty(lines) && error("empty csv: $path")
    header = strip.(split(strip(lines[1]), ','))
    rows = Vector{Dict{String,String}}()
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = strip.(split(s, ','))
        length(cols) == length(header) || continue
        row = Dict{String,String}()
        for i in eachindex(header)
            row[header[i]] = cols[i]
        end
        push!(rows, row)
    end
    return rows
end

function _kcoeffs_from_qp(qp)
    Gu = Main.EffectiveCouplings.calculate_G_from_A(qp.A.u, qp.m.u)
    Gs = Main.EffectiveCouplings.calculate_G_from_A(qp.A.s, qp.m.s)
    return Main.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,
        Main.Constants_PNJL.K_fm5,
        Gu,
        Gs,
    )
end

function _min_residual_over_gamma(meson::Symbol, m_fm::Float64, qp, tp, Kc)
    best_gamma = 0.0
    best_resid = Inf
    g = 0.0
    while g <= 2.0 + 1e-12
        f = Main.MesonMass.meson_mass_equation(meson, m_fm, g, 0.0, qp, tp, Kc)
        r = hypot(real(f), imag(f))
        if isfinite(r) && r < best_resid
            best_resid = r
            best_gamma = g
        end
        g += 0.01
    end
    return (gamma=best_gamma, resid=best_resid)
end

function _build_qp_tp_from_workflow_result(res)
    qp_in = (
        m=(
            u=Float64(res.quark_params.m.u),
            d=Float64(res.quark_params.m.d),
            s=Float64(res.quark_params.m.s),
        ),
        μ=(
            u=Float64(res.quark_params.μ.u),
            d=Float64(res.quark_params.μ.d),
            s=Float64(res.quark_params.μ.s),
        ),
    )
    tp = (
        T=Float64(res.thermo_params.T),
        Φ=Float64(res.thermo_params.Φ),
        Φbar=Float64(res.thermo_params.Φbar),
        ξ=Float64(res.thermo_params.ξ),
    )
    qp = Main.MesonMass.ensure_quark_params_has_A(qp_in, tp)
    return qp, tp
end

function run_literature_targets()
    path = joinpath(@__DIR__, "..", "..", "..", "..", "tests", "validation", "data", "targets", "relaxtime", "literature", "meson", "relaxtime_meson_mass_literature_targets_v1.csv")
    rows = _read_csv_dict_rows(path)

    failed = 0
    near_root = 0
    println("[literature]")
    println("id,meson,T_MeV,muB_MeV,expected_MeV,actual_MeV,solver_resid,expected_gamma0_resid,min_resid_fixed_mass,min_gamma")

    for row in rows
        meson = Symbol(row["meson"])
        T_MeV = parse(Float64, row["T_MeV"])
        muB_MeV = parse(Float64, row["muB_MeV"])
        xi = parse(Float64, row["xi"])
        expected_MeV = parse(Float64, row["expected_mass_MeV"])
        rtol = parse(Float64, row["rtol"])

        T_fm = T_MeV / HBARC
        muq_fm = (muB_MeV / 3.0) / HBARC
        res = WF.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=xi,
            mesons=(meson,),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
        )
        mr = res.meson_results[meson]
        actual_MeV = Float64(mr.mass) * HBARC

        if !(isfinite(actual_MeV) && isapprox(actual_MeV, expected_MeV; rtol=rtol, atol=0.0))
            failed += 1
            expected_fm = expected_MeV / HBARC
            qp, tp = _build_qp_tp_from_workflow_result(res)
            Kc = _kcoeffs_from_qp(qp)

            f0 = Main.MesonMass.meson_mass_equation(meson, expected_fm, 0.0, 0.0, qp, tp, Kc)
            r0 = hypot(real(f0), imag(f0))
            best = _min_residual_over_gamma(meson, expected_fm, qp, tp, Kc)
            if best.resid < 1e-2
                near_root += 1
            end

            println(string(
                row["target_id"], ",", meson, ",", T_MeV, ",", muB_MeV, ",", expected_MeV, ",",
                actual_MeV, ",", Float64(mr.residual), ",", r0, ",", best.resid, ",", best.gamma,
            ))
        end
    end

    println("failed_rows=", failed)
    println("failed_rows_with_near_root_fixed_mass(min_resid<1e-2)=", near_root)
end

function run_legacy_all8_targets()
    paths = (
        joinpath(@__DIR__, "..", "..", "..", "..", "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB0_v1.csv"),
        joinpath(@__DIR__, "..", "..", "..", "..", "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB600_v1.csv"),
    )

    targets = Dict{Tuple{Float64,Float64,Symbol},Float64}()
    for path in paths
        for row in _read_csv_dict_rows(path)
            get(row, "source_impl", "") == "fortran" || continue
            xi = parse(Float64, row["xi"])
            isapprox(xi, 0.0; atol=1e-12) || continue
            if get(row, "solver_status", "") == "excluded_low_quality_nonconverged"
                continue
            end
            meson = Symbol(row["meson"])
            meson in _MESONS || continue

            muB = parse(Float64, row["muB_MeV"])
            T = parse(Float64, row["T_MeV"])
            mass = parse(Float64, row["mass_MeV"])
            targets[(muB, T, meson)] = mass
        end
    end

    expected_mus = (0.0, 600.0)
    expected_Ts = (120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0)

    failed = 0
    near_root = 0
    println("[legacy-all8]")
    println("muB,T,legacy_meson,julia_meson,legacy_MeV,actual_MeV,solver_resid,expected_gamma0_resid,min_resid_fixed_mass,min_gamma")

    for muB in expected_mus
        for T in expected_Ts
            T_fm = T / HBARC
            muq_fm = (muB / 3.0) / HBARC
            res = WF.solve_gap_and_meson_point(
                T_fm,
                muq_fm;
                xi=0.0,
                mesons=_MESONS,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=25,),
                mass_kwargs=(iterations=25,),
            )

            qp, tp = _build_qp_tp_from_workflow_result(res)
            Kc = _kcoeffs_from_qp(qp)

            for meson in _MESONS
                key = (muB, T, meson)
                haskey(targets, key) || continue
                legacy_mass = targets[key]
                julia_meson = _julia_symbol_for_legacy(meson)
                mr = res.meson_results[julia_meson]
                actual_MeV = Float64(mr.mass) * HBARC

                if !(isfinite(actual_MeV) && isapprox(actual_MeV, legacy_mass; rtol=_RTOL_BY_MESON[meson], atol=0.0))
                    failed += 1
                    expected_fm = legacy_mass / HBARC

                    f0 = Main.MesonMass.meson_mass_equation(julia_meson, expected_fm, 0.0, 0.0, qp, tp, Kc)
                    r0 = hypot(real(f0), imag(f0))
                    best = _min_residual_over_gamma(julia_meson, expected_fm, qp, tp, Kc)
                    if best.resid < 1e-2
                        near_root += 1
                    end

                    println(string(
                        muB, ",", T, ",", meson, ",", julia_meson, ",", legacy_mass, ",",
                        actual_MeV, ",", Float64(mr.residual), ",", r0, ",", best.resid, ",", best.gamma,
                    ))
                end
            end
        end
    end

    println("failed_rows=", failed)
    println("failed_rows_with_near_root_fixed_mass(min_resid<1e-2)=", near_root)
end

run_literature_targets()
run_legacy_all8_targets()
