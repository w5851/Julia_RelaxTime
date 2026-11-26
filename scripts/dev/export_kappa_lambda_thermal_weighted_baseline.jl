"""Generate a candidate for the kappa/lambda thermal-weighted regression baseline.

The recorded baseline `tests/baselines/relaxtime/baseline_kappa_lambda_thermal_weighted_v2.csv`
was produced by a one-off export and never had a repository generator, so a
physics correction (`#365`, the `w.usbar -> w.subar` strange-quark partner fix in
`_species_omega_s`) left it stale.  This exporter reproduces both recorded
columns from the current source:

* `expected` — the fixture configuration of `literature_validation_helpers.jl`
  (48/36/36 nodes, 24 sigma points, sampling rtol 1e-3);
* `refined`  — the tighter configuration recorded in the baseline metadata
  (72/40/48 nodes, 48 sigma points, sampling rtol 2.5e-4).

It never writes into `tests/baselines/`: the target must not exist yet, so a
candidate is always a new file to be reviewed before any consumer moves.  The
`relative_velocity=:cm` convention of the recorded fixture is preserved.
"""

using CSV
using JSON3
using SHA

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Base.include(Main, joinpath(_ROOT, "tests", "validation", "relaxtime", "literature_validation_helpers.jl"))

const POINTS = ((210.0, 0.0), (250.0, 0.0), (210.0, 100.0), (250.0, 100.0))
# `lambda` is a net-charge response and is only defined in the recorded fixture
# at finite muB; at muB = 0 it is degenerate, so it must not be emitted there.
const KAPPA_FIELDS = ((:kappa_BB, :kappa_BB), (:kappa_QQ, :kappa_QQ), (:kappa_SS, :kappa_SS))
const LAMBDA_FIELDS = ((:lambda, :lambda),)
fields_for(muB::Float64) = muB == 0.0 ? KAPPA_FIELDS : (KAPPA_FIELDS..., LAMBDA_FIELDS...)

function _tau_kwargs(; p_nodes, angle_nodes, phi_nodes, n_sigma_points, sampling_rtol)
    return (
        p_nodes=p_nodes, angle_nodes=angle_nodes, phi_nodes=phi_nodes, n_sigma_points=n_sigma_points,
        sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
        sampling_config=Main.AverageScatteringRate.ThermalSamplingConfig(
            rtol=sampling_rtol, atol=1e-12, weighted_window=0.6, seed_order=48,
            max_evaluations=80_000, max_nodes=16_000, max_rounds=50,
            max_events=60_000_000, relative_velocity=:cm,
        ),
    )
end

const CONFIGS = (
    expected=_tau_kwargs(p_nodes=48, angle_nodes=36, phi_nodes=36, n_sigma_points=24, sampling_rtol=1e-3),
    refined=_tau_kwargs(p_nodes=72, angle_nodes=40, phi_nodes=48, n_sigma_points=48, sampling_rtol=2.5e-4),
)

const P_NUM = 12
const T_NUM = 6
const MAX_ITER = 40

const HBARC = Main.Constants_PNJL.ħc_MeV_fm

function transport_point(T_MeV::Float64, muB_MeV::Float64, tau_kwargs)
    T_fm = T_MeV / HBARC
    muq_fm = (muB_MeV / 3.0) / HBARC
    equilibrium = Main._solve_relaxtime_literature_validation_equilibrium(T_fm, muq_fm, 0.0)
    Bool(equilibrium.converged) || error("equilibrium solve failed at T=$(T_MeV) muB=$(muB_MeV)")
    masses = (u=Float64(equilibrium.masses[1]), d=Float64(equilibrium.masses[2]), s=Float64(equilibrium.masses[3]))
    ktmp = Main._build_relaxtime_literature_validation_K_coeffs(
        T_fm, muq_fm, masses, Float64(equilibrium.x_state[4]), Float64(equilibrium.x_state[5]))
    result = Main.RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW.solve_gap_and_transport(
        T_fm, muq_fm;
        xi=0.0, equilibrium=equilibrium, compute_tau=true, K_coeffs=ktmp.K_coeffs,
        compute_bulk=false, tau_kwargs=tau_kwargs, p_num=P_NUM, t_num=T_NUM,
        seed_state=Vector(equilibrium.x_state), solver_kwargs=(iterations=MAX_ITER,),
    )
    return result.transport
end

function main(args)
    length(args) == 1 || error("usage: export_kappa_lambda_thermal_weighted_baseline.jl <new-candidate.csv>")
    target = abspath(args[1])
    ispath(target) && error("refusing to overwrite existing target: $target")
    isdir(dirname(target)) || error("target directory does not exist: $(dirname(target))")

    rows = NamedTuple[]
    for (T, muB) in POINTS
        record = Dict{String, Float64}()
        for (column, mode) in ((:expected, :expected), (:refined, :refined))
            println("  computing T=$(T) muB=$(muB) config=$(mode)")
            transport = transport_point(T, muB, CONFIGS[mode])
            for (transport_field, csv_field) in fields_for(muB)
                record["$(csv_field).$(column)"] = Float64(getproperty(transport, transport_field))
            end
        end
        for (_, csv_field) in fields_for(muB)
            push!(rows, (record_id="$(csv_field)_mu$(Int(round(muB)))_t$(Int(T))",
                         T_MeV=T, muB_MeV=muB, field=String(csv_field),
                         expected=record["$(csv_field).expected"], refined=record["$(csv_field).refined"]))
        end
    end
    sort!(rows, by=row -> (row.muB_MeV, row.T_MeV, row.record_id))
    for row in rows
        isfinite(row.expected) && row.expected > 0 || error("non-finite expected for $(row.record_id)")
        isfinite(row.refined) && row.refined > 0 || error("non-finite refined for $(row.record_id)")
        isfinite(row.T_MeV) && isfinite(row.muB_MeV) || error("bad coordinates for $(row.record_id)")
    end
    CSV.write(target, rows; header=[:record_id, :T_MeV, :muB_MeV, :field, :expected, :refined])
    println("wrote $(length(rows)) candidate rows to $target")
    println("candidate sha256 = ", bytes2hex(sha256(read(target))))
end

main(ARGS)
