#!/usr/bin/env julia
# Isolate fixed-background sensitivity from a change of the production solver.
using CSV, JSON, SHA
include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
include(joinpath(@__DIR__, "ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon: start_run, finish_run, write_rows, write_json

function diagnose_background_replay(base, replaydir, outdir, point, process)
    summary = JSON.parsefile(joinpath(base, "joint_complete_summary_v1", "summary.json"))
    matches = NamedTuple[]
    for source in summary["sources"]
        directory = source["directory"]
        backgrounds = filter(r -> r.point_id == point && r.configuration == "candidate",
            collect(CSV.File(joinpath(directory, "background.csv"))))
        isempty(backgrounds) && continue
        length(backgrounds) == 1 || error("duplicate qualified background")
        rates = filter(r -> r.point_id == point && r.configuration == "candidate" &&
            r.channel == String(process), collect(CSV.File(joinpath(directory, "rates.csv"))))
        push!(matches, (directory=directory, background=only(backgrounds), rate=only(rates)))
    end
    original = only(matches)
    folder = joinpath(replaydir, "points", point)
    production = only(CSV.File(joinpath(folder, "phase_guided_transport_scan.csv"); comment="#"))
    production_rate = only(filter(r -> r.channel == String(process),
        collect(CSV.File(joinpath(folder, "phase_guided_transport_scan_sampling.csv")))))
    profile_api = Models.TransportSamplingProfiles
    profile = profile_api.load_profile("thermal_weighted_v1")
    runtime = profile_api.build_runtime(profile)
    inputs = [joinpath(original.directory, "background.csv"),
        joinpath(original.directory, "rates.csv"),
        joinpath(folder, "phase_guided_transport_scan.csv"),
        joinpath(folder, "phase_guided_transport_scan_sampling.csv")]
    start_run(outdir, joinpath(base, "convergence_selection_v1"), @__FILE__,
        (purpose="same native channel on frozen qualification and production backgrounds",
         point_id=point, channel=String(process), profile=profile_api.effective_config(profile),
         sources=[(path=abspath(p), sha256=bytes2hex(sha256(read(p)))) for p in inputs]))
    for path in inputs
        cp(path, joinpath(outdir, "inputs", basename(path)))
    end

    function state(row; equal_light_masses=false)
        masses = Float64.([row.m_u, row.m_d, row.m_s])
        equal_light_masses && (masses[1] = masses[2] = (masses[1] + masses[2]) / 2)
        equilibrium = (x_state=[0.0, 0.0, 0.0, Float64(row.Phi), Float64(row.Phibar)],
            masses=masses, mu_vec=fill(Float64(production.muq_fm), 3))
        return profile_api.point_state(profile, equilibrium, production.T_fm,
            production.muq_fm; xi=production.xi)
    end
    settings = [
        (name="qualified_background", state=state(original.background), expected=Float64(original.rate.rate)),
        (name="production_background", state=state(production), expected=Float64(production_rate.rate)),
        (name="production_background_repeat", state=state(production), expected=Float64(production_rate.rate)),
        (name="qualified_background_equal_ud", state=state(original.background; equal_light_masses=true), expected=NaN),
    ]
    rows = NamedTuple[]
    for setting in settings
        println("START fixed background ", setting.name); flush(stdout)
        GC.gc()
        Main.RelaxTime.PolarizationCache.reset_cache!()
        s = setting.state
        measured = @timed profile_api.AS.build_thermal_weighted_cache(process,
            s.quark_params, s.thermo_params, s.K_coeffs;
            runtime.tau_kwargs..., propagator_quark_params=s.propagator_quark_params)
        cache = measured.value
        error = isfinite(setting.expected) ? cache.rate / setting.expected - 1 : NaN
        push!(rows, merge((variant=setting.name, rate=cache.rate, expected_rate=setting.expected,
            replay_relative_change=error, seconds=measured.time,
            m_u=s.quark_params.m.u, m_d=s.quark_params.m.d, m_s=s.quark_params.m.s,
            Phi=s.thermo_params.Φ, Phibar=s.thermo_params.Φbar,
            estimated_error=cache.estimated_error), cache.statistics))
        write_rows(joinpath(outdir, "rates.csv"), rows)
        write_rows(joinpath(outdir, setting.name * "_table.csv"),
            [(s=x, sigma=y) for (x, y) in zip(cache.table.nodes, cache.table.sigma)])
        write_json(joinpath(outdir, setting.name * "_state.json"), s)
        println("DONE fixed background ", setting.name, " relative=", error); flush(stdout)
    end
    original_replay = abs(rows[1].replay_relative_change) <= 1e-10
    production_replay = all(abs(r.replay_relative_change) <= 1e-10 for r in rows[2:3])
    repeat_exact = rows[2].rate == rows[3].rate
    finish_run(outdir; points=1, channel=String(process), variants=length(rows),
        qualified_background_reproduced=original_replay,
        production_background_reproduced=production_replay,
        repeated_production_rate_identical=repeat_exact,
        diagnosis_passed=original_replay && production_replay && repeat_exact)
    original_replay && production_replay && repeat_exact ||
        error("frozen-background replays differ; inspect rates and retained states")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    3 <= length(ARGS) <= 5 || error("Usage: run_transport_background_replay_diagnostic.jl BASE REPLAY_DIR NEW_OUTPUT_DIR [POINT_ID CHANNEL]")
    point = length(ARGS) >= 4 ? ARGS[4] : "muB0.0_alpha1.0_xi+0.3700"
    process = length(ARGS) >= 5 ? Symbol(ARGS[5]) : :uubar_to_ddbar
    diagnose_background_replay(abspath(ARGS[1]), abspath(ARGS[2]), abspath(ARGS[3]), point, process)
end
