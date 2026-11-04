#!/usr/bin/env julia
# Reuse the qualified multistart pole diagnostic on the frozen transport grid.
# Arguments match run_transport_mott_probe.jl. Execute only on quote-remote.
include(joinpath(@__DIR__, "run_transport_mott_probe.jl"))

function eta_survey_thresholds(row, state)
    out = NamedTuple[]
    qp, tp, kc = state.qp, state.tp, state.kc
    for meson in (:eta, :eta_prime, :pi, :K)
        cuts = meson in (:eta, :eta_prime) ? ((:uu, 2qp.m.u), (:ss, 2qp.m.s)) :
            ((:pair, MT.mott_threshold_mass(meson, qp)),)
        for (cut, threshold) in cuts, offset in (-1e-4, -1e-6)
            f = MP.meson_mass_equation(meson, threshold + offset, 0.0, 0.0, qp, tp, kc)
            push!(out, merge(point_fields(row), (meson=String(meson), cut=String(cut),
                threshold_MeV=threshold*HC, offset_inv_fm=offset,
                inverse_re=real(f), inverse_im=imag(f), inverse_abs=abs(f),
                inverse_units=meson in (:eta, :eta_prime) ? "fm2" : "dimensionless")))
        end
    end
    return out
end

function run_eta_survey()
    profile = Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    ap = length(ARGS) == 6 ? parse(Int, ARGS[5]) : profile.a_integral.p_nodes
    ac = length(ARGS) == 6 ? parse(Int, ARGS[6]) : profile.a_integral.cos_nodes
    rows = collect(CSV.File(joinpath(INPUT_DIR, "source_rows.csv")))
    wanted = length(ARGS) >= 4 && ARGS[4] != "all" ? Set(split(ARGS[4], ',')) : nothing
    selected = filter(r -> wanted === nothing || String(r.point_id) in wanted, rows)
    isempty(selected) && error("no selected points")
    wanted === nothing || length(selected) == length(wanted) || error("unknown point id")
    states = Dict(String(r.point_id) => frozen_state(r, profile, ap, ac) for r in selected)
    CSV.write(joinpath(OUTPUT_DIR, "threshold_inverses.csv"),
        reduce(vcat, [eta_survey_thresholds(r, states[String(r.point_id)]) for r in selected]))
    CSV.write(joinpath(OUTPUT_DIR, "background_A.csv"), [merge(point_fields(r),
        (A_u=states[String(r.point_id)].qp.A.u, A_s=states[String(r.point_id)].qp.A.s,
         m_u_MeV=Float64(r.m_u)*HC, m_s_MeV=Float64(r.m_s)*HC)) for r in selected])
    candidates, summaries = NamedTuple[], NamedTuple[]
    for direction in ("forward", "reverse"), segment in unique(String(r.continuation_segment) for r in selected)
        scan = sort(filter(r -> String(r.continuation_segment) == segment, selected);
            by=r -> Float64(r.xi), rev=direction == "reverse")
        previous = nothing
        started = time()
        for row in scan
            cs, sm, previous = probe_roots(row, states[String(row.point_id)], :eta, direction, previous)
            append!(candidates, cs)
            push!(summaries, sm)
        end
        CSV.write(joinpath(OUTPUT_DIR, "pole_candidates.csv"), candidates)
        CSV.write(joinpath(OUTPUT_DIR, "poles.csv"), summaries)
        println("completed ", direction, " ", segment, " points=", length(scan),
            " seconds=", round(time()-started; digits=2))
        flush(stdout)
    end
    files = Dict(name => digest(joinpath(OUTPUT_DIR, name)) for name in readdir(OUTPUT_DIR))
    write_json(joinpath(OUTPUT_DIR, "summary.json"), (
        schema="transport_eta_mott_survey_v1", status="completed_diagnostic",
        timestamp_utc=string(now(UTC)), julia_version=string(VERSION),
        julia_threads=Threads.nthreads(), blas_threads=BLAS.get_num_threads(),
        hbarc_MeV_fm=HC, source_root=SOURCE_ROOT,
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        input_rows_sha256=digest(joinpath(INPUT_DIR, "source_rows.csv")),
        analysis_script_sha256=digest(@__FILE__),
        helper_script_sha256=digest(joinpath(@__DIR__, "run_transport_mott_probe.jl")),
        effective_profile=Models.TransportSamplingProfiles.effective_config(profile),
        effective_A=(p_nodes=ap, cos_nodes=ac, scale_inv_fm=profile.a_integral.scale_inv_fm),
        input_points=length(rows), root_points=length(selected),
        curves=length(unique(String(r.window) for r in selected)),
        pole_rows=length(summaries), qualified_pole_rows=count(r->r.qualified, summaries),
        ambiguous_pole_rows=count(r->r.unique_qualified_roots>1, summaries),
        root_residual_limit=ROOT_TOL,
        root_method="existing real B0 plus width-coupled pole equation",
        mixed_identity="eta algebraic minus eigenbranch; all multistart candidates retained",
        continuation="reset at reviewed first-order gaps; forward and reverse independently",
        transport_recomputed=false, output_sha256=files,
    ))
end

run_eta_survey()
