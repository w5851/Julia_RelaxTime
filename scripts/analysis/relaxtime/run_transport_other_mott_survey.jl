#!/usr/bin/env julia
# SOURCE_ROOT INPUT_DIR OUTPUT_DIR. All numerical work runs on quote-remote.
include(joinpath(@__DIR__,"run_transport_mott_probe.jl"))
include(joinpath(@__DIR__,"transport_mott_thresholds.jl"))

function run_other_survey()
    profile=Models.TransportSamplingProfiles.load_profile(joinpath(INPUT_DIR,"base.toml"))
    rows=collect(CSV.File(joinpath(INPUT_DIR,"source_rows.csv")))
    states=Dict(String(r.point_id)=>frozen_state(r,profile,32,16) for r in rows)
    thresholds=NamedTuple[]
    for row in rows, meson in MESONS
        state=states[String(row.point_id)]; qp=state.qp
        cuts=meson in MIXED ? ((:uu,2qp.m.u),(:ss,2qp.m.s)) : ((:pair,MT.mott_threshold_mass(meson,qp)),)
        for (cut,threshold) in cuts
            for evaluator in ("diagnostic_endpoint","original_offset_1e-6","original_offset_1e-4")
                offset=evaluator=="original_offset_1e-4" ? -1e-4 : -1e-6
                value = evaluator=="diagnostic_endpoint" ? diagnostic_threshold_inverse(state,meson,cut) :
                    MP.meson_mass_equation(meson,threshold+offset,0.0,0.0,state.qp,state.tp,state.kc)
                push!(thresholds,merge(point_fields(row),(meson=String(meson),cut=String(cut),
                    continuation_segment=String(row.continuation_segment),evaluator=evaluator,
                    actual_mass_offset_inv_fm=evaluator=="diagnostic_endpoint"&&meson!=:K ? 0.0 : offset,
                    threshold_MeV=threshold*HC,inverse_re=real(value),inverse_im=imag(value),
                    inverse_abs=abs(value),inverse_units=meson in MIXED ? "fm2" : "dimensionless")))
            end
        end
    end
    CSV.write(joinpath(OUTPUT_DIR,"thresholds.csv"),thresholds)
    println("THRESHOLDS_DONE rows=",length(thresholds));flush(stdout)
    all_summaries=NamedTuple[];candidate_count=0
    mkpath(joinpath(OUTPUT_DIR,"candidates"))
    for meson in MESONS
        summaries=NamedTuple[]
        for direction in ("forward","reverse"), (si,segment) in enumerate(unique(String(r.continuation_segment) for r in rows))
            scan=sort(filter(r->String(r.continuation_segment)==segment,rows);by=r->Float64(r.xi),rev=direction=="reverse")
            previous=nothing;candidates=NamedTuple[];started=time()
            for row in scan
                cs,sm,previous=probe_roots(row,states[String(row.point_id)],meson,direction,previous)
                append!(candidates,cs);push!(summaries,sm)
            end
            CSV.write(joinpath(OUTPUT_DIR,"candidates",string(meson)*"_"*direction*"_"*lpad(si,2,'0')*".csv"),candidates)
            candidate_count+=length(candidates)
            println("POLES_DONE ",meson," ",direction," ",segment," seconds=",round(time()-started;digits=2));flush(stdout)
        end
        CSV.write(joinpath(OUTPUT_DIR,"poles_"*string(meson)*".csv"),summaries)
        append!(all_summaries,summaries)
    end
    CSV.write(joinpath(OUTPUT_DIR,"poles.csv"),all_summaries)
    files=Dict(replace(relpath(joinpath(root,f),OUTPUT_DIR),'\\'=>'/')=>digest(joinpath(root,f)) for (root,_,fs) in walkdir(OUTPUT_DIR) for f in fs)
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(schema="transport_other_mott_survey_v1",
        status="completed_diagnostic",calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        source_root=SOURCE_ROOT,input_rows_sha256=digest(joinpath(INPUT_DIR,"source_rows.csv")),
        analysis_script_sha256=digest(@__FILE__),helper_sha256=digest(joinpath(@__DIR__,"run_transport_mott_probe.jl")),
        threshold_script_sha256=digest(joinpath(@__DIR__,"transport_mott_thresholds.jl")),
        points=length(rows),mesons=string.(MESONS),pole_rows=length(all_summaries),candidate_rows=candidate_count,
        effective_profile=Models.TransportSamplingProfiles.effective_config(profile),
        root_residual_limit=ROOT_TOL,qualified_pole_rows=count(r->r.qualified,all_summaries),
        ambiguous_pole_rows=count(r->r.unique_qualified_roots>1,all_summaries),
        continuation_segments=20,production_kernel_changed=false,output_sha256=files))
end

run_other_survey()
