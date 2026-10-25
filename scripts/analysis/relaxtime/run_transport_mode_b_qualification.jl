#!/usr/bin/env julia
# Apply the existing independent protocols to explicitly selected mode-B points.
using JSON, SHA

function run_mode_b_qualification(inputdir,outdir,stage;qualified_dir=nothing)
    selection=JSON.parsefile(joinpath(inputdir,"selection.json"))
    selection["mode"]=="B" || throw(ArgumentError("mode-B selection is required"))
    for name in ("review_points","channel_plan")
        bytes2hex(sha256(read(joinpath(inputdir,name*".csv"))))==selection[name*"_sha256"] ||
            error("mode-B qualification input changed: $name")
    end
    ids=String.(selection["points"])
    if stage in ("cheap","kernel","outer")
        cases=stage=="cheap" ? nothing : [(String(row[1]),Symbol(row[2])) for row in selection[stage*"_cases"]]
        chosen=stage=="cheap" ? ids : sort!(unique(first.(cases)))
        full_rs_layer_validation(inputdir,outdir;layer=stage,point_ids=chosen,channel_cases=cases,extra_sources=[@__FILE__]) ||
            error("mode-B $stage qualification failed")
    elseif stage in ("pilot","joint")
        chosen=stage=="pilot" ? [String(selection["pilot_point_id"])] : ids
        joint_validation(inputdir,outdir;scope="representatives",candidate_p=48,point_ids=chosen,extra_sources=[@__FILE__])
        status=JSON.parsefile(joinpath(outdir,"status.json"))
        status["joint_target_passed"] || error("mode-B joint qualification failed")
    elseif stage=="production"
        qualified_dir===nothing && throw(ArgumentError("production requires the reviewed joint summary directory"))
        production_validation(inputdir,outdir;scope="all",point_ids=ids,qualified_dir=qualified_dir,extra_sources=[@__FILE__])
    else
        throw(ArgumentError("stage must be cheap, kernel, outer, pilot, joint, or production"))
    end
    cp(joinpath(inputdir,"selection.json"),joinpath(outdir,"inputs","mode_b_selection.json"))
    manifest=JSON.parsefile(joinpath(outdir,"manifest.json"))
    push!(manifest["input_files"],Dict("path"=>joinpath(inputdir,"selection.json"),
        "sha256"=>bytes2hex(sha256(read(joinpath(inputdir,"selection.json")))),
        "snapshot"=>joinpath("inputs","mode_b_selection.json")))
    write_json(joinpath(outdir,"manifest.json"),manifest)
    return nothing
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    3<=length(ARGS)<=4 || error("Usage: run_transport_mode_b_qualification.jl INPUT_DIR NEW_OUTPUT_DIR STAGE [JOINT_SUMMARY]")
    stage=ARGS[3]
    if stage=="production"
        include(joinpath(@__DIR__,"run_transport_production_validation.jl"))
    else
        include(joinpath(@__DIR__,"run_transport_full_rs_layer_validation.jl"))
    end
    run_mode_b_qualification(abspath(ARGS[1]),abspath(ARGS[2]),stage;
        qualified_dir=length(ARGS)==4 ? abspath(ARGS[4]) : nothing)
end
