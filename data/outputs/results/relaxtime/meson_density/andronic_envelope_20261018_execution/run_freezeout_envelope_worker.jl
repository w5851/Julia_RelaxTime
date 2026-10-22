#!/usr/bin/env julia
# Keep the existing numerical runner loaded across sequential research batches.
include(joinpath(@__DIR__, "run_charged_gbu_freezeout_comparison.jl"))
using JSON3

function run_envelope_worker(args)
    length(args) == 2 || throw(ArgumentError("usage: EXECUTION_DIR PART(0|1)"))
    execution, part = abspath(args[1]), args[2]
    part in ("0", "1") || throw(ArgumentError("part must be 0 or 1"))
    isdir(execution) || throw(ArgumentError("execution directory does not exist"))
    runner = Main.ChargedGBUFreezeoutComparison
    sequence = 1
    while true
        stem = "worker_$(part)_$(lpad(sequence, 3, '0'))"
        request = joinpath(execution, stem * ".json")
        # The coordinator publishes each complete JSON by atomic rename.
        while !isfile(request)
            sleep(1)
        end
        job = JSON3.read(read(request, String))
        job.schema == "freezeout_envelope_worker_request_v1" ||
            throw(ArgumentError("unsupported worker request"))
        job.action == "stop" && break
        job.action == "scan" || throw(ArgumentError("action must be scan or stop"))
        manifest = runner.run_scan(; output=String(job.output),
            profiles=String.(job.profiles), energies=Float64.(job.energies_GeV),
            resolution=String(job.resolution), profile_dir=String(job.profile_dir),
            parameter_policy="literature_envelope_exploration")
        done = joinpath(execution, stem * "_done.json")
        runner.Workflow.writejson(done, (; schema="freezeout_envelope_worker_completion_v1",
            request_sha256=runner.Workflow.hashfile(request),
            output=String(job.output), scan_identity=manifest.scan_identity,
            point_count=manifest.point_count, successful_points=manifest.successful_points))
        sequence += 1
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_envelope_worker(ARGS)
end
