using Test
using HTTP
using JSON3
using Logging: NullLogger, with_logger
using Dates: DateTime, Second, Day

const PROJECT_ROOT_STJ = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_STJ, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_STJ, "src", "simulation", "FullServerApp.jl"))
end

const STJ = Main.FullServerApp

function _script_task_req(method::String, path::String, payload=Dict())
    body = method == "POST" ? JSON3.write(payload) : UInt8[]
    return HTTP.Request(method, path, ["Content-Type" => "application/json"], body)
end

function _script_task_body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

@testset "FullServer script task catalog" begin
    @test length(STJ.SCRIPT_TASK_CATALOG) >= 17

    catalog = STJ._public_script_task_catalog()
    ids = Set(String(task["id"]) for task in catalog["tasks"])
    @test "run-unified-scan" in ids
    @test "run-gap-transport-scan" in ids
    @test "run-gap-meson-mass-scan" in ids
    @test "run-mott-phase-scan" in ids
    @test "run-offline-transport-patch" in ids
    @test catalog["default_allowed_tier"] == "smoke"
    @test catalog["heavy_requires_confirmation"] == true
    @test haskey(catalog["non_default_frontend_policy"], "scripts/dev")

    phase_plots = STJ._script_task_by_id("run-phase-guided-transport-plots")
    @test phase_plots["default_preset"] == "custom"
    @test collect(keys(phase_plots["presets"])) == ["custom"]
    @test get(phase_plots, "requires_existing_input", false) == true

    resp = STJ.handle_script_task_catalog(_script_task_req("GET", "/api/modules/script-tasks"))
    body = _script_task_body(resp)
    @test resp.status == 200
    @test body.status == "ok"
    @test length(body.tasks) >= 17
end

@testset "FullServer script task request contract" begin
    job_id = "unit-script-job"

    smoke = STJ._script_task_resolve_request(
        Dict{Symbol, Any}(:task_id => "run-gap-transport-scan", :preset => "smoke"),
        job_id,
    )
    @test smoke.task_id == "run-gap-transport-scan"
    @test smoke.preset == "smoke"
    @test any(arg -> occursin(job_id, arg), smoke.args)
    @test any(arg -> occursin("frontend_jobs", arg), smoke.args)

    @test_throws ArgumentError STJ._script_task_resolve_request(
        Dict{Symbol, Any}(:task_id => "run-gap-transport-scan", :preset => "canonical"),
        job_id,
    )

    canonical = STJ._script_task_resolve_request(
        Dict{Symbol, Any}(
            :task_id => "run-gap-transport-scan",
            :preset => "canonical",
            :confirm_heavy => true,
        ),
        job_id,
    )
    @test canonical.heavy == true

    custom = STJ._script_task_resolve_request(
        Dict{Symbol, Any}(
            :task_id => "run-mott-phase-derived-csv",
            :preset => "custom",
            :confirm_heavy => true,
            :custom_args => ["--in", "data/outputs/in.csv", "--out", "data/outputs/out.csv"],
        ),
        job_id,
    )
    @test custom.args == ["--in", "data/outputs/in.csv", "--out", "data/outputs/out.csv"]

    magnetic_eb = STJ._script_task_resolve_request(
        Dict{Symbol, Any}(:task_id => "run-magnetic-eb-scan", :preset => "smoke"),
        job_id,
    )
    @test any(arg -> occursin(job_id, arg), magnetic_eb.args)
    @test any(arg -> occursin("pnjl_magnetic_eb_scan.csv", arg), magnetic_eb.args)

    magnetic_stability = STJ._script_task_resolve_request(
        Dict{Symbol, Any}(:task_id => "run-magnetic-stability-scan", :preset => "smoke"),
        job_id,
    )
    @test any(arg -> occursin("pnjl_magnetic_stability_scan.csv", arg), magnetic_stability.args)
    @test any(arg -> occursin("pnjl_magnetic_stability_failures.csv", arg), magnetic_stability.args)

    @test_throws ArgumentError STJ._script_task_resolve_request(
        Dict{Symbol, Any}(
            :task_id => "run-mott-phase-derived-csv",
            :preset => "custom",
            :confirm_heavy => true,
            :custom_args => ["--bad\narg"],
        ),
        job_id,
    )
end

@testset "concurrent script reservations respect the pending limit" begin
    lock(STJ._SCRIPT_TASK_JOBS_LOCK) do
        empty!(STJ._SCRIPT_TASK_JOBS)
    end
    request = Dict{Symbol,Any}(:task_id => "run-unified-scan", :preset => "smoke")
    try
        attempts = STJ._SCRIPT_TASK_MAX_PENDING + 8
        results = with_logger(NullLogger()) do
            fetch.([Threads.@spawn STJ._reserve_script_task_job(request) for _ in 1:attempts])
        end
        @test count(!isnothing, results) == STJ._SCRIPT_TASK_MAX_PENDING
        @test count(isnothing, results) == 8
        @test STJ._script_task_queue_snapshot().queued == STJ._SCRIPT_TASK_MAX_PENDING
        @test STJ._script_task_queue_snapshot().running == 0
    finally
        lock(STJ._SCRIPT_TASK_JOBS_LOCK) do
            for job in values(STJ._SCRIPT_TASK_JOBS)
                path = abspath(job["output_dir"])
                dirname(path) == abspath(STJ._SCRIPT_TASK_JOB_ROOT) || error("unexpected test output directory")
                rm(path; recursive=true, force=true)
            end
            empty!(STJ._SCRIPT_TASK_JOBS)
        end
    end
end

@testset "FullServer script task HTTP guards" begin
    req = _script_task_req("POST", "/api/modules/script-tasks/jobs", Dict(
        "task_id" => "run-gap-transport-scan",
        "preset" => "canonical",
    ))
    resp = STJ.handle_script_task_job_create(req)
    body = _script_task_body(resp)
    @test resp.status == 400
    @test body.status == "error"
    @test body.error_code == "INVALID_REQUEST"
    @test occursin("confirm_heavy", String(body.error))

    missing_status = STJ.handle_script_task_job_status("missing-script-job")
    missing_body = _script_task_body(missing_status)
    @test missing_status.status == 404
    @test missing_body.error_code == "JOB_NOT_FOUND"

    route_resp = STJ.route_request(
        HTTP.Request("GET", "/api/modules/script-tasks", ["Content-Type" => "application/json"], UInt8[]),
        PROJECT_ROOT_STJ,
    )
    route_body = _script_task_body(route_resp)
    @test route_resp.status == 200
    @test route_body.status == "ok"

    broken_task = Dict{String, Any}(
        "id" => "unit-broken-task",
        "name" => "Broken task",
        "script" => "scripts/pnjl/missing.jl",
        "default_preset" => "smoke",
        "presets" => Dict("smoke" => Dict("label" => "smoke", "heavy" => false)),
    )
    push!(STJ.SCRIPT_TASK_CATALOG, broken_task)
    try
        @test_logs (:error, "script task create failed") begin
            broken_resp = STJ.handle_script_task_job_create(_script_task_req("POST", "/api/modules/script-tasks/jobs", Dict(
                "task_id" => "unit-broken-task",
                "preset" => "smoke",
            )))
            broken_body = _script_task_body(broken_resp)
            @test broken_resp.status == 500
            @test broken_body.error_code == "COMPUTATION_ERROR"
            @test broken_body.error == "script task create failed"
        end
    finally
        filter!(task -> String(get(task, "id", "")) != "unit-broken-task", STJ.SCRIPT_TASK_CATALOG)
    end
end

@testset "FullServer script task command preview" begin
    task = STJ._script_task_by_id("run-unified-scan")
    args = ["scan", "tmu", "--output_path=data/outputs/frontend_jobs/unit/out.csv"]
    preview = STJ._script_task_command_preview(task, args)
    parts, _wrapper = STJ._script_task_command_parts(task, args)
    cmd = Cmd(Cmd(parts); dir=STJ._PROJECT_ROOT)
    @test occursin("run_unified_scan.jl", preview["display"])
    @test preview["wrapper"] in ("run_with_sysimage.ps1", "run_with_sysimage.sh", "julia --project fallback")
    @test preview["argv"] isa Vector
    @test cmd isa Cmd
end

@testset "FullServer script task log tail" begin
    mktempdir() do dir
        path = joinpath(dir, "long.log")
        write(path, repeat("a", 4096) * "tail-window")
        tail = STJ._tail_file(path, 64)
        @test !occursin(repeat("a", 100), tail)
        @test occursin("tail-window", tail)
        @test length(tail) <= 64
    end
end

@testset "script task terminal history retention" begin
    now = DateTime(2026, 11, 13, 12)
    fixture(id, status; ended_at=now, seq=0, worker_active=false, process=nothing) = Dict{String,Any}(
        "job_id" => id, "status" => status, "ended_at" => ended_at, "seq" => seq,
        "worker_active" => worker_active, "process" => process,
    )
    lock(STJ._SCRIPT_TASK_JOBS_LOCK) do
        saved = copy(STJ._SCRIPT_TASK_JOBS)
        empty!(STJ._SCRIPT_TASK_JOBS)
        try
            @test STJ._SCRIPT_TASK_FINISHED_KEEP_MAX == 64
            @test STJ._SCRIPT_TASK_FINISHED_TTL_SECONDS == 86400
            @test_throws ArgumentError STJ._prune_finished_script_tasks_locked!(-1; now)
            @test_throws ArgumentError STJ._prune_finished_script_tasks_locked!(1; now, ttl_seconds=-1)

            STJ._SCRIPT_TASK_JOBS["expired"] = fixture("expired", "succeeded"; ended_at=string(now - Second(11)))
            STJ._SCRIPT_TASK_JOBS["boundary"] = fixture("boundary", "failed"; ended_at=now - Second(10))
            STJ._SCRIPT_TASK_JOBS["queued"] = fixture("queued", "queued"; ended_at=now - Day(2))
            STJ._SCRIPT_TASK_JOBS["running"] = fixture("running", "running"; ended_at=now - Day(2))
            STJ._SCRIPT_TASK_JOBS["cancelling"] = fixture("cancelling", "cancelled"; ended_at=now - Day(2), worker_active=true)
            STJ._SCRIPT_TASK_JOBS["process-live"] = fixture("process-live", "cancelled"; ended_at=now - Day(2), process=:occupied)
            @test STJ._prune_finished_script_tasks_locked!(10; now, ttl_seconds=10) == 1
            @test !haskey(STJ._SCRIPT_TASK_JOBS, "expired")
            @test Set(keys(STJ._SCRIPT_TASK_JOBS)) == Set(("boundary", "queued", "running", "cancelling", "process-live"))

            # Retain most recently completed tasks, regardless of submission order.
            STJ._SCRIPT_TASK_JOBS["recent"] = fixture("recent", "succeeded"; seq=1)
            STJ._SCRIPT_TASK_JOBS["older"] = fixture("older", "cancelled"; ended_at=now - Second(2), seq=100)
            @test STJ._prune_finished_script_tasks_locked!(1; now, ttl_seconds=10) == 2
            @test haskey(STJ._SCRIPT_TASK_JOBS, "recent")
            @test !haskey(STJ._SCRIPT_TASK_JOBS, "older")
            STJ._SCRIPT_TASK_JOBS["cancelling"]["worker_active"] = false
            STJ._SCRIPT_TASK_JOBS["process-live"]["process"] = nothing
            @test STJ._prune_finished_script_tasks_locked!(0; now, ttl_seconds=10) == 3
            @test Set(keys(STJ._SCRIPT_TASK_JOBS)) == Set(("queued", "running"))

            mktempdir() do output
                result = joinpath(output, "result.csv")
                log = joinpath(output, "stdout.log")
                write(result, "value\n42\n")
                write(log, "completed")
                stale = fixture("stale", "succeeded"; ended_at=STJ.Dates.now() - Day(2))
                merge!(stale, Dict("output_dir" => output, "stdout_path" => log))
                STJ._SCRIPT_TASK_JOBS["stale"] = stale
                withenv("SCRIPT_TASK_FINISHED_KEEP_MAX" => "2", "SCRIPT_TASK_FINISHED_TTL_SECONDS" => "5") do
                    @test STJ._script_task_history_policy() == (keep_max=2, ttl_seconds=5)
                    for handler in (STJ.handle_script_task_job_status, STJ.handle_script_task_job_result, STJ.handle_script_task_job_cancel)
                        STJ._SCRIPT_TASK_JOBS["stale"] = copy(stale)
                        response = handler("stale")
                        @test response.status == 404
                        @test _script_task_body(response).error_code == "JOB_NOT_FOUND"
                        @test !haskey(STJ._SCRIPT_TASK_JOBS, "stale")
                    end
                end
                @test read(result, String) == "value\n42\n"
                @test read(log, String) == "completed"
            end
        finally
            empty!(STJ._SCRIPT_TASK_JOBS)
            merge!(STJ._SCRIPT_TASK_JOBS, saved)
        end
    end
end
