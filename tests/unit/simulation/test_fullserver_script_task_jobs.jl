using Test
using HTTP
using JSON3

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
