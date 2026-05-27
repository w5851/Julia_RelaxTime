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
