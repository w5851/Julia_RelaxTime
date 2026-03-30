using Test
using HTTP
using JSON3

const PROJECT_ROOT_TCE2E = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_TCE2E, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_TCE2E, "src", "simulation", "FullServerApp.jl"))
end

const TCE2E = Main.FullServerApp

function _body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

function _post_scan_request(payload::Dict)
    return HTTP.Request(
        "POST",
        "/api/modules/pnjl-scan/jobs",
        ["Content-Type" => "application/json"],
        JSON3.write(payload),
    )
end

function _reset_scan_jobs!()
    lock(TCE2E._PNJL_SCAN_JOBS_LOCK) do
        empty!(TCE2E._PNJL_SCAN_JOBS)
        TCE2E._PNJL_SCAN_JOB_SEQ[] = 0
        if isdefined(TCE2E, :_PNJL_SCAN_IDEMPOTENCY_CACHE)
            empty!(TCE2E._PNJL_SCAN_IDEMPOTENCY_CACHE)
        end
    end
    return nothing
end

@testset "Task center e2e smoke" begin
    _reset_scan_jobs!()

    payload = Dict(
        "kind" => "tmu",
        "params" => Dict(
            "T_values" => [150.0],
            "mu_values" => [0.0],
            "xi" => 0.0,
            "max_retries" => 0,
            "p_num" => 16,
            "t_num" => 8,
        ),
    )

    create_resp = TCE2E.handle_pnjl_scan_job_create(_post_scan_request(payload))
    create_body = _body(create_resp)
    @test create_resp.status == 202
    @test create_body.status == "accepted"
    @test haskey(create_body, :job_id)

    job_id = String(create_body.job_id)
    terminal_status = ""
    deadline = time() + 180.0
    while time() < deadline
        status_resp = TCE2E.handle_pnjl_scan_job_status(job_id)
        status_body = _body(status_resp)
        @test status_resp.status == 200
        terminal_status = String(status_body.job_status)
        if terminal_status in ("succeeded", "failed", "cancelled")
            break
        end
        sleep(0.5)
    end

    @test terminal_status == "succeeded"

    result_resp = TCE2E.handle_pnjl_scan_job_result(job_id)
    result_body = _body(result_resp)
    @test result_resp.status == 200
    @test result_body.status == "ok"
    @test result_body.job_status == "succeeded"
    @test haskey(result_body, :result)
    @test haskey(result_body.result, :output_path)
    @test isfile(String(result_body.result.output_path))

    _reset_scan_jobs!()
end
