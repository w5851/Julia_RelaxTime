using Test
using HTTP
using JSON3

const PROJECT_ROOT_SMRC = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_SMRC, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_SMRC, "src", "simulation", "FullServerApp.jl"))
end

const SMRC = Main.FullServerApp

@testset "service module registry contract" begin
    resp = SMRC.handle_modules_list()
    @test resp.status == 200

    body = JSON3.read(String(resp.body))
    @test length(body) >= 3

    modules_by_id = Dict(String(item.id) => item for item in body)
    @test haskey(modules_by_id, "pnjl-gap")
    @test haskey(modules_by_id, "pnjl-scan")
    @test haskey(modules_by_id, "transport-point")

    gap = modules_by_id["pnjl-gap"]
    @test gap.invocation_style == "sync"
    @test gap.service_surface == "point"
    @test gap.default_client_surface == "service"
    @test gap.stable_entrypoint == "Models.solve_pnjl_point"
    @test gap.http.method == "POST"
    @test gap.http.path == "/api/modules/pnjl-gap/run"

    scan = modules_by_id["pnjl-scan"]
    @test scan.invocation_style == "async"
    @test scan.service_surface == "job"
    @test scan.default_client_surface == "service"
    @test scan.stable_entrypoint == "Models.run_scan_pipeline"
    @test scan.http.create.path == "/api/modules/pnjl-scan/jobs"
    @test scan.http.status.path == "/api/modules/pnjl-scan/jobs/{job_id}"
    @test scan.http.result.path == "/api/modules/pnjl-scan/jobs/{job_id}/result"
    @test scan.http.cancel.path == "/api/modules/pnjl-scan/jobs/{job_id}/cancel"

    transport = modules_by_id["transport-point"]
    @test transport.invocation_style == "sync"
    @test transport.service_surface == "point"
    @test transport.default_client_surface == "service"
    @test transport.stable_entrypoint == "Models.solve_transport_from_equilibrium"
    @test transport.http.method == "POST"
    @test transport.http.path == "/api/modules/transport-point/run"
end
