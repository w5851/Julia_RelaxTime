using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "Scan cancellation preserves flushed rows for resume" begin
    for (scan, grid) in ((Models.run_tmu_scan, (; mu_values=[0.0, 5.0, 10.0])),
                        (Models.run_trho_scan, (; rho_values=[0.2, 0.3, 0.4], reverse_rho=false)))
        mktempdir() do dir
            output = joinpath(dir, "partial.csv")
            common = (; T_values=[150.0], xi_values=[0.0], grid...,
                continuation=:single_pass, p_num=12, t_num=4, iterations=80)
            @test_throws Models.ScanStopRequested scan(; common..., output_path=output,
                overwrite=true, stop_requested=() -> true)
            @test !ispath(output)

            points = Float64[]
            @test_throws Models.ScanStopRequested scan(; common..., output_path=output,
                overwrite=true, stop_requested=() -> !isempty(points),
                progress_cb=(point, _) -> push!(points, hasproperty(point, :mu) ? point.mu : point.rho))
            @test length(points) == 1
            partial = readlines(output)
            @test length(partial) == 2

            stats = scan(; common..., output_path=output, resume=true, overwrite=false)
            @test (stats.skipped, stats.success, stats.failure) == (1, 2, 0)
            rows = readlines(output)
            @test rows[1:2] == partial
            @test [parse(Float64, split(row, ',')[2]) for row in rows[2:end]] ==
                (haskey(grid, :mu_values) ? grid.mu_values : grid.rho_values)
        end
    end
end

@testset "Continuation propagates stop requests" begin
    model = Models.create_model(:PNJL)
    hbarc = Main.Constants_PNJL.ħc_MeV_fm
    mu_path = Models.FixedMuPath(150 / hbarc, [0., 20., 40.] ./ hbarc; p_num=12, t_num=4)
    rho_path = Models.FixedRhoPath(150 / hbarc, [0.2, 0.3, 0.4]; p_num=12, t_num=4)
    palc = Models.PALCWorkTelemetry()
    @test_throws Models.ScanStopRequested Models.solve_path(model, mu_path;
        continuation_strategy=Models.PALCContinuation(), path_work=palc,
        stop_requested=() -> palc.continuation_residual_calls > 0)
    @test palc.continuation_points == 0  # stopped inside BifurcationKit before a curve returned

    natural = Models.PALCWorkTelemetry()
    @test_throws Models.ScanStopRequested Models.solve_path(model, rho_path;
        continuation_strategy=Models.NaturalContinuation(), path_work=natural,
        stop_requested=() -> natural.anchor_solves > 0)
    @test natural.continuation_residual_calls == 0  # no PALC recovery on cancellation
end
