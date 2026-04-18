using Test

const _SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_gap_transport_scan.jl")
if !isdefined(Main, :ScanOptions)
    include(_SCAN_SCRIPT)
end

@testset "run_gap_transport_scan solver entry" begin
    @test isdefined(Main, :_solve_fixedmu_via_models_solve)
end

@testset "prefer Models.solve over solve_constraint" begin
    opts = Main.ScanOptions(
        "tmp.csv",
        nothing,
        nothing,
        nothing,
        [-0.5],
        150.0,
        150.0,
        1.0,
        0.0,
        0.0,
        1.0,
        true,
        false,
        false,
        12,
        6,
        40,
        28,
        8,
        8,
        6,
        true,
        0.6,
        8,
        10,
        :linear,
        128,
        :finite_15,
        5,
        24,
        8.0,
    )

    eq = Main.solve_models_equilibrium(
        150 / Main.Constants_PNJL.ħc_MeV_fm,
        0.0,
        -0.5,
        Main.TransportWorkflow.PNJL.HADRON_SEED_5,
        opts,
    )
    @test eq !== nothing
    @test eq.solver_backend == :models
end

@testset "failed point sidecar" begin
    opts = Main.parse_args([
        "--output", "tmp.csv",
        "--failed-points-output", "tmp_failed.csv",
        "--tmin", "150",
        "--tmax", "150",
        "--tstep", "1",
        "--mubmin", "0",
        "--mubmax", "0",
        "--mubstep", "1",
        "--xi-list", "-0.5",
    ])
    @test hasproperty(opts, :failed_points_output)
    @test opts.failed_points_output == "tmp_failed.csv"

    io = IOBuffer()
    diag = (seed_source="seed", phase_prev=:unknown, phase_curr=:hadron)
    Main.write_failed_point_row!(io, 150.0, 0.0, -0.5, diag, ArgumentError("boom"))
    line = String(take!(io))

    @test occursin("150.0,0.0,-0.5", line)
    @test occursin("seed", line)
    @test occursin("ArgumentError", line)
    @test occursin("boom", line)
end
