using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_phase_shift_meson_density_scan.jl")
const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "test_outputs")
const OUTFILE = joinpath(OUTDIR, "phase_shift_meson_density_scan_smoke.csv")

@testset "phase-shift meson density scan smoke" begin
    mkpath(OUTDIR)
    isfile(OUTFILE) && rm(OUTFILE)

    cmd = `julia --project=. $SCRIPT --output $OUTFILE --overwrite --tmin 210 --tmax 210 --tstep 2 --q-nodes 8 --omega-nodes 8 --qmax 12 --omega-max 10 --p-num 8 --t-num 4 --max-iter 20`
    run(cmd)

    @test isfile(OUTFILE)
    text = read(OUTFILE, String)
    @test occursin("# workflow_entry: Models.solve_gap_and_phase_shift_meson_density_point", text)
    @test occursin("T_MeV,muB_MeV,xi", text)
    @test occursin("n_pi", text)
    @test occursin("pi_q_integral_estimate", text)

    data_lines = [
        line for line in split(text, '\n')
        if !isempty(strip(line)) && !startswith(strip(line), "#") && !startswith(strip(line), "T_MeV,")
    ]
    @test length(data_lines) == 1

    cols = split(data_lines[1], ',')
    @test length(cols) == 29
    @test cols[19] == "phase_shift_current"
end

@testset "phase-shift meson density scan smoke (gbu reference)" begin
    mkpath(OUTDIR)
    isfile(OUTFILE) && rm(OUTFILE)

    cmd = `julia --project=. $SCRIPT --output $OUTFILE --overwrite --scheme gbu --tmin 210 --tmax 210 --tstep 2 --q-nodes 8 --omega-nodes 8 --qmax 12 --omega-max 10 --p-num 8 --t-num 4 --max-iter 20`
    run(cmd)

    @test isfile(OUTFILE)
    text = read(OUTFILE, String)
    @test occursin("# phase_shift_scheme: gbu_reference", text)

    data_lines = [
        line for line in split(text, '\n')
        if !isempty(strip(line)) && !startswith(strip(line), "#") && !startswith(strip(line), "T_MeV,")
    ]
    @test length(data_lines) == 1

    cols = split(data_lines[1], ',')
    @test length(cols) == 29
    @test cols[19] == "phase_shift_gbu_reference"
end
