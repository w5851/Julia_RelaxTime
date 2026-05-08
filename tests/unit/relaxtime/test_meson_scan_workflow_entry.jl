using Test

const _GAP_MESON_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_gap_meson_mass_scan.jl")
const _MESON_DENSITY_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_meson_density_scan.jl")
const _STRICT_BW_MESON_DENSITY_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_strict_bw_meson_density_scan.jl")
const _PHASE_SHIFT_MESON_DENSITY_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_phase_shift_meson_density_scan.jl")
const _CROSSOVER_MESON_DENSITY_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_crossover_meson_density_scan.jl")
const _FREEZEOUT_MESON_MASS_PATH_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_freezeout_meson_mass_scan.jl")
const _ISENTROPIC_MESON_MASS_PATH_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_isentropic_meson_mass_scan.jl")
const _MOTT_PHASE_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_mott_phase_scan.jl")

@testset "meson scan scripts use Models entrypoint" begin
    gap_scan = read(_GAP_MESON_SCAN_SCRIPT, String)
    density_scan = read(_MESON_DENSITY_SCAN_SCRIPT, String)
    strict_bw_density_scan = read(_STRICT_BW_MESON_DENSITY_SCAN_SCRIPT, String)
    phase_shift_density_scan = read(_PHASE_SHIFT_MESON_DENSITY_SCAN_SCRIPT, String)
    crossover_density_scan = read(_CROSSOVER_MESON_DENSITY_SCAN_SCRIPT, String)
    freezeout_meson_mass_path_scan = read(_FREEZEOUT_MESON_MASS_PATH_SCAN_SCRIPT, String)
    isentropic_meson_mass_path_scan = read(_ISENTROPIC_MESON_MASS_PATH_SCAN_SCRIPT, String)
    mott_scan = read(_MOTT_PHASE_SCAN_SCRIPT, String)

    @test occursin("using .Models: solve_gap_and_meson_point", gap_scan)
    @test !occursin("using .MesonMassWorkflow: solve_gap_and_meson_point", gap_scan)

    @test occursin("using .Models: solve_gap_and_meson_density_point", density_scan)
    @test occursin("workflow_entry\" => \"Models.solve_gap_and_meson_density_point", density_scan)
    @test occursin("continuation_state = nothing", density_scan)
    @test occursin("continuation_state=continuation_state", density_scan)
    @test occursin("continuation_state = res.continuation_state", density_scan)

    @test occursin("using .Models: solve_gap_and_strict_bw_meson_density_point", strict_bw_density_scan)
    @test occursin("workflow_entry\" => \"Models.solve_gap_and_strict_bw_meson_density_point", strict_bw_density_scan)
    @test occursin("continuation_state = nothing", strict_bw_density_scan)
    @test occursin("continuation_state=continuation_state", strict_bw_density_scan)
    @test occursin("continuation_state = res.continuation_state", strict_bw_density_scan)

    @test occursin("using .Models: solve_gap_and_phase_shift_meson_density_point", phase_shift_density_scan)
    @test occursin("workflow_entry\" => \"Models.solve_gap_and_phase_shift_meson_density_point", phase_shift_density_scan)
    @test occursin("continuation_state = nothing", phase_shift_density_scan)
    @test occursin("continuation_state=continuation_state", phase_shift_density_scan)
    @test occursin("continuation_state = res.continuation_state", phase_shift_density_scan)
    @test occursin("current Phase-E3 phase-shift meson-density scan only supports xi = 0", phase_shift_density_scan)

    @test occursin("using .Models: run_crossover_meson_density_scan", crossover_density_scan)
    @test occursin("println(\"workflow_entry=\$(result.workflow_entry)\")", crossover_density_scan)

    @test occursin("using .Models: run_freezeout_meson_mass_scan", freezeout_meson_mass_path_scan)
    @test occursin("default_baseline_freezeout_xi0_loggrid_1to200_n30", freezeout_meson_mass_path_scan)
    @test occursin("build_log10_sqrts_grid(1.0, 200.0, 30)", freezeout_meson_mass_path_scan)

    @test occursin("using .Models: run_isentropic_meson_mass_scan", isentropic_meson_mass_path_scan)
    @test occursin("missing required option: --sigma-target", isentropic_meson_mass_path_scan)

    @test occursin("using .Models: solve_gap_and_meson_point", mott_scan)
    @test !occursin("using .MesonMassWorkflow: solve_gap_and_meson_point", mott_scan)
    @test occursin("continuation_state = nothing", mott_scan)
    @test occursin("continuation_state=continuation_state", mott_scan)
end
