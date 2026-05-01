using Test

const _GAP_MESON_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_gap_meson_mass_scan.jl")
const _MESON_DENSITY_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_meson_density_scan.jl")
const _MOTT_PHASE_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_mott_phase_scan.jl")

@testset "meson scan scripts use Models entrypoint" begin
    gap_scan = read(_GAP_MESON_SCAN_SCRIPT, String)
    density_scan = read(_MESON_DENSITY_SCAN_SCRIPT, String)
    mott_scan = read(_MOTT_PHASE_SCAN_SCRIPT, String)

    @test occursin("using .Models: solve_gap_and_meson_point", gap_scan)
    @test !occursin("using .MesonMassWorkflow: solve_gap_and_meson_point", gap_scan)

    @test occursin("using .Models: solve_gap_and_meson_density_point", density_scan)
    @test occursin("workflow_entry\" => \"Models.solve_gap_and_meson_density_point", density_scan)
    @test occursin("equilibrium_seed_state = nothing", density_scan)
    @test occursin("meson_seed_state = nothing", density_scan)
    @test occursin("seed_state=(equilibrium_seed_state === nothing ? Models.HADRON_SEED_5 : equilibrium_seed_state)", density_scan)
    @test occursin("meson_seed_state=meson_seed_state", density_scan)

    @test occursin("using .Models: solve_gap_and_meson_point", mott_scan)
    @test !occursin("using .MesonMassWorkflow: solve_gap_and_meson_point", mott_scan)
end
