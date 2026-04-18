using Test

const _RUN_SCAN_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_scan.jl"))
if !isdefined(Main, :SingleEntryScanCLI)
    Base.include(Main, _RUN_SCAN_PATH)
end

using Main.SingleEntryScanCLI

@testset "run_scan subcommand parsing" begin
    sub, rest = SingleEntryScanCLI.parse_entry_args(["gap-transport", "--tmin", "120"])
    @test sub == :gap_transport
    @test rest == ["--tmin", "120"]

    sub2, rest2 = SingleEntryScanCLI.parse_entry_args(["tau-vs-t"])
    @test sub2 == :tau_vs_t
    @test isempty(rest2)

    sub3, rest3 = SingleEntryScanCLI.parse_entry_args(["manual-workflow", "--overwrite"])
    @test sub3 == :manual_workflow
    @test rest3 == ["--overwrite"]

    @test_throws ArgumentError SingleEntryScanCLI.parse_entry_args(String[])
    @test_throws ArgumentError SingleEntryScanCLI.parse_entry_args(["unknown-subcmd"])
end

@testset "run_scan subcommand target resolution" begin
    gap_script = SingleEntryScanCLI.resolve_target_script(:gap_transport)
    tau_script = SingleEntryScanCLI.resolve_target_script(:tau_vs_t)
    manual_script = SingleEntryScanCLI.resolve_target_script(:manual_workflow)

    @test endswith(gap_script, joinpath("scripts", "relaxtime", "run_gap_transport_scan.jl"))
    @test endswith(tau_script, joinpath("scripts", "relaxtime", "scan_relaxation_times_vs_T.jl"))
    @test endswith(manual_script, joinpath("scripts", "relaxtime", "run_manual_relaxation_scan_workflow.jl"))
end

@testset "run_scan help and manual wrapper availability" begin
    io = IOBuffer()
    SingleEntryScanCLI.print_usage(io)
    help_text = String(take!(io))
    @test occursin("gap-transport", help_text)
    @test occursin("tau-vs-t", help_text)
    @test occursin("manual-workflow", help_text)

    manual_script = SingleEntryScanCLI.resolve_target_script(:manual_workflow)
    @test isfile(manual_script)
end
