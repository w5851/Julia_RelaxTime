# ScanResultFinalize 模块单元测试

using Test
using StaticArrays
using CSV

const PROJECT_ROOT_SRF = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_SRF, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ScanResultFinalize 通过 Models 公共接口访问 SolverResult
using Main.Models.ScanResultFinalize: is_success, promote_near_converged, scan_acceptance_diagnostics

const PNJL_SRF = Models.pnjl_module()
const SolverResult_srf = getproperty(PNJL_SRF, :SolverResult)
const FixedMu_srf = getproperty(PNJL_SRF, :FixedMu)

@testset "ScanResultFinalize" begin
    @testset "is_success 判定" begin
        # SolverResult 位置参数: (mode, converged, solution, x_state, mu_vec, omega, pressure, rho_norm, entropy, energy, masses, iterations, residual_norm, xi)
        result = SolverResult_srf(
            FixedMu_srf(),            # mode
            true,                     # converged
            zeros(5),                 # solution
            SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2),  # x_state
            SVector{3, Float64}(0.0, 0.0, 0.0),                   # mu_vec
            0.0,                      # omega
            0.0,                      # pressure
            0.0,                      # rho_norm
            0.0,                      # entropy
            0.0,                      # energy
            SVector{3, Float64}(0.3, 0.31, 0.5),                  # masses
            5,                        # iterations
            1e-12,                    # residual_norm
            0.0                       # xi
        )
        @test is_success(result; acceptable_residual=1e-6)
    end

    @testset "未收敛结果" begin
        result = SolverResult_srf(
            FixedMu_srf(),
            false,
            zeros(5),
            SVector{5, Float64}(0.0, 0.0, 0.0, 0.0, 0.0),
            SVector{3, Float64}(0.0, 0.0, 0.0),
            0.0, 0.0, 0.0, 0.0, 0.0,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            100,
            1.0,
            0.0
        )
        @test !is_success(result; acceptable_residual=1e-6)
    end
end

function _scan_diagnostic_fixture(; converged=false, residual=5e-7,
    x_state=SVector(-0.03, -0.03, -0.04, 0.2, 0.2),
    convergence=Models.SolverConvergence())
    return SolverResult_srf(FixedMu_srf(), converged, collect(x_state), x_state,
        SVector(0.0, 0.0, 0.0), 0.0, 0.0, 0.25, 0.0, 0.0,
        SVector(0.3, 0.31, 0.5), 5, residual, 0.0; convergence)
end

@testset "Convergence evidence and scan decisions" begin
    legacy = _scan_diagnostic_fixture(converged=true)
    view = Models.solver_convergence_view(legacy)
    @test view.solver_converged === nothing
    @test view.residual_ok === nothing
    @test view.physical_valid === nothing
    @test view.residual_norm_max === nothing
    @test view.source === :not_recorded
    @test !view.scan_promoted
    @test Models.solver_convergence_view((convergence=nothing,)) == view
    @test scan_acceptance_diagnostics(legacy; acceptable_residual=1e-6).acceptance_reason === :legacy_converged

    evidence = Models.SolverConvergence(solver_converged=false, residual_ok=false,
        physical_valid=true, residual_norm_max=1e-8, source=:joint_nlsolve)
    near = _scan_diagnostic_fixture(convergence=evidence)
    diagnostic = scan_acceptance_diagnostics(near; acceptable_residual=1e-6)
    @test diagnostic.solver_converged === false
    @test diagnostic.solver_residual_ok === false
    @test diagnostic.scan_accepted
    @test diagnostic.acceptance_reason === :near_residual
    promoted, message = promote_near_converged(near; acceptable_residual=1e-6)
    @test promoted.converged && !isempty(message)
    @test promoted.solution == near.solution
    @test promoted.residual_norm == near.residual_norm
    @test promoted.convergence.scan_promoted
    @test promoted.convergence.solver_converged === false
    @test !promoted.convergence.residual_ok
    @test scan_acceptance_diagnostics(promoted; acceptable_residual=1e-6).acceptance_reason === :near_residual
    @test promote_near_converged(promoted; acceptable_residual=1e-6) == (promoted, "")

    raw = Models.solver_result_view(promoted)
    roundtrip = Models.coerce_solver_result(FixedMu_srf(), raw)
    @test Models.solver_convergence_view(roundtrip) == raw.convergence
    old_raw = (; (key => value for (key, value) in pairs(raw) if key !== :convergence)...)
    @test Models.solver_convergence_view(Models.coerce_solver_result(FixedMu_srf(), old_raw)) == view
    converted = Models.TmuScan._fixedmu_scan_units(promoted)
    @test Models.solver_convergence_view(converted) == raw.convergence
    @test converted.rho_norm ≈ promoted.rho_norm * Main.Constants_PNJL.ρ0_inv_fm3

    good = _scan_diagnostic_fixture(converged=true, residual=1e-10,
        convergence=Models.SolverConvergence(solver_converged=true, residual_ok=true,
            physical_valid=true, residual_norm_max=1e-6, source=:gap_nlsolve))
    @test scan_acceptance_diagnostics(good; acceptable_residual=1e-6).acceptance_reason === :solver_converged
    @test scan_acceptance_diagnostics(good; acceptable_residual=1e-6,
        accepted=false).acceptance_reason === :scan_rejected
    for (residual, reason) in ((1.0, :residual_rejected), (Inf, :nonfinite_residual), (NaN, :nonfinite_residual))
        rejected = _scan_diagnostic_fixture(; residual)
        diagnostic = scan_acceptance_diagnostics(rejected; acceptable_residual=1e-6)
        @test !diagnostic.scan_accepted
        @test diagnostic.acceptance_reason === reason
        @test first(promote_near_converged(rejected; acceptable_residual=1e-6)) === rejected
    end
    absent = scan_acceptance_diagnostics(nothing; acceptable_residual=1e-6)
    @test !absent.scan_accepted
    @test absent.solver_converged === nothing
    @test absent.acceptance_reason === :no_result

    unphysical = _scan_diagnostic_fixture(converged=true,
        x_state=SVector(-0.03, -0.03, -0.04, 1.2, 0.2))
    @test !Models.TmuScan._is_success(unphysical)
    io = IOBuffer()
    println(io, Models.TmuScan.HEADER)
    Models.TmuScan._write_row(io, 150.0, 0.0, 0.0, unphysical, "nonphysical")
    row = only(CSV.File(take!(io)))
    @test row.converged
    @test !row.scan_accepted
    @test !row.physical_valid
    @test row.acceptance_reason == "physical_invalid"
end

@testset "Scan CSV diagnostics and legacy resume" begin
    common = Models.TmuScan.ScanCommon
    for scan in (Models.TmuScan, Models.TrhoScan)
        mktempdir() do directory
            file = joinpath(directory, "scan.csv")
            result = _scan_diagnostic_fixture()
            # Generate an actual legacy row, including a quoted comma in its message.
            io = IOBuffer()
            scan._write_row(io, 150.0, 0.25, 0.0, result, "seed failed, recovered")
            current_row = chomp(String(take!(io)))
            legacy_row = first(rsplit(current_row, ','; limit=length(common.CONVERGENCE_COLUMNS) + 1))
            write(file, scan.LEGACY_HEADER * "\n" * legacy_row * "\n\n")
            original = only(CSV.File(file))
            common.ensure_convergence_header!(file, scan.LEGACY_HEADER, scan.HEADER)
            upgraded = only(CSV.File(file))
            for name in propertynames(original)
                @test isequal(getproperty(upgraded, name), getproperty(original, name))
            end
            @test upgraded.solver_converged == "unknown"
            @test upgraded.physical_valid == "unknown"
            @test upgraded.scan_accepted == "unknown"
            @test upgraded.acceptance_reason == "not_recorded"
            @test length(propertynames(upgraded)) == length(split(scan.HEADER, ','))
            bytes = read(file)
            common.ensure_convergence_header!(file, scan.LEGACY_HEADER, scan.HEADER)
            @test read(file) == bytes
            open(file, "a") do output
                scan._write_row(output, 151.0, 0.25, 0.0, nothing, "all attempts failed")
            end
            rows = collect(CSV.File(file))
            @test length(rows) == 2
            @test string(rows[2].scan_accepted) == "false"
            @test rows[2].acceptance_reason == "no_result"
            @test rows[2].convergence_source == "not_recorded"

            write(file, "unsupported,header\n1,2\n")
            bytes = read(file)
            @test_throws ArgumentError common.ensure_convergence_header!(file, scan.LEGACY_HEADER, scan.HEADER)
            @test read(file) == bytes
            @test readdir(directory) == ["scan.csv"]
        end
    end
end
