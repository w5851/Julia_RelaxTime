using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(
    PROJECT_ROOT,
    "tests", "baselines", "relaxtime", "baseline_meson_density_freezeout_phase_shift_gbu_path_v1.csv",
)

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _read_csv_rows(path::String)
    isfile(path) || error("CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    rows = NamedTuple[]
    for line in lines[2:end]
        row = strip(line)
        isempty(row) && continue
        cols = split(row, ','; keepempty=true)
        length(cols) == length(header) || error("invalid CSV row: $line")
        push!(rows, (; (Symbol(header[i]) => cols[i] for i in eachindex(header))...))
    end
    return rows
end

function _run_freezeout_phase_shift_gbu_scan(output_path::String, sqrt_s_values::Vector{Float64})
    return Models.run_freezeout_meson_density_scan(
        sqrt_s_NN_values=sqrt_s_values,
        xi_values=[0.0],
        freezeout_profile_name="default",
        path_profile_name="baseline_freezeout",
        flavor_chemical_profile_name="default",
        meson_chemical_profile_name="blaschke2019_mu_pi_100",
        regime=:phase_shift_gbu,
        output_path=output_path,
        traversal=:as_given,
        overwrite=true,
        resume=false,
        p_num=24,
        t_num=8,
        strict_bw_qmax=12.0,
        strict_bw_q_nodes=48,
        strict_bw_omega_max=10.0,
        strict_bw_omega_nodes=48,
        phase_shift_qmax=12.0,
        phase_shift_q_nodes=48,
        phase_shift_omega_min=0.05,
        phase_shift_omega_max=10.0,
        phase_shift_omega_nodes=48,
        phase_shift_eta=1e-6,
        solver_kwargs=(; iterations=40),
        mass_kwargs=(; iterations=40),
    )
end

@testset "Meson density freezeout phase_shift_gbu full-path regression" begin
    baseline_rows = _read_csv_rows(BASELINE_PATH)
    sqrt_s_values = Float64[parse(Float64, row.sqrt_s_NN_GeV) for row in baseline_rows]

    outdir = mktempdir()
    scan_path = joinpath(outdir, "freezeout_phase_shift_gbu_scan.csv")
    summary = _run_freezeout_phase_shift_gbu_scan(scan_path, sqrt_s_values)
    actual_rows = _read_csv_rows(scan_path)

    @test summary.total == length(baseline_rows)
    @test summary.success == length(baseline_rows)
    @test summary.failure == 0
    @test summary.skipped == 0
    @test length(actual_rows) == length(baseline_rows)

    rtol = 5e-6
    atol = 1e-8

    for (baseline, actual) in zip(baseline_rows, actual_rows)
        @testset "sqrt(s_NN)=$(baseline.sqrt_s_NN_GeV)" begin
            @test parse(Float64, actual.sqrt_s_NN_GeV) == parse(Float64, baseline.sqrt_s_NN_GeV)
            @test strip(actual.equilibrium_converged) == strip(baseline.equilibrium_converged)
            @test isapprox(parse(Float64, actual.T_MeV), parse(Float64, baseline.T_MeV); rtol=rtol, atol=atol)
            @test isapprox(parse(Float64, actual.muB_MeV), parse(Float64, baseline.muB_MeV); rtol=rtol, atol=atol)
            @test isapprox(parse(Float64, actual.n_pi), parse(Float64, baseline.n_pi); rtol=rtol, atol=atol)
            @test isapprox(parse(Float64, actual.n_K), parse(Float64, baseline.n_K); rtol=rtol, atol=atol)
            @test isapprox(parse(Float64, actual.kpi_ratio), parse(Float64, baseline.kpi_ratio); rtol=rtol, atol=atol)
        end
    end
end
