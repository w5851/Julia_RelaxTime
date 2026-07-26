using Test

const CRITICALITY_FEASIBILITY_SCRIPT = normpath(joinpath(
    @__DIR__,
    "..",
    "..",
    "..",
    "scripts",
    "analysis",
    "pnjl_criticality_feasibility.jl",
))

if !isdefined(Main, :PNJLCriticalityFeasibility)
    include(CRITICALITY_FEASIBILITY_SCRIPT)
end

const CriticalityFeasibility = Main.PNJLCriticalityFeasibility

@testset "PNJL criticality feasibility diagnostic" begin
    cases = Dict(case.name => case for case in CriticalityFeasibility.default_cases())

    visible = CriticalityFeasibility.evaluate_case(cases["resolved_visible_s"])
    @test visible.baseline_has_s_shape
    @test visible.fit_s_topology_all
    @test visible.validation_has_s_topology
    @test visible.observed_status == "resolved_s_shape"

    hidden = CriticalityFeasibility.evaluate_case(cases["hidden_s_warning"])
    @test !hidden.baseline_has_s_shape
    @test !hidden.validation_available
    @test hidden.fit_min_derivative_median < 0
    @test hidden.observed_status == "near_critical"

    validated = CriticalityFeasibility.evaluate_case(cases["hidden_s_targeted_validation"])
    @test !validated.baseline_has_s_shape
    @test validated.validation_available
    @test validated.validation_min_secant_slope < 0
    @test validated.validation_negative_secant_run >= 3
    @test validated.observed_status == "resolved_s_shape"

    critical = CriticalityFeasibility.evaluate_case(cases["critical_cusp"])
    @test abs(critical.fit_min_derivative_median) <= 1e-12
    @test critical.observed_status == "near_critical"

    unverified = CriticalityFeasibility.evaluate_case(cases["monotone_without_certificate"])
    @test !unverified.validation_available
    @test unverified.fit_min_derivative_median > 0
    @test unverified.observed_status == "unresolved"

    certified = CriticalityFeasibility.evaluate_case(cases["monotone_with_certificate"])
    @test certified.validation_available
    @test certified.validation_min_secant_slope > 0
    @test certified.observed_status == "resolved_monotone"

    noisy = CriticalityFeasibility.evaluate_case(cases["noisy_monotone_minimum_bias"])
    @test noisy.baseline_has_s_shape
    @test !noisy.fit_stable
    @test noisy.observed_status == "unresolved"

    hidden_case = cases["hidden_s_warning"]
    coarse_rho, coarse_mu = CriticalityFeasibility._case_curve(hidden_case)
    validation_x = collect(range(-1.0, 1.0; length=401))
    validation_rho = validation_x .+ 1.0
    validation_mu = validation_x .^ 3 .+ 0.08 .* validation_x
    validation_mu[201] -= 0.01
    isolated_negative = CriticalityFeasibility.analyze_curve(
        coarse_rho,
        coarse_mu;
        validation_rho=validation_rho,
        validation_mu=validation_mu,
    )
    @test isolated_negative.validation_min_secant_slope < 0
    @test isolated_negative.validation_negative_secant_run < 3
    @test !isolated_negative.validation_has_s_topology
    @test isolated_negative.status == :near_critical

    decreasing_rho = collect(range(0.0, 2.0; length=13))
    decreasing_validation_rho = collect(range(0.0, 2.0; length=401))
    decreasing = CriticalityFeasibility.analyze_curve(
        decreasing_rho,
        .-decreasing_rho;
        validation_rho=decreasing_validation_rho,
        validation_mu=.-decreasing_validation_rho,
    )
    @test decreasing.fit_min_derivative_median < 0
    @test !decreasing.fit_s_topology_all
    @test !decreasing.validation_has_s_topology
    @test decreasing.status == :unresolved

    suite = CriticalityFeasibility.run_suite(repetitions=1)
    @test all(row.matches_expected for row in suite.rows)
    @test count(row -> row.comparable, suite.performance) == 2
    @test count(row -> !row.comparable, suite.performance) == 2

    mktempdir() do output_dir
        CriticalityFeasibility.write_evidence_package(output_dir; repetitions=1)
        case_path = joinpath(output_dir, "tables", "case_results.csv")
        manifest_path = joinpath(output_dir, "manifest.json")
        @test isfile(case_path)
        @test isfile(manifest_path)
        @test !occursin(r"\b(?:NaN|Inf)\b", read(case_path, String))
        manifest = CriticalityFeasibility.JSON3.read(read(manifest_path, String))
        @test manifest.case_count == 7
        @test manifest.all_cases_match_expected
        @test length(manifest.producer_scripts) == 2
        CriticalityFeasibility._refresh_manifest(output_dir)
        @test isfile(manifest_path)
    end
end

@testset "PNJL criticality feasibility validation contract" begin
    bad_config = CriticalityFeasibility.CriticalityConfig(
        near_slope_tol=0.02,
        resolved_slope_margin=0.01,
    )
    @test_throws ArgumentError CriticalityFeasibility._validate_config(bad_config)

    rho = collect(range(0.0, 2.0; length=13))
    mu = rho .^ 3
    assessment = CriticalityFeasibility.analyze_curve(rho, mu)
    @test assessment.status != :resolved_monotone
    @test !assessment.validation_available

    @test_throws ArgumentError CriticalityFeasibility.analyze_curve(rho[1:5], mu[1:5])
    @test_throws ArgumentError CriticalityFeasibility.analyze_curve(rho, mu[1:end-1])

    too_short_rho = collect(range(0.0, 2.0; length=9))
    @test_throws ArgumentError CriticalityFeasibility.analyze_curve(
        too_short_rho,
        too_short_rho .^ 3,
    )
    @test_throws ErrorException CriticalityFeasibility._assert_finite_rows(
        [(value=Inf,)],
        "test rows",
    )
end
