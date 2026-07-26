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

    cascade_rows = Dict(
        row.case_name => row for row in
        CriticalityFeasibility.evaluate_cascade_case.(values(cases))
    )
    cascade_visible = cascade_rows["resolved_visible_s"]
    @test cascade_visible.observed_status == "resolved_s_shape"
    @test cascade_visible.stage == "fast_topology"
    @test cascade_visible.targeted_function_evaluations == 0
    @test cascade_visible.polynomial_fits == 0

    cascade_hidden = cascade_rows["hidden_s_warning"]
    @test cascade_hidden.observed_status == "resolved_s_shape"
    @test cascade_hidden.stage == "rho_support_targeted"
    @test 0 < cascade_hidden.targeted_function_evaluations <= 12
    @test cascade_hidden.polynomial_fits == 1
    @test cascade_hidden.targeted_min_secant_slope < 0
    @test cascade_hidden.spinodal_rho_gap > 0

    cascade_critical = cascade_rows["critical_cusp"]
    @test cascade_critical.observed_status == "near_critical"
    @test abs(cascade_critical.fit_min_derivative) <= 1e-12

    @test cascade_rows["monotone_without_certificate"].observed_status == "unresolved"
    @test cascade_rows["monotone_with_certificate"].observed_status == "unresolved"
    @test cascade_rows["noisy_monotone_minimum_bias"].observed_status == "unresolved"
    @test all(row.matches_expected for row in values(cascade_rows))

    strict_support = CriticalityFeasibility.RhoSupportConfig(support_slope_tol=1e-4)
    strict_rho, strict_mu = CriticalityFeasibility._case_curve(cases["hidden_s_warning"])
    strict_sampler = rho -> CriticalityFeasibility._case_mu(
        cases["hidden_s_warning"],
        rho,
    )
    strict_without_prior = CriticalityFeasibility.analyze_curve_cascade(
        strict_rho,
        strict_mu;
        sample_mu=strict_sampler,
        config=strict_support,
    )
    @test strict_without_prior.status == :unresolved
    strict_with_prior = CriticalityFeasibility.analyze_curve_cascade(
        strict_rho,
        strict_mu;
        sample_mu=strict_sampler,
        prior=CriticalityFeasibility.RhoSupportPrior(0.8, 1.0, 0.4),
        config=strict_support,
    )
    @test strict_with_prior.status == :resolved_s_shape
    @test strict_with_prior.support_origin == :temperature_prior
    @test strict_with_prior.extra_point_count <= 12

    noisy_case = cases["noisy_monotone_minimum_bias"]
    noisy_rho, noisy_mu = CriticalityFeasibility._case_curve(noisy_case)
    noisy_with_prior = CriticalityFeasibility.analyze_curve_cascade(
        noisy_rho,
        noisy_mu;
        sample_mu=rho -> CriticalityFeasibility._case_mu(noisy_case, rho),
        prior=CriticalityFeasibility.RhoSupportPrior(0.8, 1.0, 0.4),
    )
    @test noisy_with_prior.status == :unresolved
    @test noisy_with_prior.fit_rmse > 0.02

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

    decreasing_cascade = CriticalityFeasibility.analyze_curve_cascade(
        decreasing_rho,
        .-decreasing_rho;
        sample_mu=rho -> -rho,
    )
    @test decreasing_cascade.status == :unresolved
    @test decreasing_cascade.extra_point_count == 0

    sequence = CriticalityFeasibility.run_temperature_sequence()
    @test length(sequence.rows) == 10
    @test all(row.matches_expected for row in sequence.rows)
    @test maximum(row.targeted_function_evaluations for row in sequence.rows) <= 12
    @test first(sequence.rows).targeted_function_evaluations == 0
    @test all(row.prior_used for row in sequence.rows[2:end])
    resolved_widths = [
        row.spinodal_rho_gap for row in sequence.rows if !ismissing(row.spinodal_rho_gap)
    ]
    @test issorted(resolved_widths; rev=true)
    @test sequence.extrapolation.found
    @test abs(sequence.extrapolation.T_cep - 1.0) <= 0.05

    exact_T = [0.5, 0.7, 0.9]
    exact_widths = [2 * sqrt(0.12 * (1.0 - T) / 3) for T in exact_T]
    exact_cep = CriticalityFeasibility.extrapolate_cep_from_density_width(
        exact_T,
        exact_widths,
    )
    @test exact_cep.found
    @test exact_cep.T_cep ≈ 1.0 atol=1e-12

    suite = CriticalityFeasibility.run_suite(repetitions=1)
    @test all(row.matches_expected for row in suite.rows)
    @test all(row.matches_expected for row in suite.cascade_rows)
    @test all(row.matches_expected for row in suite.temperature_sequence)
    @test suite.cep_extrapolation.within_synthetic_gate
    @test count(row -> row.comparable, suite.performance) == 4
    @test count(row -> !row.comparable, suite.performance) == 2
    cascade_performance = only(filter(
        row -> row.method == "rho_support_cascade",
        suite.performance,
    ))
    @test cascade_performance.targeted_function_evaluations <= 12 * length(suite.cascade_rows)
    @test cascade_performance.classification_matches == length(suite.cascade_rows)
    multi_evidence_performance = only(filter(
        row -> row.method == "multi_evidence_criticality",
        suite.performance,
    ))
    @test multi_evidence_performance.performance_verdict == "performance_rejected"

    mktempdir() do output_dir
        CriticalityFeasibility.write_evidence_package(output_dir; repetitions=1)
        case_path = joinpath(output_dir, "tables", "case_results.csv")
        cascade_path = joinpath(output_dir, "tables", "cascade_case_results.csv")
        sequence_path = joinpath(output_dir, "tables", "temperature_sequence_results.csv")
        cep_path = joinpath(output_dir, "tables", "cep_extrapolation_summary.csv")
        manifest_path = joinpath(output_dir, "manifest.json")
        @test isfile(case_path)
        @test isfile(cascade_path)
        @test isfile(sequence_path)
        @test isfile(cep_path)
        @test isfile(manifest_path)
        @test !occursin(r"\b(?:NaN|Inf)\b", read(case_path, String))
        @test !occursin(r"\b(?:NaN|Inf)\b", read(cascade_path, String))
        @test !occursin(r"\b(?:NaN|Inf)\b", read(sequence_path, String))
        @test !occursin(r"\b(?:NaN|Inf)\b", read(cep_path, String))
        manifest = CriticalityFeasibility.JSON3.read(read(manifest_path, String))
        @test manifest.case_count == 7
        @test manifest.all_cases_match_expected
        @test manifest.cascade_case_count == 7
        @test manifest.all_cascade_cases_match_expected
        @test manifest.temperature_sequence_count == 10
        @test manifest.all_temperature_sequence_match_expected
        @test manifest.cep_extrapolation_within_synthetic_gate
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
    @test_throws ArgumentError CriticalityFeasibility._validate_rho_support_config(
        CriticalityFeasibility.RhoSupportConfig(target_point_count=4),
    )
    @test_throws ArgumentError CriticalityFeasibility._validate_rho_support_config(
        CriticalityFeasibility.RhoSupportConfig(max_extra_points=8),
    )

    rho = collect(range(0.0, 2.0; length=13))
    mu = rho .^ 3
    assessment = CriticalityFeasibility.analyze_curve(rho, mu)
    @test assessment.status != :resolved_monotone
    @test !assessment.validation_available
    @test_throws ArgumentError CriticalityFeasibility.analyze_curve_cascade(
        rho,
        mu;
        sample_mu=value -> value^3,
        prior=CriticalityFeasibility.RhoSupportPrior(0.0, 1.0, 0.0),
    )

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
    @test_throws ArgumentError CriticalityFeasibility.extrapolate_cep_from_density_width(
        [0.0, 0.5],
        [0.2],
    )
    @test_throws ArgumentError CriticalityFeasibility.extrapolate_cep_from_density_width(
        [0.0, 0.5, 0.9],
        [0.2, 0.1, 0.0],
    )
    @test_throws ArgumentError CriticalityFeasibility.extrapolate_cep_from_density_width(
        [0.5, 0.5, 0.5],
        [0.3, 0.2, 0.1],
    )
    too_few = CriticalityFeasibility.extrapolate_cep_from_density_width(
        [0.0, 0.5],
        [0.3, 0.2],
    )
    @test !too_few.found
end
