using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const ANALYSIS_SCRIPT = joinpath(
    PROJECT_ROOT,
    "scripts",
    "analysis",
    "pnjl_criticality_feasibility.jl",
)

function parse_cli(args)
    output = joinpath(PROJECT_ROOT, "data", "outputs", "perf", "pnjl_criticality_feasibility_probe.csv")
    repetitions = 100
    for arg in args
        if startswith(arg, "--output=")
            output = abspath(split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--repetitions=")
            repetitions = parse(Int, split(arg, "="; limit=2)[2])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/perf/pnjl_criticality_feasibility_probe.jl [--output=PATH] [--repetitions=N]")
            return nothing
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    repetitions >= 1 || throw(ArgumentError("repetitions must be positive"))
    return (; output, repetitions)
end

function main(args=ARGS)
    options = parse_cli(args)
    options === nothing && return 0
    models_preloaded = isdefined(Main, :Models)

    include_start = time_ns()
    if !isdefined(Main, :PNJLCriticalityFeasibility)
        include(ANALYSIS_SCRIPT)
    end
    include_elapsed = time_ns() - include_start
    feasibility = Base.invokelatest(() -> getfield(Main, :PNJLCriticalityFeasibility))
    suite_runner = Base.invokelatest(() -> getfield(feasibility, :run_suite))
    case_factory = Base.invokelatest(() -> getfield(feasibility, :default_cases))
    case_evaluator = Base.invokelatest(() -> getfield(feasibility, :evaluate_case))
    cascade_evaluator = Base.invokelatest(() -> getfield(feasibility, :evaluate_cascade_case))

    first_start = time_ns()
    first_suite = Base.invokelatest(suite_runner; repetitions=1)
    first_elapsed = time_ns() - first_start
    first_case_count = length(first_suite.rows)
    first_fit_count = sum(row.fit_count for row in first_suite.rows)
    first_holdout_count = sum(row.fit_holdout_evaluation_count for row in first_suite.rows)
    first_targeted_count = sum(
        row.targeted_function_evaluations for row in first_suite.cascade_rows
    )
    first_cascade_fit_count = sum(row.polynomial_fits for row in first_suite.cascade_rows)

    cases = Base.invokelatest(case_factory)
    warm_start = time_ns()
    for _ in 1:options.repetitions, case in cases
        Base.invokelatest(case_evaluator, case)
    end
    warm_elapsed = time_ns() - warm_start
    case_count = length(cases)

    cascade_start = time_ns()
    for _ in 1:options.repetitions, case in cases
        Base.invokelatest(cascade_evaluator, case)
    end
    cascade_elapsed = time_ns() - cascade_start

    rows = [
        (
            phase="analysis_module_include",
            evidence_scope="local_process_startup",
            models_preloaded=models_preloaded,
            repetitions=1,
            case_evaluations=0,
            targeted_function_evaluations=0,
            polynomial_fits=0,
            fit_holdout_evaluations=0,
            elapsed_ns=include_elapsed,
            mean_time_ns=Float64(include_elapsed),
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            note="package/module loading only; no PNJL equilibrium solve",
        ),
        (
            phase="first_suite_call",
            evidence_scope="synthetic_first_call",
            models_preloaded=models_preloaded,
            repetitions=1,
            case_evaluations=missing,
            targeted_function_evaluations=missing,
            polynomial_fits=missing,
            fit_holdout_evaluations=missing,
            elapsed_ns=first_elapsed,
            mean_time_ns=missing,
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            note="aggregate suite orchestration; includes Models load and compilation when models_preloaded=false; no PNJL equilibrium solve",
        ),
        (
            phase="warm_suite",
            evidence_scope="synthetic_same_process",
            models_preloaded=true,
            repetitions=options.repetitions,
            case_evaluations=options.repetitions * case_count,
            targeted_function_evaluations=0,
            polynomial_fits=options.repetitions * first_fit_count,
            fit_holdout_evaluations=options.repetitions * first_holdout_count,
            elapsed_ns=warm_elapsed,
            mean_time_ns=warm_elapsed / (options.repetitions * case_count),
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            note="same-process full-window analytic classification; not a PNJL or PALC benchmark",
        ),
        (
            phase="warm_rho_support_cascade",
            evidence_scope="synthetic_same_process",
            models_preloaded=true,
            repetitions=options.repetitions,
            case_evaluations=options.repetitions * case_count,
            targeted_function_evaluations=options.repetitions * first_targeted_count,
            polynomial_fits=options.repetitions * first_cascade_fit_count,
            fit_holdout_evaluations=0,
            elapsed_ns=cascade_elapsed,
            mean_time_ns=cascade_elapsed / (options.repetitions * case_count),
            pnjl_residual_evaluations=missing,
            pnjl_jacobian_evaluations=missing,
            pnjl_newton_iterations=missing,
            pnjl_anchor_solves=missing,
            pnjl_branch_points=missing,
            note="same-process rho-support cascade; targeted counts are fixed-rho analytic evaluations, not PNJL solve costs",
        ),
    ]

    finite_checker = Base.invokelatest(() -> getfield(feasibility, :_assert_finite_rows))
    Base.invokelatest(finite_checker, rows, "process startup timings")
    mkpath(dirname(options.output))
    CSV.write(options.output, rows)
    evidence_dir = dirname(dirname(options.output))
    if basename(dirname(options.output)) == "tables" &&
       isfile(joinpath(evidence_dir, "manifest.json"))
        manifest_refresher = Base.invokelatest(() -> getfield(feasibility, :_refresh_manifest))
        Base.invokelatest(manifest_refresher, evidence_dir)
    end
    println("criticality feasibility performance probe: ", options.output)
    return 0
end

exit(main())
