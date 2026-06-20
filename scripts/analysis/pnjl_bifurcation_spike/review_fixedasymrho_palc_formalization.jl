const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !(REPO_ROOT in LOAD_PATH)
    insert!(LOAD_PATH, min(2, length(LOAD_PATH) + 1), REPO_ROOT)
end

include(joinpath(@__DIR__, "run_fixedasymrho_multibranch_palc.jl"))

module FixedAsymRhoPALCFormalizationReview

const M = Main.FixedAsymRhoMultiBranchPALCSpike
const U = Main.FixedMuPALCSpike
const DEFAULT_OUTPUT_REL = joinpath("data", "outputs", "results", "analysis", "palc_fixedasymrho_spike")

Base.@kwdef struct ReviewConfig
    repo_root::String
    output_dir::String
    run_id::String
    p_num::Int = 8
    t_num::Int = 4
    max_steps::Int = 80
    ds_rho::Float64 = 0.01
    dsmax_rho::Float64 = 0.04
    production_source_rho::Float64 = 1.0
    production_source_step::Float64 = 0.05
end

function _usage()
    return """
    Usage:
      julia --project=scripts/analysis/pnjl_bifurcation_spike scripts/analysis/pnjl_bifurcation_spike/review_fixedasymrho_palc_formalization.jl [options]

    Options:
      --p-num=<int>                       Momentum nodes. Default: 8.
      --t-num=<int>                       Theta nodes. Default: 4.
      --max-steps=<int>                   PALC max continuation steps. Default: 80.
      --ds-rho=<value>                    PALC initial pseudo-arclength scale. Default: 0.01.
      --dsmax-rho=<value>                 PALC max pseudo-arclength scale. Default: 0.04.
      --production-source-rho=<value>     Production-like anchor source rho/rho0. Default: 1.0.
      --production-source-step=<value>    Production-like anchor step. Default: 0.05.
      --output-dir=<path>                 Aggregate output directory.
    """
end

function _parse_args(args::Vector{String}; repo_root::String)
    raw = Dict{String, String}()
    for arg in args
        arg in ("-h", "--help") && (println(_usage()); exit(0))
        startswith(arg, "--") || throw(ArgumentError("unexpected positional argument: $(arg)"))
        parts = split(arg[3:end], "="; limit=2)
        length(parts) == 2 || throw(ArgumentError("expected --key=value, got $(arg)"))
        raw[replace(lowercase(parts[1]), "-" => "_")] = parts[2]
    end

    get_float(key, default) = haskey(raw, key) ? parse(Float64, raw[key]) : default
    get_int(key, default) = haskey(raw, key) ? parse(Int, raw[key]) : default

    run_id = U._run_id()
    output_dir = get(raw, "output_dir", nothing)
    output_dir = output_dir === nothing ?
        joinpath(repo_root, DEFAULT_OUTPUT_REL, "formalization_$(run_id)") :
        (isabspath(output_dir) ? output_dir : joinpath(repo_root, output_dir))

    return ReviewConfig(
        repo_root=repo_root,
        output_dir=output_dir,
        run_id=run_id,
        p_num=get_int("p_num", 8),
        t_num=get_int("t_num", 4),
        max_steps=get_int("max_steps", 80),
        ds_rho=get_float("ds_rho", 0.01),
        dsmax_rho=get_float("dsmax_rho", 0.04),
        production_source_rho=get_float("production_source_rho", 1.0),
        production_source_step=get_float("production_source_step", 0.05),
    )
end

function _scenarios()
    return (
        (
            name="t120_rho035",
            T_MeV=120.0,
            rho_anchor=0.35,
            rho_min=0.30,
            rho_max=0.40,
            rho_step=0.05,
        ),
        (
            name="t130_rho080",
            T_MeV=130.0,
            rho_anchor=0.80,
            rho_min=0.75,
            rho_max=0.85,
            rho_step=0.05,
        ),
    )
end

function _scenario_args(cfg::ReviewConfig, scenario)
    outdir = joinpath(cfg.output_dir, scenario.name)
    return String[
        "--T-MeV=$(scenario.T_MeV)",
        "--rho-anchor=$(scenario.rho_anchor)",
        "--rho-min=$(scenario.rho_min)",
        "--rho-max=$(scenario.rho_max)",
        "--rho-step=$(scenario.rho_step)",
        "--production-source-rho=$(cfg.production_source_rho)",
        "--production-source-step=$(cfg.production_source_step)",
        "--p-num=$(cfg.p_num)",
        "--t-num=$(cfg.t_num)",
        "--ds-rho=$(cfg.ds_rho)",
        "--dsmax-rho=$(cfg.dsmax_rho)",
        "--max-steps=$(cfg.max_steps)",
        "--run-phase3-review=true",
        "--output-dir=$(outdir)",
    ]
end

function _scenario_review(scenario, summary)
    phase3 = summary.phase3_review
    decision = phase3 === nothing ? nothing : phase3.decision
    return (
        name=scenario.name,
        T_MeV=scenario.T_MeV,
        rho_anchor=scenario.rho_anchor,
        output_dir=summary.config["output_dir"],
        distinct_root_count=length(summary.anchor_roots),
        groundstate_sample_count=summary.groundstate_sample_count,
        path_branch_count=summary.path_contract.branch_count,
        path_selection_count=summary.path_contract.selection_count,
        path_failure_count=summary.path_contract.failure_count,
        path_branch_jump_count=summary.path_contract.branch_jump_count,
        experimental_backend_candidate=summary.experimental_backend_candidate,
        phase3_decision=decision,
    )
end

function _aggregate_decision(scenario_reviews)
    all_have_multiroot = all(row -> row.distinct_root_count >= 2, scenario_reviews)
    all_have_selection = all(row -> row.path_selection_count > 0, scenario_reviews)
    any_failures = any(row -> row.path_failure_count > 0, scenario_reviews)
    ready_as_diagnostic = all_have_multiroot && all_have_selection && !any_failures
    return (
        recommendation=:keep_palc_as_isolated_analysis_backend,
        add_bifurcationkit_to_root_project=false,
        ready_for_production_replacement=false,
        ready_for_opt_in_diagnostic_backend=ready_as_diagnostic,
        reason=ready_as_diagnostic ?
            "Both target multiroot scenarios produced distinct PALC branches and explicit ground-state selections, but BifurcationKit should remain isolated until broader opt-in regression and precompile evidence justify a root dependency." :
            "The target scenarios did not yet satisfy all multibranch diagnostic gates, so PALC remains an analysis spike only.",
    )
end

function _write_aggregate_report(path::String, cfg::ReviewConfig, scenario_reviews, decision)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# FixedAsymmetricRho PALC Formalization Aggregate Review")
        println(io)
        println(io, "## Scenarios")
        println(io, "| scenario | T_MeV | rho_anchor | roots | selections | failures | jumps | diagnostic candidate |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---|")
        for row in scenario_reviews
            println(io, "| $(row.name) | $(row.T_MeV) | $(row.rho_anchor) | $(row.distinct_root_count) | $(row.path_selection_count) | $(row.path_failure_count) | $(row.path_branch_jump_count) | $(row.experimental_backend_candidate) |")
        end
        println(io)
        println(io, "## Decision")
        println(io, "- recommendation: $(decision.recommendation)")
        println(io, "- add_bifurcationkit_to_root_project: $(decision.add_bifurcationkit_to_root_project)")
        println(io, "- ready_for_production_replacement: $(decision.ready_for_production_replacement)")
        println(io, "- ready_for_opt_in_diagnostic_backend: $(decision.ready_for_opt_in_diagnostic_backend)")
        println(io, "- reason: $(decision.reason)")
        println(io)
        println(io, "## Output")
        println(io, "- aggregate_output_dir: $(cfg.output_dir)")
    end
    return path
end

function run(args::Vector{String}; repo_root::String)
    cfg = _parse_args(args; repo_root=repo_root)
    mkpath(cfg.output_dir)

    scenario_summaries = NamedTuple[]
    scenario_reviews = NamedTuple[]
    for scenario in _scenarios()
        summary = M.run(_scenario_args(cfg, scenario); repo_root=repo_root)
        push!(scenario_summaries, summary)
        push!(scenario_reviews, _scenario_review(scenario, summary))
    end

    decision = _aggregate_decision(scenario_reviews)
    summary_path = joinpath(cfg.output_dir, "fixedasymrho_palc_formalization_review.json")
    report_path = joinpath(cfg.output_dir, "fixedasymrho_palc_formalization_review.md")
    aggregate = (
        config=(
            output_dir=cfg.output_dir,
            p_num=cfg.p_num,
            t_num=cfg.t_num,
            max_steps=cfg.max_steps,
            ds_rho=cfg.ds_rho,
            dsmax_rho=cfg.dsmax_rho,
        ),
        scenario_reviews=scenario_reviews,
        decision=decision,
        artifacts=(
            report=report_path,
        ),
    )
    U._write_json(summary_path, aggregate)
    _write_aggregate_report(report_path, cfg, scenario_reviews, decision)

    println("FixedAsymmetricRho PALC formalization review output: $(cfg.output_dir)")
    println("recommendation: $(decision.recommendation)")
    println("add_bifurcationkit_to_root_project: $(decision.add_bifurcationkit_to_root_project)")
    return aggregate
end

function main_run(args::Vector{String}; repo_root::String)
    try
        run(args; repo_root=repo_root)
    catch err
        println(stderr, "FixedAsymmetricRho PALC formalization review failed.")
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        exit(1)
    end
    return nothing
end

end # module

using .FixedAsymRhoPALCFormalizationReview

if abspath(PROGRAM_FILE) == @__FILE__
    FixedAsymRhoPALCFormalizationReview.main_run(ARGS; repo_root=REPO_ROOT)
end
