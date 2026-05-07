#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(Main, :ScanOptions)
    include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_gap_transport_scan.jl"))
end

include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "phase_guided_transport_scan_cli.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "phase_guided_transport_scan_plan.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "phase_guided_transport_scan_assets.jl"))

const PhaseGuidedCLI = PhaseGuidedTransportScanCLI
const PhaseGuidedPlan = PhaseGuidedTransportScanPlan
const PhaseGuidedAssets = PhaseGuidedTransportScanAssets

function _default_output_paths(opts::PhaseGuidedCLI.PhaseGuidedScanOptions)
    result_csv = joinpath(opts.outdir, "phase_guided_transport_scan.csv")
    plan_csv = joinpath(opts.outdir, "sampling_plan.csv")
    readme = joinpath(opts.outdir, "README.md")
    effective_config = joinpath(opts.outdir, "effective_config.json")
    return (; result_csv, plan_csv, readme, effective_config)
end

function _build_runtime_scan_opts(result_csv::String, opts::PhaseGuidedCLI.PhaseGuidedScanOptions)
    base_args = String[
        "--output", result_csv,
        "--provenance-dir", opts.outdir,
        "--failed-points-output", joinpath(opts.outdir, "failed_points.csv"),
    ]
    opts.overwrite && push!(base_args, "--overwrite")
    opts.resume && push!(base_args, "--resume")
    return Main.parse_args(base_args)
end

function _point_meta(point)
    return (
        mode=point.mode,
        phase_reference_kind=point.phase_reference_kind,
        scan_group=point.scan_group,
        group_label=point.group_label,
        T_phase_base_MeV=point.T_phase_base_MeV,
        alpha_T=point.alpha_T,
    )
end

function run_phase_guided_scan(opts::PhaseGuidedCLI.PhaseGuidedScanOptions, ctx)
    paths = _default_output_paths(opts)
    mkpath(opts.outdir)
    mkpath(joinpath(opts.outdir, "figures"))

    if opts.overwrite
        for path in (paths.result_csv, paths.plan_csv, paths.readme, paths.effective_config)
            isfile(path) && rm(path)
        end
    end

    plan = PhaseGuidedPlan.build_plan(opts)
    PhaseGuidedAssets.write_plan_csv(paths.plan_csv, plan)
    PhaseGuidedAssets.write_readme(paths.readme, opts, plan; result_csv_name=basename(paths.result_csv))
    PhaseGuidedAssets.write_effective_config(paths.effective_config, PhaseGuidedAssets.build_effective_config(opts, paths.result_csv, paths.plan_csv))

    if opts.dry_run
        Main.ProvenanceMetadata.write_run_sidecars(
            opts.outdir;
            ctx=ctx,
            effective_config=PhaseGuidedAssets.build_effective_config(opts, paths.result_csv, paths.plan_csv),
            artifacts=[paths.plan_csv, paths.readme, paths.effective_config],
            summary=Dict{String,Any}("points_total" => plan.total, "dry_run" => true),
        )
        println("Phase-guided transport plan prepared (dry-run). Output: $(opts.outdir)")
        return plan
    end

    scan_opts = _build_runtime_scan_opts(paths.result_csv, opts)
    runtime = Main.build_scan_runtime(scan_opts)
    existing = scan_opts.resume && isfile(paths.result_csv) && !scan_opts.overwrite ?
        Main.read_existing_keys(paths.result_csv) :
        Set{Tuple{Float64,Float64,Float64}}()

    if scan_opts.overwrite && isfile(paths.result_csv)
        rm(paths.result_csv)
    end

    io = open(paths.result_csv, "a")
    failed_io = open(joinpath(opts.outdir, "failed_points.csv"), "a")
    try
        if filesize(paths.result_csv) == 0
            Main.ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "phase_guided_transport_scan",
                "script" => "scripts/relaxtime/run_phase_guided_transport_scan.jl",
                "git_commit" => Main.current_git_commit(),
            ))
            Main.write_header_if_needed(io)
        end
        if filesize(joinpath(opts.outdir, "failed_points.csv")) == 0
            Main.write_failed_points_header_if_needed(failed_io)
        end

        stats = PhaseGuidedPlan.execute_plan!(
            (point, previous_solution, previous_phase) ->
                Main.execute_gap_transport_scan_point!(
                    io,
                    nothing,
                    failed_io,
                    point.T_MeV,
                    point.muB_MeV,
                    point.xi,
                    scan_opts,
                    ctx,
                    runtime;
                    previous_solution=previous_solution,
                    previous_phase=previous_phase,
                    point_meta=_point_meta(point),
                ),
            opts,
            plan,
            existing,
        )

        Main.ProvenanceMetadata.write_run_sidecars(
            opts.outdir;
            ctx=ctx,
            effective_config=PhaseGuidedAssets.build_effective_config(opts, paths.result_csv, paths.plan_csv),
            artifacts=[paths.result_csv, paths.plan_csv, paths.readme, paths.effective_config, joinpath(opts.outdir, "failed_points.csv")],
            summary=Dict{String,Any}(
                "points_total" => stats.total,
                "success_count" => stats.success,
                "error_count" => stats.error,
                "skipped_count" => stats.skipped,
            ),
        )
    finally
        close(io)
        close(failed_io)
    end

    println("Phase-guided transport scan finished. Output: $(paths.result_csv)")
    return plan
end

function main()
    opts = PhaseGuidedCLI.parse_args(copy(ARGS))
    ctx = Main.ProvenanceMetadata.new_run_context("scripts/relaxtime/run_phase_guided_transport_scan.jl", copy(ARGS))
    run_phase_guided_scan(opts, ctx)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
