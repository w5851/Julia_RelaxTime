"""
批量扫描 (T, μ) 网格，串联：各向异性 PNJL 平衡求解 → τ → RTA 输运系数。

输出：CSV（默认 data/outputs/results/relaxtime/scan/gap_transport_scan.csv）

单位约定：
- CLI 的温度/化学势默认使用 MeV（更符合扫描习惯）；脚本内部会换算到 fm⁻¹。
- 输出同时包含 MeV 与 fm⁻¹ 的关键量。

示例：
    julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 120 --tmax 400 --tstep 10 --mubmin 0 --mubmax 800 --mubstep 800 --xi-list -0.6,-0.4,-0.2,0.0,0.2,0.4,0.6 --mode finite_15 --overwrite --output data/outputs/results/relaxtime/scan/gap_transport_scan_xi-0p6to0p6.csv

注意：
- compute_bulk 默认开启；若需要更快的预览扫描，可显式传 `--no-compute-bulk`。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "scan_quality.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "provenance_metadata.jl"))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm, ρ0_inv_fm3
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, gauleg
using StaticArrays
using Dates
using .ScanCSV: ScanCSV
using .ScanQuality: assess_point_quality
using .ProvenanceMetadata: ProvenanceMetadata

const PNJL_MODEL = Models.create_model(:PNJL)
const TransportWorkflow = Models.transport_workflow_module()
const RT_ASR = Main.AverageScatteringRate
const RT_TCS = Main.TotalCrossSection
const REQUIRED_PROCESSES = TransportWorkflow.RelaxationTime.REQUIRED_PROCESSES

function preferred_phase_reference_path(dense_name::String, legacy_name::String)
    dense_path = joinpath(PROJECT_ROOT, "data", "reference", "pnjl", dense_name)
    return isfile(dense_path) ? dense_path : joinpath(PROJECT_ROOT, "data", "reference", "pnjl", legacy_name)
end

const DEFAULT_PHASE_BOUNDARY_PATH = preferred_phase_reference_path("boundary_dense.csv", "boundary.csv")
const DEFAULT_PHASE_CEP_PATH = preferred_phase_reference_path("cep_dense.csv", "cep.csv")
const DEFAULT_PHASE_CROSSOVER_PATH = preferred_phase_reference_path("crossover_dense.csv", "crossover.csv")

const MODULE_DEFAULT_P_NODES = RT_ASR.DEFAULT_P_NODES           # 20
const MODULE_DEFAULT_ANGLE_NODES = RT_ASR.DEFAULT_ANGLE_NODES   # 4
const MODULE_DEFAULT_PHI_NODES = RT_ASR.DEFAULT_PHI_NODES       # 8
const MODULE_DEFAULT_SIGMA_GRID_N = RT_ASR.DEFAULT_SIGMA_GRID_N # 60
const PHASE_BOUNDARY_XI_CACHE = Ref{Union{Nothing, Vector{Float64}}}(nothing)
const PHASE_CROSSOVER_XI_CACHE = Ref{Union{Nothing, Vector{Float64}}}(nothing)
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "gap_transport_scan_cli.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "gap_transport_scan_io.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "gap_transport_scan_orchestration.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "gap_transport_scan_plan.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "gap_transport_scan_phase_equilibrium.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "gap_transport_scan_provenance.jl"))

using .GapTransportScanCLI: ScanOptions, parse_args, print_usage
using .GapTransportScanIO: write_channel_diagnostics_header_if_needed,
    write_failed_points_header_if_needed,
    write_failed_point_row!,
    ensure_parent_dir,
    current_git_commit,
    ensure_output_header_compatible,
    write_header_if_needed
using .GapTransportScanOrchestration: build_scan_runtime, execute_gap_transport_scan_point!
using .GapTransportScanPlan: build_scan_plan, execute_scan_plan!
using .GapTransportScanPhaseEquilibrium: LocalPhaseTracker,
    tracker_phase,
    phase_structure,
    describe_seed_source,
    classify_phase_from_solution,
    solve_models_equilibrium,
    solve_equilibrium_with_diagnostics,
    is_phase_transition
using .GapTransportScanProvenance: build_effective_config, build_summary, collect_artifacts, write_scan_sidecars


@inline function _rate_with_alias(rates, key::Symbol)
    if hasproperty(rates, key)
        return getproperty(rates, key)
    end
    if key == :dubar_to_dubar && hasproperty(rates, :udbar_to_udbar)
        return getproperty(rates, :udbar_to_udbar)
    elseif key == :subar_to_subar && hasproperty(rates, :usbar_to_usbar)
        return getproperty(rates, :usbar_to_usbar)
    elseif key == :ubardbar_to_ubardbar && hasproperty(rates, :ud_to_ud)
        return getproperty(rates, :ud_to_ud)
    elseif key == :ubarubar_to_ubarubar && hasproperty(rates, :uu_to_uu)
        return getproperty(rates, :uu_to_uu)
    elseif key == :ubarsbar_to_ubarsbar && hasproperty(rates, :us_to_us)
        return getproperty(rates, :us_to_us)
    elseif key == :sbarsbar_to_sbarsbar && hasproperty(rates, :ss_to_ss)
        return getproperty(rates, :ss_to_ss)
    end
    return 0.0
end

function _fallback_relaxation_rate_contribution_rows(dens, rates)
    rows = NamedTuple[]

    function add_row(species::Symbol, channel::Symbol, density_key::Symbol, multiplicity::Float64)
        density = getproperty(dens, density_key)
        rate = Float64(_rate_with_alias(rates, channel))
        contribution = multiplicity * density * rate
        push!(rows, (
            species=species,
            channel=channel,
            density_key=density_key,
            multiplicity=multiplicity,
            density=density,
            rate=rate,
            contribution=contribution,
            total=0.0,
        ))
    end

    # u / d (isospin)
    add_row(:u, :uu_to_uu, :u, 1.0)
    add_row(:u, :ud_to_ud, :u, 1.0)
    add_row(:u, :uubar_to_uubar, :ubar, 1.0)
    add_row(:u, :uubar_to_ddbar, :ubar, 1.0)
    add_row(:u, :uubar_to_ssbar, :ubar, 1.0)
    add_row(:u, :udbar_to_udbar, :ubar, 1.0)
    add_row(:u, :us_to_us, :s, 1.0)
    add_row(:u, :usbar_to_usbar, :sbar, 1.0)

    add_row(:d, :uu_to_uu, :d, 1.0)
    add_row(:d, :ud_to_ud, :d, 1.0)
    add_row(:d, :uubar_to_uubar, :dbar, 1.0)
    add_row(:d, :uubar_to_ddbar, :dbar, 1.0)
    add_row(:d, :uubar_to_ssbar, :dbar, 1.0)
    add_row(:d, :udbar_to_udbar, :dbar, 1.0)
    add_row(:d, :us_to_us, :s, 1.0)
    add_row(:d, :usbar_to_usbar, :sbar, 1.0)

    # s
    add_row(:s, :us_to_us, :u, 2.0)
    add_row(:s, :usbar_to_usbar, :ubar, 2.0)
    add_row(:s, :ss_to_ss, :s, 1.0)
    add_row(:s, :ssbar_to_ssbar, :sbar, 1.0)
    add_row(:s, :ssbar_to_uubar, :sbar, 2.0)

    # anti-u / anti-d (isospin)
    add_row(:ubar, :uubar_to_uubar, :u, 1.0)
    add_row(:ubar, :uubar_to_ddbar, :u, 1.0)
    add_row(:ubar, :uubar_to_ssbar, :u, 1.0)
    add_row(:ubar, :dubar_to_dubar, :u, 1.0)
    add_row(:ubar, :ubardbar_to_ubardbar, :ubar, 1.0)
    add_row(:ubar, :ubarubar_to_ubarubar, :ubar, 1.0)
    add_row(:ubar, :subar_to_subar, :s, 1.0)
    add_row(:ubar, :ubarsbar_to_ubarsbar, :sbar, 1.0)

    add_row(:dbar, :uubar_to_uubar, :d, 1.0)
    add_row(:dbar, :uubar_to_ddbar, :d, 1.0)
    add_row(:dbar, :uubar_to_ssbar, :d, 1.0)
    add_row(:dbar, :dubar_to_dubar, :d, 1.0)
    add_row(:dbar, :ubardbar_to_ubardbar, :dbar, 1.0)
    add_row(:dbar, :ubarubar_to_ubarubar, :dbar, 1.0)
    add_row(:dbar, :subar_to_subar, :s, 1.0)
    add_row(:dbar, :ubarsbar_to_ubarsbar, :sbar, 1.0)

    # anti-s
    add_row(:sbar, :usbar_to_usbar, :u, 2.0)
    add_row(:sbar, :ubarsbar_to_ubarsbar, :ubar, 2.0)
    add_row(:sbar, :sbarsbar_to_sbarsbar, :sbar, 1.0)
    add_row(:sbar, :ssbar_to_ssbar, :s, 1.0)
    add_row(:sbar, :ssbar_to_uubar, :s, 2.0)

    totals = Dict{Symbol, Float64}()
    for row in rows
        totals[row.species] = get(totals, row.species, 0.0) + row.contribution
    end

    out = NamedTuple[]
    for row in rows
        push!(out, (
            species=row.species,
            channel=row.channel,
            density_key=row.density_key,
            multiplicity=row.multiplicity,
            density=row.density,
            rate=row.rate,
            contribution=row.contribution,
            total=get(totals, row.species, 0.0),
        ))
    end
    return out
end

function write_channel_diagnostics_rows!(io, T_mev::Float64, muq_mev::Float64, muB_mev::Float64, xi::Float64,
    dens, rates, tauinv, eq_backend, diag)
    rows = if isdefined(TransportWorkflow.RelaxationTime, :relaxation_rate_contribution_rows)
        TransportWorkflow.RelaxationTime.relaxation_rate_contribution_rows(dens, rates)
    else
        _fallback_relaxation_rate_contribution_rows(dens, rates)
    end
    for row in rows
        species = row.species
        tauinv_species = getproperty(tauinv, species)
        line = join([
            string(T_mev), string(muq_mev), string(muB_mev), string(xi),
            string(species), string(row.channel), string(row.density_key), string(row.multiplicity),
            string(row.density), string(row.rate), string(row.contribution), string(row.total), string(tauinv_species),
            string(eq_backend), string(diag.phase_curr), string(diag.phase_structure),
        ], ',')
        println(io, line)
    end
    flush(io)
end

function read_existing_keys(path::AbstractString)
    return ScanCSV.read_existing_keys(path, ["T_MeV", "muB_MeV", "xi"])
end

function csv_bool(x::Bool)
    return x ? "true" : "false"
end

function build_K_data(T_fm::Float64, mu_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    nodes = DEFAULT_MOMENTUM_NODES
    weights = DEFAULT_MOMENTUM_WEIGHTS
    A_u = A(masses.u, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    A_s = A(masses.s, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function integration_grids(opts::ScanOptions)
    if opts.integration_mode == :finite_15
        pg, pw = gauleg(0.0, 15.0, opts.tau_p_nodes)
        return (pg, pw, Λ_inv_fm)
    elseif opts.integration_mode == :finite_lambda
        pg, pw = gauleg(0.0, Λ_inv_fm, opts.tau_p_nodes)
        return (pg, pw, Λ_inv_fm)
    else
        return (nothing, nothing, nothing)
    end
end

function safe_total_cross_section(process::Symbol, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int, max_tries::Int=4)
    s_try = s
    last_err = nothing
    for _ in 1:max_tries
        try
            σ = RT_TCS.total_cross_section(process, s_try, quark_params, thermo_params, K_coeffs; n_points=n_points)
            isfinite(σ) && return σ
        catch err
            last_err = err
        end
        s_try = s_try * (1.0 + 1e-6) + 1e-10
    end
    @warn "failed to compute sigma; returning NaN" process=process s=s last_error=last_err
    return NaN
end

function build_sigma_caches(processes::Tuple, quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_sigma_points::Int, sigma_grid_n::Int)
    caches = Dict{Symbol,RT_ASR.CrossSectionCache}()
    for process in processes
        s_grid = RT_ASR.design_w0cdf_s_grid(process, quark_params, thermo_params;
            N=sigma_grid_n, p_cutoff=Λ_inv_fm)

        cache = RT_ASR.CrossSectionCache(process)
        n_ok, n_bad = 0, 0
        for s in s_grid
            σ = safe_total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
            if !isfinite(σ)
                n_bad += 1
                continue
            end
            RT_ASR.insert_sigma!(cache, s, σ)
            n_ok += 1
        end
        n_bad > 0 && @warn "sigma grid had non-finite points" process=process n_ok=n_ok n_bad=n_bad
        n_ok >= 2 || error("sigma cache has too few valid points for $process (n_ok=$n_ok)")
        caches[process] = cache
    end
    return caches
end

function run_scan(opts::ScanOptions, ctx::ProvenanceMetadata.RunContext)
    ensure_parent_dir(opts.output)
    if opts.channel_diagnostics_output !== nothing
        ensure_parent_dir(opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing
        ensure_parent_dir(opts.failed_points_output)
    end
    provenance_dir = isnothing(opts.provenance_dir) ? dirname(opts.output) : opts.provenance_dir
    mkpath(provenance_dir)

    if opts.resume && isfile(opts.output) && !opts.overwrite
        ensure_output_header_compatible(opts.output)
    end

    existing = opts.resume && isfile(opts.output) && !opts.overwrite ? read_existing_keys(opts.output) : Set{Tuple{Float64,Float64,Float64}}()

    if opts.overwrite && isfile(opts.output)
        rm(opts.output)
    end
    if opts.channel_diagnostics_output !== nothing && opts.overwrite && isfile(opts.channel_diagnostics_output)
        rm(opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing && opts.overwrite && isfile(opts.failed_points_output)
        rm(opts.failed_points_output)
    end

    new_file = !isfile(opts.output) || filesize(opts.output) == 0
    runtime = build_scan_runtime(opts)
    io = open(opts.output, "a")
    channel_io = opts.channel_diagnostics_output === nothing ? nothing : open(opts.channel_diagnostics_output, "a")
    failed_io = opts.failed_points_output === nothing ? nothing : open(opts.failed_points_output, "a")
    stats_success = 0
    stats_error = 0
    stats_skipped = 0
    try
        if new_file
            ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "gap_transport_scan",
                "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                "git_commit" => current_git_commit(),
                "provenance.entrypoint" => "workflow",
                "provenance.equilibrium_backend" => "models.solve_constraint(FixedMu) with phase-aware xi/T continuity",
                "provenance.tau_path" => "TransportWorkflow.solve_gap_and_transport",
                "provenance.integration_mode" => string(opts.integration_mode),
                "sigma_grid_n" => string(opts.sigma_grid_n),
                "integration_mode" => string(opts.integration_mode),
                "gc_every_n" => string(opts.gc_every_n),
                "tau_p_nodes" => string(opts.tau_p_nodes),
                "tau_angle_nodes" => string(opts.tau_angle_nodes),
                "tau_phi_nodes" => string(opts.tau_phi_nodes),
                "tau_n_sigma_points" => string(opts.tau_n_sigma_points),
                "tau_threshold_subtraction" => string(opts.tau_threshold_subtraction),
                "tau_asym_window" => string(opts.tau_asym_window),
                "tau_asym_fit_min_points" => string(opts.tau_asym_fit_min_points),
                "tau_asym_extra_points" => string(opts.tau_asym_extra_points),
                "tau_interpolation_mode" => string(opts.tau_interpolation_mode),
                "note.tau_threshold_hint" => "for near-threshold sharp channels, linear+threshold_subtraction often more robust than pchip",
                "tr_p_nodes" => string(opts.tr_p_nodes),
                "tr_p_max_fm" => string(opts.tr_p_max_fm),

                # labels for plotting convenience
                "y_label.sigma_over_T" => "σ/T",
                "y_label.sigma_over_T_over_eta_over_s" => "(σ/T)/(η/s)",
                "y_label.zeta_over_s_over_eta_over_s" => "(ζ/s)/(η/s)",
            ))
            write_header_if_needed(io)
        end
        if channel_io !== nothing
            channel_new_file = !isfile(opts.channel_diagnostics_output) || filesize(opts.channel_diagnostics_output) == 0
            if channel_new_file
                ScanCSV.write_metadata(channel_io, Dict(
                    "schema" => "scan_csv_v1",
                    "title" => "gap_transport_scan_channel_diagnostics",
                    "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                    "git_commit" => current_git_commit(),
                    "source_csv" => opts.output,
                ))
                write_channel_diagnostics_header_if_needed(channel_io)
            end
        end
        if failed_io !== nothing
            failed_new_file = !isfile(opts.failed_points_output) || filesize(opts.failed_points_output) == 0
            if failed_new_file
                ScanCSV.write_metadata(failed_io, Dict(
                    "schema" => "scan_csv_v1",
                    "title" => "gap_transport_scan_failed_points",
                    "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                    "git_commit" => current_git_commit(),
                    "source_csv" => opts.output,
                ))
                write_failed_points_header_if_needed(failed_io)
            end
        end

        plan = build_scan_plan(opts)
        stats = execute_scan_plan!(
            (T_mev, muB_mev, xi, previous_solution, previous_phase) ->
                execute_gap_transport_scan_point!(
                    io,
                    channel_io,
                    failed_io,
                    T_mev,
                    muB_mev,
                    xi,
                    opts,
                    ctx,
                    runtime;
                    previous_solution=previous_solution,
                    previous_phase=previous_phase,
                ),
            opts,
            plan,
            existing,
        )
        stats_success += stats.success
        stats_error += stats.error
        stats_skipped += stats.skipped
    finally
        close(io)
        if channel_io !== nothing
            close(channel_io)
        end
        if failed_io !== nothing
            close(failed_io)
        end
    end

    write_scan_sidecars(provenance_dir, ctx, opts, stats_success, stats_error, stats_skipped)

    println("Scan finished. Output: $(opts.output)")
end

function main()
    opts = parse_args(copy(ARGS))
    ctx = ProvenanceMetadata.new_run_context("scripts/relaxtime/run_gap_transport_scan.jl", copy(ARGS))
    run_scan(opts, ctx)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
