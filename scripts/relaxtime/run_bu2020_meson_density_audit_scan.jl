"""
BU2020/temp7 mainline meson-density audit scan.

This script intentionally stays on the main `Models` / `MesonDensityWorkflow`
entrypoints. It emits a compact multi-prescription CSV plus README for auditing
real-axis mode, charged chemical potentials, and Bose-domain policy metadata.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using Dates
using .Constants_PNJL: ħc_MeV_fm
using .Models

const DEFAULT_OUTPUT_DIR = joinpath(
    "data", "outputs", "results", "relaxtime", "meson_density",
    "bu2020_temp7_mainline_audit",
)

const OUTPUT_COLUMNS = [
    "T_MeV", "muq_MeV", "muB_MeV", "xi",
    "mu_u_MeV", "mu_d_MeV", "mu_s_MeV",
    "pi_channel", "k_channel", "charge_resolved",
    "mu_pi_MeV", "mu_K_MeV", "chemical_profile", "flavor_profile",
    "branch_policy", "phase_side_label",
    "regime", "phase_scheme", "real_axis_mode", "eta",
    "density_policy", "noanom_policy", "phase_convention",
    "qmax", "q_nodes", "omega_min", "omega_max", "omega_nodes",
    "n_pi", "n_K", "kpi_ratio",
    "unsafe_bose_count", "min_E_minus_mu", "bose_x_min",
    "status", "message",
]

struct AuditOptions
    output_dir::String
    T_MeV::Float64
    muq_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    max_iter::Int
    qmax::Float64
    q_nodes::Int
    omega_min::Float64
    omega_max::Float64
    omega_nodes::Int
    overwrite::Bool
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_bu2020_meson_density_audit_scan.jl [options]\n")
    println("Options:")
    println("  --output-dir <path>         输出目录")
    println("  --T <MeV>                   温度 (default 90)")
    println("  --muq <MeV>                 light quark chemical potential (default 350)")
    println("  --xi <value>                anisotropy, currently phase-shift helper supports 0 (default 0)")
    println("  --p-num/--t-num <int>       equilibrium integration nodes (default 8/4)")
    println("  --max-iter <int>            solver/mass iterations (default 20)")
    println("  --qmax <fm^-1>              density q cutoff (default 4)")
    println("  --q-nodes <int>             density q nodes (default 4)")
    println("  --omega-min <fm^-1>         phase-shift omega lower bound (default 0.05)")
    println("  --omega-max <fm^-1>         phase-shift omega upper bound (default 3)")
    println("  --omega-nodes <int>         phase-shift omega nodes (default 4)")
    println("  --overwrite                 overwrite existing audit outputs")
    println("  -h, --help                  显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :output_dir => DEFAULT_OUTPUT_DIR,
        :T_MeV => 90.0,
        :muq_MeV => 350.0,
        :xi => 0.0,
        :p_num => 8,
        :t_num => 4,
        :max_iter => 20,
        :qmax => 4.0,
        :q_nodes => 4,
        :omega_min => 0.05,
        :omega_max => 3.0,
        :omega_nodes => 4,
        :overwrite => false,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            val = args[i + 1]
            i += 1
            return val
        end

        if arg == "--output-dir"
            opts[:output_dir] = require_value()
        elseif arg == "--T"
            opts[:T_MeV] = parse(Float64, require_value())
        elseif arg == "--muq"
            opts[:muq_MeV] = parse(Float64, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--qmax"
            opts[:qmax] = parse(Float64, require_value())
        elseif arg == "--q-nodes"
            opts[:q_nodes] = parse(Int, require_value())
        elseif arg == "--omega-min"
            opts[:omega_min] = parse(Float64, require_value())
        elseif arg == "--omega-max"
            opts[:omega_max] = parse(Float64, require_value())
        elseif arg == "--omega-nodes"
            opts[:omega_nodes] = parse(Int, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return AuditOptions(
        String(opts[:output_dir]),
        Float64(opts[:T_MeV]),
        Float64(opts[:muq_MeV]),
        Float64(opts[:xi]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Float64(opts[:qmax]),
        Int(opts[:q_nodes]),
        Float64(opts[:omega_min]),
        Float64(opts[:omega_max]),
        Int(opts[:omega_nodes]),
        Bool(opts[:overwrite]),
    )
end

fmt(x) = x === nothing ? "" : string(x)
fmt(x::Symbol) = String(x)
fmt(x::Bool) = x ? "true" : "false"

function row_values(row::Dict{String, Any})
    return [fmt(get(row, col, "")) for col in OUTPUT_COLUMNS]
end

function getprop_or(x, field::Symbol, fallback)
    return hasproperty(x, field) ? getproperty(x, field) : fallback
end

function density_row(base::Dict{String, Any}, regime::Symbol, density)
    return merge(copy(base), Dict{String, Any}(
        "regime" => regime,
        "phase_scheme" => getprop_or(density, :scheme, ""),
        "real_axis_mode" => getprop_or(density, :real_axis_mode, ""),
        "eta" => getprop_or(density, :eta, ""),
        "density_policy" => getprop_or(density, :density_policy, :not_applicable),
        "noanom_policy" => :none,
        "phase_convention" => getprop_or(density, :phase_convention, ""),
        "qmax" => getprop_or(density, :qmax, ""),
        "q_nodes" => getprop_or(density, :q_nodes, getprop_or(density, :num_q_nodes, "")),
        "omega_min" => getprop_or(density, :omega_min, ""),
        "omega_max" => getprop_or(density, :omega_max, ""),
        "omega_nodes" => getprop_or(density, :omega_nodes, ""),
        "n_pi" => getprop_or(density, :n_pi, NaN),
        "n_K" => getprop_or(density, :n_K, NaN),
        "kpi_ratio" => getprop_or(density, :kpi_ratio, NaN),
        "unsafe_bose_count" => getprop_or(density, :unsafe_bose_count, 0),
        "min_E_minus_mu" => getprop_or(density, :min_E_minus_mu, ""),
        "bose_x_min" => getprop_or(density, :bose_x_min, ""),
        "status" => getprop_or(density, :status, :ok),
        "message" => getprop_or(density, :message, ""),
    ))
end

function failure_row(base::Dict{String, Any}, regime::Symbol, err)
    return merge(copy(base), Dict{String, Any}(
        "regime" => regime,
        "phase_scheme" => "",
        "real_axis_mode" => "",
        "eta" => "",
        "density_policy" => "",
        "noanom_policy" => :none,
        "phase_convention" => "",
        "qmax" => "",
        "q_nodes" => "",
        "omega_min" => "",
        "omega_max" => "",
        "omega_nodes" => "",
        "n_pi" => "NaN",
        "n_K" => "NaN",
        "kpi_ratio" => "NaN",
        "unsafe_bose_count" => "",
        "min_E_minus_mu" => "",
        "bose_x_min" => "",
        "status" => :failed,
        "message" => replace(sprint(showerror, err), ',' => ';'),
    ))
end

function append_density!(rows, base, regime::Symbol, thunk)
    try
        push!(rows, density_row(base, regime, thunk()))
    catch err
        push!(rows, failure_row(base, regime, err))
    end
end

append_density!(thunk::Function, rows, base, regime::Symbol) = append_density!(rows, base, regime, thunk)

function write_summary(path::String, opts::AuditOptions, csv_path::String, rows)
    statuses = Dict{String, Int}()
    for row in rows
        key = fmt(row["status"])
        statuses[key] = get(statuses, key, 0) + 1
    end

    open(path, "w") do io
        println(io, "# BU2020/temp7 mainline meson-density audit")
        println(io)
        println(io, "date: $(Dates.today())")
        println(io)
        println(io, "## Scope")
        println(io)
        println(io, "- temp7 evidence: `D:\\Desktop\\Julia_RelaxTime\\temp7\\particles2020_bu_kpi_independent_audit`")
        println(io, "- output CSV: `$(relpath(csv_path, PROJECT_ROOT))`")
        println(io, "- state point: `T=$(opts.T_MeV) MeV`, `mu_q=$(opts.muq_MeV) MeV`, `xi=$(opts.xi)`")
        println(io, "- flavor profile: `bu2020_mu_s_0p2`")
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "- `real_axis_mode=:pv_b0_eta0` is emitted as a separate branch from finite-eta smoothing.")
        println(io, "- `density_policy=:strict_normal_domain` marks Bose-domain unsafe rows instead of silently skipping `omega <= mu_M`.")
        println(io, "- `density_policy=:excitation_only_E_gt_mu` is an explicit diagnostic continuation, not a literature fact.")
        println(io, "- no-anomalous subtraction is not enabled in this mainline audit output.")
        println(io)
        println(io, "## Status Counts")
        println(io)
        for key in sort(collect(keys(statuses)))
            println(io, "- `$(key)`: `$(statuses[key])`")
        end
    end
end

function main()
    opts = parse_args(ARGS)
    outdir = normpath(joinpath(PROJECT_ROOT, opts.output_dir))
    mkpath(outdir)
    csv_path = joinpath(outdir, "bu2020_temp7_meson_density_audit.csv")
    summary_path = joinpath(outdir, "README.md")

    if !opts.overwrite && (isfile(csv_path) || isfile(summary_path))
        error("output exists; pass --overwrite to replace $(relpath(outdir, PROJECT_ROOT))")
    end

    T_fm = opts.T_MeV / ħc_MeV_fm
    muq_fm = opts.muq_MeV / ħc_MeV_fm
    flavor_profile = Models.FlavorChemicalProfiles.load_flavor_chemical_profile(profile="bu2020_mu_s_0p2")
    flavor_mev = Models.FlavorChemicalProfiles.flavor_mu_profile_MeV(flavor_profile, opts.muq_MeV)
    flavor_fm = Models.FlavorChemicalProfiles.flavor_mu_profile_fm(flavor_profile, opts.muq_MeV)
    flavor_override = (
        flavor_fm.mu_u_fm,
        flavor_fm.mu_d_fm,
        flavor_fm.mu_s_fm,
    )

    profile_names = [
        "bu2020_kplus_over_piplus_mu_pi_100",
        "bu2020_kminus_over_piminus_mu_pi_100",
        "bu2020_kplus_over_piplus_mu_pi_134p5",
        "bu2020_kminus_over_piminus_mu_pi_134p5",
    ]

    rows = Dict{String, Any}[]
    for profile_name in profile_names
        chem_profile = Models.MesonChemicalProfiles.load_meson_chemical_profile(profile=profile_name)
        chem_fm = Models.MesonChemicalProfiles.meson_chemical_profile_fm(chem_profile; flavor_mev=flavor_mev)
        mu_K_MeV = chem_fm.mu_K_fm * ħc_MeV_fm
        base = Dict{String, Any}(
            "T_MeV" => opts.T_MeV,
            "muq_MeV" => opts.muq_MeV,
            "muB_MeV" => 3.0 * opts.muq_MeV,
            "xi" => opts.xi,
            "mu_u_MeV" => flavor_mev.mu_u_MeV,
            "mu_d_MeV" => flavor_mev.mu_d_MeV,
            "mu_s_MeV" => flavor_mev.mu_s_MeV,
            "pi_channel" => chem_fm.pi_channel,
            "k_channel" => chem_fm.k_channel,
            "charge_resolved" => chem_fm.charge_resolved,
            "mu_pi_MeV" => chem_profile.mu_pi_MeV,
            "mu_K_MeV" => mu_K_MeV,
            "chemical_profile" => chem_profile.profile_name,
            "flavor_profile" => flavor_profile.profile_name,
            "branch_policy" => "Models.solve_gap_and_meson_point default selector",
            "phase_side_label" => "not_forced",
        )

        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=opts.xi,
            mesons=(chem_fm.pi_channel, chem_fm.k_channel),
            mixed_branch_align=:strict_sign_binding,
            flavor_mu_override=flavor_override,
            p_num=opts.p_num,
            t_num=opts.t_num,
            solver_kwargs=(; iterations=opts.max_iter),
            mass_kwargs=(; iterations=opts.max_iter),
        )

        common_density = (
            pi_channel=chem_fm.pi_channel,
            k_channel=chem_fm.k_channel,
            μ_pi=chem_fm.mu_pi_fm,
            μ_K=chem_fm.mu_K_fm,
            d_pi=chem_fm.d_pi,
            d_K=chem_fm.d_K,
        )

        append_density!(rows, base, :stable) do
            Models.solve_meson_density_from_meson_point(meson_point; common_density..., num_q_nodes=opts.q_nodes)
        end
        append_density!(rows, base, :strict_bw_stage1_reduced) do
            Models.solve_strict_bw_meson_density_from_meson_point(
                meson_point;
                common_density...,
                stage=:stage1_reduced,
                qmax=opts.qmax,
                q_nodes=opts.q_nodes,
                omega_max=opts.omega_max,
                omega_nodes=opts.omega_nodes,
            )
        end
        append_density!(rows, base, :phase_shift_current_strict_pv) do
            Models.solve_phase_shift_meson_density_from_meson_point(
                meson_point;
                common_density...,
                scheme=:current,
                qmax=opts.qmax,
                q_nodes=opts.q_nodes,
                omega_min=opts.omega_min,
                omega_max=opts.omega_max,
                omega_nodes=opts.omega_nodes,
                real_axis_mode=:pv_b0_eta0,
                phase_convention=:arg_inverse_propagator,
                density_policy=:strict_normal_domain,
            )
        end
        append_density!(rows, base, :phase_shift_gbu_strict_pv) do
            Models.solve_phase_shift_meson_density_from_meson_point(
                meson_point;
                common_density...,
                scheme=:gbu_reference,
                qmax=opts.qmax,
                q_nodes=opts.q_nodes,
                omega_min=opts.omega_min,
                omega_max=opts.omega_max,
                omega_nodes=opts.omega_nodes,
                real_axis_mode=:pv_b0_eta0,
                phase_convention=:arg_inverse_propagator,
                density_policy=:strict_normal_domain,
            )
        end
        append_density!(rows, base, :phase_shift_gbu_excitation_only_pv) do
            Models.solve_phase_shift_meson_density_from_meson_point(
                meson_point;
                common_density...,
                scheme=:gbu_reference,
                qmax=opts.qmax,
                q_nodes=opts.q_nodes,
                omega_min=opts.omega_min,
                omega_max=opts.omega_max,
                omega_nodes=opts.omega_nodes,
                real_axis_mode=:pv_b0_eta0,
                phase_convention=:arg_inverse_propagator,
                density_policy=:excitation_only_E_gt_mu,
            )
        end
    end

    open(csv_path, "w") do io
        println(io, join(OUTPUT_COLUMNS, ','))
        for row in rows
            println(io, join(row_values(row), ','))
        end
    end
    write_summary(summary_path, opts, csv_path, rows)
    println("csv=$(relpath(csv_path, PROJECT_ROOT))")
    println("summary=$(relpath(summary_path, PROJECT_ROOT))")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
