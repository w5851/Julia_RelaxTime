"""
Combined meson-density scan entrypoint.

This script keeps the scan-path axis separate from the density-regime axis:

- path strategy: currently `tmu`, a fixed-mu temperature scan in `(T, mu)` space;
- density regime: stable, strict BW, current phase shift, generalized BU.

It is intentionally a thin orchestrator over `Models` / `MesonDensityWorkflow`.
It does not reimplement gap solving, meson masses, widths, or density kernels.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using Dates
using .Constants_PNJL: ħc_MeV_fm
using .Models

const DEFAULT_OUTPUT_DIR = joinpath(
    "data", "outputs", "results", "relaxtime", "meson_density",
    "combined_tmu_mu0_temperature_scan",
)

const DEFAULT_REGIMES = [
    :stable,
    :strict_bw_stage1,
    :phase_shift_current,
    :phase_shift_gbu_reference,
]

const OUTPUT_COLUMNS = [
    "path_strategy", "path_point_index",
    "T_MeV", "muq_MeV", "muB_MeV", "xi",
    "flavor_profile", "meson_profile",
    "pi_channel", "k_channel", "charge_resolved",
    "mu_pi_MeV", "mu_K_MeV", "d_pi", "d_K",
    "regime", "phase_scheme", "strict_bw_stage",
    "real_axis_mode", "eta", "density_policy", "noanom_policy", "phase_convention", "phase_display",
    "qmax", "q_nodes", "omega_min", "omega_max", "omega_nodes",
    "gamma_zero_tol",
    "Phi", "Phibar", "m_u", "m_d", "m_s",
    "m_pi_MeV", "m_K_MeV", "gamma_pi_MeV", "gamma_K_MeV",
    "n_pi", "n_K", "kpi_ratio",
    "unsafe_bose_count", "min_E_minus_mu", "bose_x_min",
    "noanom_removed_component_count",
    "status", "message",
]

struct CombinedOptions
    output_dir::String
    figure_dir::String
    path_strategy::Symbol
    regimes::Vector{Symbol}
    tmin_MeV::Float64
    tmax_MeV::Float64
    tstep_MeV::Float64
    muq_MeV::Float64
    muq_values_MeV::Vector{Float64}
    xi::Float64
    flavor_profile::String
    meson_profile::String
    p_num::Int
    t_num::Int
    max_iter::Int
    stable_q_nodes::Int
    qmax::Float64
    q_nodes::Int
    omega_min::Float64
    omega_max::Float64
    omega_nodes::Int
    eta::Float64
    real_axis_mode::Symbol
    phase_convention::Symbol
    phase_display::Symbol
    density_policy::Symbol
    noanom_policy::Symbol
    strict_bw_stage::Symbol
    gamma_zero_tol::Float64
    overwrite::Bool
    plot::Bool
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_combined_meson_density_scan.jl [options]\n")
    println("Options:")
    println("  --output-dir <path>         Result output directory")
    println("  --figure-dir <path>         Figure output directory")
    println("                              default replaces data/outputs/results with data/outputs/figures")
    println("  --path <tmu>                Scan path strategy (default tmu)")
    println("  --regimes <list>            Comma-separated regimes")
    println("                              default stable,strict_bw_stage1,phase_shift_current,phase_shift_gbu_reference")
    println("  --tmin/--tmax/--tstep <MeV> Temperature range")
    println("  --muq <MeV>                 Quark chemical potential (default 0)")
    println("  --mub <MeV>                 Baryon chemical potential; sets muq=mub/3")
    println("  --muq-values <list>         Comma-separated fixed-muq paths")
    println("  --mumin/--mumax/--mustep    Fixed-muq path grid")
    println("  --xi <value>                Anisotropy, phase-shift density currently requires 0")
    println("  --flavor-profile <name>     config/physics/flavor_chemical profile (default default)")
    println("  --meson-profile <name>      config/physics/meson_chemical profile (default default)")
    println("  --p-num/--t-num <int>       Equilibrium integration nodes (default 8/4)")
    println("  --max-iter <int>            Solver/mass iterations (default 20)")
    println("  --stable-q-nodes <int>      Stable-particle q nodes (default 48)")
    println("  --qmax <fm^-1>              BW/phase-shift q cutoff (default 4)")
    println("  --q-nodes <int>             BW/phase-shift q nodes (default 6)")
    println("  --omega-min <fm^-1>         Phase-shift omega lower bound (default 0.05)")
    println("  --omega-max <fm^-1>         BW/phase-shift omega upper bound (default 4)")
    println("  --omega-nodes <int>         BW/phase-shift omega nodes (default 6)")
    println("  --real-axis-mode <mode>     finite_eta or pv_b0_eta0 (default pv_b0_eta0)")
    println("  --phase-convention <mode>   arg_propagator or arg_inverse_propagator (default arg_inverse_propagator)")
    println("  --phase-display <mode>      unwrapped or fold_0_pi (default unwrapped)")
    println("  --density-policy <policy>   strict_normal_domain/excitation_only_E_gt_mu/x_min_cut")
    println("  --noanom-policy <policy>    none or low_energy_branch_subtraction")
    println("  --eta <float>               finite_eta width parameter (ignored by pv_b0_eta0)")
    println("  --strict-bw-stage <stage>   stage1 or stage2 (default stage1)")
    println("  --gamma-zero-tol <float>    Strict BW stable fallback threshold (default 1e-12)")
    println("  --no-plot                   Do not emit SVG plot")
    println("  --overwrite                 Replace existing outputs")
    println("  -h, --help                  Show help")
end

function _split_symbols(raw::AbstractString)
    return Symbol[
        Symbol(strip(part))
        for part in split(raw, ','; keepempty=false)
        if !isempty(strip(part))
    ]
end

function _split_floats(raw::AbstractString)
    values = Float64[
        parse(Float64, strip(part))
        for part in split(raw, ','; keepempty=false)
        if !isempty(strip(part))
    ]
    isempty(values) && throw(ArgumentError("float list must contain at least one value"))
    return values
end

function _range_values(lo::Float64, hi::Float64, step::Float64)
    step > 0.0 || throw(ArgumentError("mustep must be positive"))
    hi >= lo || throw(ArgumentError("mumax must be >= mumin"))
    values = Float64[]
    x = lo
    while x <= hi + 1e-9
        push!(values, x)
        x += step
    end
    return values
end

function _default_figure_dir(output_dir::AbstractString)
    parts = collect(splitpath(normpath(String(output_dir))))
    if length(parts) >= 3
        for i in 1:(length(parts) - 2)
            if lowercase(parts[i]) == "data" && lowercase(parts[i + 1]) == "outputs" && lowercase(parts[i + 2]) == "results"
                parts[i + 2] = "figures"
                return joinpath(parts...)
            end
        end
    end
    return joinpath(String(output_dir), "figures")
end

function _path_strategy_symbol(value::Symbol)
    value in (:tmu, :Tmu, :temperature_mu, :temperature_scan) && return :tmu
    throw(ArgumentError("unsupported path strategy $(value); currently supported: tmu"))
end

function _regime_symbol(value::Symbol)
    value in (:stable, :stable_density) && return :stable
    value in (:strict_bw, :strict_bw_stage1, :strict_bw_reduced, :stage1_reduced) && return :strict_bw_stage1
    value in (:strict_bw_stage2, :strict_bw_qpole, :stage2_qpole) && return :strict_bw_stage2
    value in (:phase_shift_current, :current, :phase_e3) && return :phase_shift_current
    value in (:phase_shift_gbu, :phase_shift_gbu_reference, :gbu, :gbu_reference, :generalized_bu) && return :phase_shift_gbu_reference
    throw(ArgumentError("unsupported density regime: $(value)"))
end

function _strict_bw_stage_symbol(value::Symbol)
    value in (:stage1, :strict_bw_stage1, :stage1_reduced, :strict_bw_reduced) && return :stage1_reduced
    value in (:stage2, :strict_bw_stage2, :stage2_qpole, :strict_bw_qpole) && return :stage2_qpole
    throw(ArgumentError("unsupported strict BW stage: $(value)"))
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :output_dir => DEFAULT_OUTPUT_DIR,
        :figure_dir => nothing,
        :path => :tmu,
        :regimes => copy(DEFAULT_REGIMES),
        :tmin => 120.0,
        :tmax => 220.0,
        :tstep => 20.0,
        :muq => 0.0,
        :muq_values => nothing,
        :mumin => nothing,
        :mumax => nothing,
        :mustep => nothing,
        :xi => 0.0,
        :flavor_profile => "default",
        :meson_profile => "default",
        :p_num => 8,
        :t_num => 4,
        :max_iter => 20,
        :stable_q_nodes => 48,
        :qmax => 4.0,
        :q_nodes => 6,
        :omega_min => 0.05,
        :omega_max => 4.0,
        :omega_nodes => 6,
        :eta => 1e-6,
        :real_axis_mode => :pv_b0_eta0,
        :phase_convention => :arg_inverse_propagator,
        :phase_display => :unwrapped,
        :density_policy => :strict_normal_domain,
        :noanom_policy => :none,
        :strict_bw_stage => :stage1_reduced,
        :gamma_zero_tol => 1e-12,
        :overwrite => false,
        :plot => true,
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
        elseif arg == "--figure-dir"
            opts[:figure_dir] = require_value()
        elseif arg == "--path"
            opts[:path] = Symbol(require_value())
        elseif arg == "--regimes"
            opts[:regimes] = _split_symbols(require_value())
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--tstep"
            opts[:tstep] = parse(Float64, require_value())
        elseif arg == "--muq"
            opts[:muq] = parse(Float64, require_value())
        elseif arg == "--mub"
            opts[:muq] = parse(Float64, require_value()) / 3.0
        elseif arg == "--muq-values"
            opts[:muq_values] = _split_floats(require_value())
        elseif arg == "--mumin"
            opts[:mumin] = parse(Float64, require_value())
        elseif arg == "--mumax"
            opts[:mumax] = parse(Float64, require_value())
        elseif arg == "--mustep"
            opts[:mustep] = parse(Float64, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--flavor-profile"
            opts[:flavor_profile] = require_value()
        elseif arg == "--meson-profile"
            opts[:meson_profile] = require_value()
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--stable-q-nodes"
            opts[:stable_q_nodes] = parse(Int, require_value())
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
        elseif arg == "--eta"
            opts[:eta] = parse(Float64, require_value())
        elseif arg == "--real-axis-mode"
            opts[:real_axis_mode] = Symbol(require_value())
        elseif arg == "--phase-convention"
            opts[:phase_convention] = Symbol(require_value())
        elseif arg == "--phase-display"
            opts[:phase_display] = Symbol(require_value())
        elseif arg == "--density-policy"
            opts[:density_policy] = Symbol(require_value())
        elseif arg == "--noanom-policy"
            opts[:noanom_policy] = Symbol(require_value())
        elseif arg == "--strict-bw-stage"
            opts[:strict_bw_stage] = Symbol(require_value())
        elseif arg == "--gamma-zero-tol"
            opts[:gamma_zero_tol] = parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-plot"
            opts[:plot] = false
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    regimes = unique(_regime_symbol.(Vector{Symbol}(opts[:regimes])))
    isempty(regimes) && throw(ArgumentError("at least one regime must be selected"))
    any_range_key = opts[:mumin] !== nothing || opts[:mumax] !== nothing || opts[:mustep] !== nothing
    opts[:muq_values] !== nothing && any_range_key && throw(ArgumentError("use either --muq-values or --mumin/--mumax/--mustep, not both"))
    muq_values = if opts[:muq_values] !== nothing
        Float64.(opts[:muq_values])
    elseif any_range_key
        opts[:mumin] === nothing && throw(ArgumentError("--mumin is required with muq range"))
        opts[:mumax] === nothing && throw(ArgumentError("--mumax is required with muq range"))
        opts[:mustep] === nothing && throw(ArgumentError("--mustep is required with muq range"))
        _range_values(Float64(opts[:mumin]), Float64(opts[:mumax]), Float64(opts[:mustep]))
    else
        [Float64(opts[:muq])]
    end
    all(isfinite, muq_values) || throw(ArgumentError("muq values must be finite"))

    tstep = Float64(opts[:tstep])
    tstep > 0.0 || throw(ArgumentError("tstep must be positive"))
    Float64(opts[:tmax]) >= Float64(opts[:tmin]) || throw(ArgumentError("tmax must be >= tmin"))
    Int(opts[:p_num]) > 0 || throw(ArgumentError("p-num must be positive"))
    Int(opts[:t_num]) > 0 || throw(ArgumentError("t-num must be positive"))
    Int(opts[:max_iter]) > 0 || throw(ArgumentError("max-iter must be positive"))
    Int(opts[:stable_q_nodes]) > 1 || throw(ArgumentError("stable-q-nodes must be > 1"))
    Int(opts[:q_nodes]) > 1 || throw(ArgumentError("q-nodes must be > 1"))
    Int(opts[:omega_nodes]) > 1 || throw(ArgumentError("omega-nodes must be > 1"))
    Float64(opts[:qmax]) > 0.0 || throw(ArgumentError("qmax must be positive"))
    Float64(opts[:omega_max]) > Float64(opts[:omega_min]) || throw(ArgumentError("omega-max must exceed omega-min"))
    Float64(opts[:gamma_zero_tol]) >= 0.0 || throw(ArgumentError("gamma-zero-tol must be nonnegative"))

    return CombinedOptions(
        String(opts[:output_dir]),
        opts[:figure_dir] === nothing ? _default_figure_dir(String(opts[:output_dir])) : String(opts[:figure_dir]),
        _path_strategy_symbol(Symbol(opts[:path])),
        regimes,
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        tstep,
        Float64(opts[:muq]),
        muq_values,
        Float64(opts[:xi]),
        String(opts[:flavor_profile]),
        String(opts[:meson_profile]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:stable_q_nodes]),
        Float64(opts[:qmax]),
        Int(opts[:q_nodes]),
        Float64(opts[:omega_min]),
        Float64(opts[:omega_max]),
        Int(opts[:omega_nodes]),
        Float64(opts[:eta]),
        Symbol(opts[:real_axis_mode]),
        Symbol(opts[:phase_convention]),
        Symbol(opts[:phase_display]),
        Symbol(opts[:density_policy]),
        Symbol(opts[:noanom_policy]),
        _strict_bw_stage_symbol(Symbol(opts[:strict_bw_stage])),
        Float64(opts[:gamma_zero_tol]),
        Bool(opts[:overwrite]),
        Bool(opts[:plot]),
    )
end

@inline _fmt(x::Real) = Models.ScanCommon.fmt(x)
@inline _fmt(x::Integer) = string(x)
@inline _fmt(x::Symbol) = String(x)
@inline _fmt(x::Bool) = x ? "true" : "false"
@inline _fmt(::Nothing) = ""
@inline _fmt(x) = string(x)
@inline _rel(path::AbstractString) = replace(relpath(String(path), PROJECT_ROOT), '\\' => '/')

function _csv_escape(value)
    s = _fmt(value)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return string('"', replace(replace(s, '"' => "\"\""), '\n' => ' ', '\r' => ' '), '"')
    end
    return s
end

@inline function _getprop_or(x, field::Symbol, fallback)
    return hasproperty(x, field) ? getproperty(x, field) : fallback
end

@inline function _maybe_mev(value)
    value === nothing && return NaN
    v = Float64(value)
    return isfinite(v) ? v * ħc_MeV_fm : NaN
end

function _temperature_grid(opts::CombinedOptions)
    values = Float64[]
    T = opts.tmin_MeV
    while T <= opts.tmax_MeV + 1e-9
        push!(values, T)
        T += opts.tstep_MeV
    end
    return values
end

function _density_kwargs_for_profile(meson_profile, flavor_mev)
    chemical = Models.MesonChemicalProfiles.meson_chemical_profile_fm(meson_profile; flavor_mev=flavor_mev)
    return chemical, (
        pi_channel=chemical.pi_channel,
        k_channel=chemical.k_channel,
        μ_pi=chemical.mu_pi_fm,
        μ_K=chemical.mu_K_fm,
        d_pi=chemical.d_pi,
        d_K=chemical.d_K,
    )
end

function _solve_density_for_regime(regime::Symbol, meson_point, common_density, opts::CombinedOptions)
    if regime === :stable
        return Models.solve_meson_density_from_meson_point(
            meson_point;
            common_density...,
            num_q_nodes=opts.stable_q_nodes,
        )
    elseif regime === :strict_bw_stage1 || regime === :strict_bw_stage2
        stage = regime === :strict_bw_stage2 ? :stage2_qpole : opts.strict_bw_stage
        return Models.solve_strict_bw_meson_density_from_meson_point(
            meson_point;
            common_density...,
            stage=stage,
            qmax=opts.qmax,
            q_nodes=opts.q_nodes,
            omega_max=opts.omega_max,
            omega_nodes=opts.omega_nodes,
            gamma_zero_tol=opts.gamma_zero_tol,
        )
    elseif regime === :phase_shift_current
        return Models.solve_phase_shift_meson_density_from_meson_point(
            meson_point;
            common_density...,
            scheme=:current,
            qmax=opts.qmax,
            q_nodes=opts.q_nodes,
            omega_min=opts.omega_min,
            omega_max=opts.omega_max,
            omega_nodes=opts.omega_nodes,
            eta=opts.eta,
            real_axis_mode=opts.real_axis_mode,
            phase_convention=opts.phase_convention,
            phase_display=opts.phase_display,
            density_policy=opts.density_policy,
            noanom_policy=opts.noanom_policy,
        )
    end

    return Models.solve_phase_shift_meson_density_from_meson_point(
        meson_point;
        common_density...,
        scheme=:gbu_reference,
        qmax=opts.qmax,
        q_nodes=opts.q_nodes,
        omega_min=opts.omega_min,
        omega_max=opts.omega_max,
        omega_nodes=opts.omega_nodes,
        eta=opts.eta,
        real_axis_mode=opts.real_axis_mode,
        phase_convention=opts.phase_convention,
        phase_display=opts.phase_display,
        density_policy=opts.density_policy,
        noanom_policy=opts.noanom_policy,
    )
end

function _base_row(
    opts::CombinedOptions,
    point_index::Int,
    muq_MeV::Float64,
    T_MeV::Float64,
    flavor_profile,
    flavor_mev,
    meson_profile,
    chemical,
    meson_point,
)
    qp = meson_point.quark_params
    tp = meson_point.thermo_params
    return Dict{String, Any}(
        "path_strategy" => opts.path_strategy,
        "path_point_index" => point_index,
        "T_MeV" => T_MeV,
        "muq_MeV" => muq_MeV,
        "muB_MeV" => 3.0 * muq_MeV,
        "xi" => opts.xi,
        "flavor_profile" => flavor_profile.profile_name,
        "meson_profile" => meson_profile.profile_name,
        "pi_channel" => chemical.pi_channel,
        "k_channel" => chemical.k_channel,
        "charge_resolved" => chemical.charge_resolved,
        "mu_pi_MeV" => chemical.mu_pi_fm * ħc_MeV_fm,
        "mu_K_MeV" => chemical.mu_K_fm * ħc_MeV_fm,
        "d_pi" => chemical.d_pi,
        "d_K" => chemical.d_K,
        "Phi" => tp.Φ,
        "Phibar" => tp.Φbar,
        "m_u" => qp.m.u,
        "m_d" => qp.m.d,
        "m_s" => qp.m.s,
        "message" => "",
    )
end

function _density_row(base::Dict{String, Any}, regime::Symbol, density)
    pi_density = _getprop_or(density, :pi_density, nothing)
    k_density = _getprop_or(density, :k_density, nothing)
    return merge(copy(base), Dict{String, Any}(
        "regime" => regime,
        "phase_scheme" => _getprop_or(density, :scheme, ""),
        "strict_bw_stage" => _getprop_or(density, :stage, ""),
        "real_axis_mode" => _getprop_or(density, :real_axis_mode, ""),
        "eta" => _getprop_or(density, :eta, ""),
        "density_policy" => _getprop_or(density, :density_policy, :not_applicable),
        "noanom_policy" => _getprop_or(density, :noanom_policy, :none),
        "phase_convention" => _getprop_or(density, :phase_convention, ""),
        "phase_display" => _getprop_or(density, :phase_display, ""),
        "qmax" => _getprop_or(density, :qmax, ""),
        "q_nodes" => _getprop_or(density, :q_nodes, _getprop_or(density, :num_q_nodes, "")),
        "omega_min" => _getprop_or(density, :omega_min, ""),
        "omega_max" => _getprop_or(density, :omega_max, ""),
        "omega_nodes" => _getprop_or(density, :omega_nodes, ""),
        "gamma_zero_tol" => _getprop_or(density, :gamma_zero_tol, ""),
        "m_pi_MeV" => _maybe_mev(_getprop_or(density, :m_pi, nothing)),
        "m_K_MeV" => _maybe_mev(_getprop_or(density, :m_K, nothing)),
        "gamma_pi_MeV" => _maybe_mev(_getprop_or(density, :gamma_pi, nothing)),
        "gamma_K_MeV" => _maybe_mev(_getprop_or(density, :gamma_K, nothing)),
        "n_pi" => _getprop_or(density, :n_pi, NaN),
        "n_K" => _getprop_or(density, :n_K, NaN),
        "kpi_ratio" => _getprop_or(density, :kpi_ratio, NaN),
        "unsafe_bose_count" => _getprop_or(density, :unsafe_bose_count, ""),
        "min_E_minus_mu" => _getprop_or(density, :min_E_minus_mu, ""),
        "bose_x_min" => _getprop_or(density, :bose_x_min, ""),
        "noanom_removed_component_count" => _getprop_or(density, :noanom_removed_component_count, ""),
        "status" => _getprop_or(density, :status, :ok),
        "message" => _getprop_or(density, :message, ""),
        "pi_q_integral_estimate" => pi_density === nothing ? "" : _getprop_or(pi_density, :q_integral_estimate, ""),
        "K_q_integral_estimate" => k_density === nothing ? "" : _getprop_or(k_density, :q_integral_estimate, ""),
    ))
end

function _failure_row(base::Dict{String, Any}, regime::Symbol, err)
    msg = Models.ScanCommon.clean_message(sprint(showerror, err))
    return merge(copy(base), Dict{String, Any}(
        "regime" => regime,
        "phase_scheme" => "",
        "strict_bw_stage" => "",
        "real_axis_mode" => "",
        "eta" => "",
        "density_policy" => "",
        "noanom_policy" => "",
        "phase_convention" => "",
        "phase_display" => "",
        "qmax" => "",
        "q_nodes" => "",
        "omega_min" => "",
        "omega_max" => "",
        "omega_nodes" => "",
        "gamma_zero_tol" => "",
        "m_pi_MeV" => "NaN",
        "m_K_MeV" => "NaN",
        "gamma_pi_MeV" => "NaN",
        "gamma_K_MeV" => "NaN",
        "n_pi" => "NaN",
        "n_K" => "NaN",
        "kpi_ratio" => "NaN",
        "unsafe_bose_count" => "",
        "min_E_minus_mu" => "",
        "bose_x_min" => "",
        "noanom_removed_component_count" => "",
        "status" => :failed,
        "message" => msg,
    ))
end

function _row_values(row::Dict{String, Any})
    return [_csv_escape(get(row, col, "")) for col in OUTPUT_COLUMNS]
end

function _write_csv(path::String, opts::CombinedOptions, rows)
    open(path, "w") do io
        println(io, "# format: combined_meson_density_scan_v1")
        println(io, "# script: scripts/relaxtime/run_combined_meson_density_scan.jl")
        println(io, "# date: $(Dates.today())")
        println(io, "# bridge: path_strategy x density_regime")
        println(io, "# path_strategy: $(opts.path_strategy)")
        println(io, "# regimes: $(join(string.(opts.regimes), ','))")
        println(io, "# muq_values_MeV: $(join(_fmt.(opts.muq_values_MeV), ','))")
        println(io, "# real_axis_mode: $(opts.real_axis_mode)")
        println(io, "# phase_display: $(opts.phase_display)")
        println(io, "# density_policy: $(opts.density_policy)")
        println(io, "# noanom_policy: $(opts.noanom_policy)")
        println(io, "# gamma_zero_tol: $(opts.gamma_zero_tol)")
        println(io, join(OUTPUT_COLUMNS, ','))
        for row in rows
            println(io, join(_row_values(row), ','))
        end
    end
end

function _finite_rows(rows, field::String)
    out = Dict{String, Any}[]
    for row in rows
        string(get(row, "status", "")) == "ok" || continue
        y = tryparse(Float64, _fmt(get(row, field, "")))
        T = tryparse(Float64, _fmt(get(row, "T_MeV", "")))
        if y !== nothing && T !== nothing && isfinite(y) && isfinite(T) && y > 0.0
            push!(out, row)
        end
    end
    return out
end

@inline function _xml_escape(s)
    return replace(replace(replace(string(s), '&' => "&amp;"), '<' => "&lt;"), '>' => "&gt;")
end

function _regime_color(regime)
    r = Symbol(regime)
    r === :stable && return "#2563eb"
    r === :strict_bw_stage1 && return "#dc2626"
    r === :strict_bw_stage2 && return "#ea580c"
    r === :phase_shift_current && return "#059669"
    r === :phase_shift_gbu_reference && return "#7c3aed"
    return "#374151"
end

function _plot_panel(io, rows, field::String, title::String, x0::Float64, y0::Float64, w::Float64, h::Float64; logy::Bool)
    finite = _finite_rows(rows, field)
    println(io, "<g>")
    println(io, "<text x=\"$(x0)\" y=\"$(y0 - 14)\" font-size=\"14\" font-weight=\"600\">$(_xml_escape(title))</text>")
    println(io, "<rect x=\"$(x0)\" y=\"$(y0)\" width=\"$(w)\" height=\"$(h)\" fill=\"#ffffff\" stroke=\"#d1d5db\"/>")
    if isempty(finite)
        println(io, "<text x=\"$(x0 + 16)\" y=\"$(y0 + 32)\" font-size=\"12\" fill=\"#6b7280\">no finite ok rows</text>")
        println(io, "</g>")
        return
    end

    xs = [Float64(row["T_MeV"]) for row in finite]
    ys_raw = [Float64(row[field]) for row in finite]
    ys = logy ? log10.(ys_raw) : ys_raw
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    if xmin == xmax
        xmin -= 1.0
        xmax += 1.0
    end
    if ymin == ymax
        delta = max(abs(ymin) * 0.1, 1e-6)
        ymin -= delta
        ymax += delta
    end

    sx(x) = x0 + (Float64(x) - xmin) / (xmax - xmin) * w
    sy(y) = y0 + h - (Float64(y) - ymin) / (ymax - ymin) * h

    println(io, "<line x1=\"$(x0)\" y1=\"$(y0 + h)\" x2=\"$(x0 + w)\" y2=\"$(y0 + h)\" stroke=\"#9ca3af\"/>")
    println(io, "<line x1=\"$(x0)\" y1=\"$(y0)\" x2=\"$(x0)\" y2=\"$(y0 + h)\" stroke=\"#9ca3af\"/>")
    println(io, "<text x=\"$(x0)\" y=\"$(y0 + h + 18)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(xmin))</text>")
    println(io, "<text x=\"$(x0 + w - 42)\" y=\"$(y0 + h + 18)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(xmax))</text>")
    println(io, "<text x=\"$(x0 - 42)\" y=\"$(y0 + h)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(logy ? 10.0^ymin : ymin))</text>")
    println(io, "<text x=\"$(x0 - 42)\" y=\"$(y0 + 10)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(logy ? 10.0^ymax : ymax))</text>")

    regimes = unique(string(row["regime"]) for row in finite)
    for regime in regimes
        series = [row for row in finite if string(row["regime"]) == regime]
        sort!(series; by=row -> Float64(row["T_MeV"]))
        pts = String[]
        for row in series
            yv = Float64(row[field])
            push!(pts, "$(sx(row["T_MeV"])),$(sy(logy ? log10(yv) : yv))")
        end
        color = _regime_color(Symbol(regime))
        length(pts) > 1 && println(io, "<polyline points=\"$(join(pts, ' '))\" fill=\"none\" stroke=\"$(color)\" stroke-width=\"2\"/>")
        for pt in pts
            x, y = split(pt, ',')
            println(io, "<circle cx=\"$(x)\" cy=\"$(y)\" r=\"3\" fill=\"$(color)\"/>")
        end
    end
    println(io, "</g>")
end

@inline function _try_float_value(x)
    x isa Real && return Float64(x)
    return tryparse(Float64, _fmt(x))
end

function _heatmap_rows(rows, regime::String, field::String)
    out = Dict{String, Any}[]
    for row in rows
        string(get(row, "status", "")) == "ok" || continue
        string(get(row, "regime", "")) == regime || continue
        y = _try_float_value(get(row, field, ""))
        T = _try_float_value(get(row, "T_MeV", ""))
        μ = _try_float_value(get(row, "muq_MeV", ""))
        if y !== nothing && T !== nothing && μ !== nothing && isfinite(y) && isfinite(T) && isfinite(μ)
            push!(out, row)
        end
    end
    return out
end

function _grid_edges(xs::Vector{Float64})
    values = sort(unique(xs))
    isempty(values) && return Float64[]
    if length(values) == 1
        return [values[1] - 0.5, values[1] + 0.5]
    end
    edges = Vector{Float64}(undef, length(values) + 1)
    for i in 2:length(values)
        edges[i] = 0.5 * (values[i - 1] + values[i])
    end
    edges[1] = values[1] - (edges[2] - values[1])
    edges[end] = values[end] + (values[end] - edges[end - 1])
    return edges
end

@inline function _svg_rgb(r::Real, g::Real, b::Real)
    return "rgb($(round(Int, clamp(r, 0, 255))),$(round(Int, clamp(g, 0, 255))),$(round(Int, clamp(b, 0, 255))))"
end

function _heat_color(value::Float64, vmin::Float64, vmax::Float64)
    isfinite(value) || return "#9ca3af"
    value < vmin && return "#d1d5db"
    vmax <= vmin && return "#2563eb"
    t = clamp((value - vmin) / (vmax - vmin), 0.0, 1.0)
    stops = (
        (0.0, (37.0, 99.0, 235.0)),
        (0.35, (16.0, 185.0, 129.0)),
        (0.70, (250.0, 204.0, 21.0)),
        (1.0, (220.0, 38.0, 38.0)),
    )
    for i in 1:(length(stops) - 1)
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1
            u = (t - t0) / (t1 - t0)
            return _svg_rgb(
                c0[1] + u * (c1[1] - c0[1]),
                c0[2] + u * (c1[2] - c0[2]),
                c0[3] + u * (c1[3] - c0[3]),
            )
        end
    end
    return _svg_rgb(stops[end][2]...)
end

function _plot_heatmap_panel(io, rows, regime::String, x0::Float64, y0::Float64, w::Float64, h::Float64, vmin::Float64, vmax::Float64)
    panel_rows = _heatmap_rows(rows, regime, "kpi_ratio")
    println(io, "<g>")
    println(io, "<text x=\"$(x0)\" y=\"$(y0 - 14)\" font-size=\"14\" font-weight=\"600\">$(_xml_escape(regime))</text>")
    println(io, "<rect x=\"$(x0)\" y=\"$(y0)\" width=\"$(w)\" height=\"$(h)\" fill=\"#f3f4f6\" stroke=\"#d1d5db\"/>")
    if isempty(panel_rows)
        println(io, "<text x=\"$(x0 + 16)\" y=\"$(y0 + 32)\" font-size=\"12\" fill=\"#6b7280\">no finite ok rows</text>")
        println(io, "</g>")
        return
    end

    μs = sort(unique(Float64(_try_float_value(row["muq_MeV"])) for row in panel_rows))
    Ts = sort(unique(Float64(_try_float_value(row["T_MeV"])) for row in panel_rows))
    μ_edges = _grid_edges(μs)
    T_edges = _grid_edges(Ts)
    μmin, μmax = first(μ_edges), last(μ_edges)
    Tmin, Tmax = first(T_edges), last(T_edges)
    sx(x) = x0 + (Float64(x) - μmin) / (μmax - μmin) * w
    sy(y) = y0 + h - (Float64(y) - Tmin) / (Tmax - Tmin) * h

    for row in panel_rows
        μ = Float64(_try_float_value(row["muq_MeV"]))
        T = Float64(_try_float_value(row["T_MeV"]))
        value = Float64(_try_float_value(row["kpi_ratio"]))
        μ_idx = findfirst(==(μ), μs)
        T_idx = findfirst(==(T), Ts)
        μ_idx === nothing && continue
        T_idx === nothing && continue
        x_left = sx(μ_edges[μ_idx])
        x_right = sx(μ_edges[μ_idx + 1])
        y_top = sy(T_edges[T_idx + 1])
        y_bottom = sy(T_edges[T_idx])
        println(io, "<rect x=\"$(x_left)\" y=\"$(y_top)\" width=\"$(x_right - x_left)\" height=\"$(y_bottom - y_top)\" fill=\"$(_heat_color(value, vmin, vmax))\" stroke=\"#ffffff\" stroke-width=\"0.25\"/>")
    end

    println(io, "<rect x=\"$(x0)\" y=\"$(y0)\" width=\"$(w)\" height=\"$(h)\" fill=\"none\" stroke=\"#6b7280\"/>")
    println(io, "<text x=\"$(x0)\" y=\"$(y0 + h + 18)\" font-size=\"11\" fill=\"#4b5563\">mu_q: $(_fmt(first(μs)))..$(_fmt(last(μs))) MeV</text>")
    println(io, "<text x=\"$(x0 - 54)\" y=\"$(y0 + h)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(first(Ts)))</text>")
    println(io, "<text x=\"$(x0 - 54)\" y=\"$(y0 + 10)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(last(Ts)))</text>")
    println(io, "</g>")
end

function _write_svg_heatmap_plot(path::String, rows)
    regimes = sort(unique(string(row["regime"]) for row in rows))
    values = Float64[]
    for row in rows
        string(get(row, "status", "")) == "ok" || continue
        value = _try_float_value(get(row, "kpi_ratio", ""))
        value !== nothing && isfinite(value) && value >= 0.0 && push!(values, value)
    end
    vmax = isempty(values) ? 1.0 : maximum(values)
    vmax = max(vmax, 0.2)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1180\" height=\"860\" viewBox=\"0 0 1180 860\">")
        println(io, "<rect width=\"1180\" height=\"860\" fill=\"#f8fafc\"/>")
        println(io, "<text x=\"40\" y=\"38\" font-size=\"22\" font-weight=\"700\">Combined Meson Density Scan: FIG3-like T-mu heatmap</text>")
        println(io, "<text x=\"40\" y=\"62\" font-size=\"12\" fill=\"#475569\">Field: K/pi ratio; gray cells are missing, failed, negative, or below the color floor.</text>")
        panel_positions = ((72.0, 112.0), (650.0, 112.0), (72.0, 480.0), (650.0, 480.0))
        for (idx, regime) in enumerate(regimes[1:min(end, length(panel_positions))])
            x0, y0 = panel_positions[idx]
            _plot_heatmap_panel(io, rows, regime, x0, y0, 480.0, 280.0, 0.0, vmax)
        end
        xbar = 1040.0
        ybar = 112.0
        hbar = 280.0
        for i in 0:80
            y = ybar + hbar - i / 80 * hbar
            value = i / 80 * vmax
            println(io, "<rect x=\"$(xbar)\" y=\"$(y)\" width=\"24\" height=\"$(hbar / 80 + 0.5)\" fill=\"$(_heat_color(value, 0.0, vmax))\"/>")
        end
        println(io, "<rect x=\"$(xbar)\" y=\"$(ybar)\" width=\"24\" height=\"$(hbar)\" fill=\"none\" stroke=\"#6b7280\"/>")
        println(io, "<text x=\"$(xbar + 32)\" y=\"$(ybar + 10)\" font-size=\"11\" fill=\"#4b5563\">$(_fmt(vmax))</text>")
        println(io, "<text x=\"$(xbar + 32)\" y=\"$(ybar + hbar)\" font-size=\"11\" fill=\"#4b5563\">0</text>")
        println(io, "<text x=\"$(xbar)\" y=\"$(ybar - 14)\" font-size=\"12\" fill=\"#111827\">K/pi</text>")
        println(io, "</svg>")
    end
end

function _write_svg_plot(path::String, rows)
    finite_mu = sort(unique(Float64(_try_float_value(row["muq_MeV"])) for row in rows if _try_float_value(get(row, "muq_MeV", "")) !== nothing))
    if length(finite_mu) > 1
        return _write_svg_heatmap_plot(path, rows)
    end
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1180\" height=\"760\" viewBox=\"0 0 1180 760\">")
        println(io, "<rect width=\"1180\" height=\"760\" fill=\"#f8fafc\"/>")
        println(io, "<text x=\"40\" y=\"38\" font-size=\"22\" font-weight=\"700\">Combined Meson Density Scan: T sweep</text>")
        println(io, "<text x=\"40\" y=\"62\" font-size=\"12\" fill=\"#475569\">Rows with status=ok; density panels use logarithmic y scale.</text>")
        _plot_panel(io, rows, "n_pi", "n_pi vs T", 72.0, 110.0, 480.0, 230.0; logy=true)
        _plot_panel(io, rows, "n_K", "n_K vs T", 650.0, 110.0, 480.0, 230.0; logy=true)
        _plot_panel(io, rows, "kpi_ratio", "K/pi ratio vs T", 72.0, 450.0, 480.0, 230.0; logy=false)
        println(io, "<g>")
        println(io, "<text x=\"650\" y=\"454\" font-size=\"14\" font-weight=\"600\">Regime</text>")
        y = 484
        for regime in unique(string(row["regime"]) for row in rows)
            color = _regime_color(Symbol(regime))
            println(io, "<rect x=\"650\" y=\"$(y - 10)\" width=\"18\" height=\"4\" fill=\"$(color)\"/>")
            println(io, "<text x=\"678\" y=\"$(y)\" font-size=\"13\" fill=\"#111827\">$(_xml_escape(regime))</text>")
            y += 24
        end
        println(io, "</g>")
        println(io, "</svg>")
    end
end

function _json_escape(value)
    s = string(value)
    s = replace(s, '\\' => "\\\\")
    s = replace(s, '"' => "\\\"")
    s = replace(s, '\n' => "\\n")
    s = replace(s, '\r' => "\\r")
    return s
end

function _write_json_string(io, key::AbstractString, value; comma::Bool=true, indent::String="  ")
    suffix = comma ? "," : ""
    println(io, "$(indent)\"$(_json_escape(key))\": \"$(_json_escape(value))\"$(suffix)")
end

function _write_plot_manifest(path::String, opts::CombinedOptions, csv_path::String, figure_paths::Vector{String})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "{")
        _write_json_string(io, "format", "combined_meson_density_plot_manifest_v1")
        _write_json_string(io, "date", Dates.today())
        _write_json_string(io, "generated_by", "scripts/relaxtime/run_combined_meson_density_scan.jl")
        _write_json_string(io, "source_csv", _rel(csv_path))
        _write_json_string(io, "path_strategy", opts.path_strategy)
        println(io, "  \"density_regimes\": [")
        for (idx, regime) in enumerate(opts.regimes)
            suffix = idx == length(opts.regimes) ? "" : ","
            println(io, "    \"$(_json_escape(regime))\"$(suffix)")
        end
        println(io, "  ],")
        println(io, "  \"figures\": [")
        for (idx, figure_path) in enumerate(figure_paths)
            suffix = idx == length(figure_paths) ? "" : ","
            println(io, "    {")
            _write_json_string(io, "path", _rel(figure_path); indent="      ")
            _write_json_string(io, "kind", "quicklook_svg"; comma=false, indent="      ")
            println(io, "    }$(suffix)")
        end
        println(io, "  ]")
        println(io, "}")
    end
end

function _write_summary(path::String, opts::CombinedOptions, csv_path::String, plot_path::Union{Nothing, String}, plot_manifest_path::Union{Nothing, String}, rows)
    counts = Dict{String, Int}()
    for row in rows
        key = string(get(row, "status", ""))
        counts[key] = get(counts, key, 0) + 1
    end

    open(path, "w") do io
        println(io, "# Combined meson-density T-mu scan")
        println(io)
        println(io, "date: $(Dates.today())")
        println(io)
        println(io, "## Scope")
        println(io)
        println(io, "- path strategy: `$(opts.path_strategy)`")
        println(io, "- density regimes: `$(join(string.(opts.regimes), "`, `"))`")
        if length(opts.muq_values_MeV) == 1
            println(io, "- fixed chemical potential: `mu_q=$(only(opts.muq_values_MeV)) MeV`, `mu_B=$(3.0 * only(opts.muq_values_MeV)) MeV`")
        else
            println(io, "- fixed-mu path values: `mu_q=$(join(_fmt.(opts.muq_values_MeV), ",")) MeV`")
        end
        println(io, "- temperature range: `$(opts.tmin_MeV):$(opts.tstep_MeV):$(opts.tmax_MeV) MeV`")
        println(io, "- flavor profile: `$(opts.flavor_profile)`")
        println(io, "- meson profile: `$(opts.meson_profile)`")
        println(io)
        println(io, "## Outputs")
        println(io)
        println(io, "- result directory: `$(_rel(dirname(csv_path)))`")
        println(io, "- CSV: `$(_rel(csv_path))`")
        if plot_path !== nothing
            println(io, "- figure directory: `$(_rel(dirname(plot_path)))`")
            println(io, "- SVG: `$(_rel(plot_path))`")
        end
        if plot_manifest_path !== nothing
            println(io, "- plot manifest: `$(_rel(plot_manifest_path))`")
        end
        println(io)
        println(io, "## Regime Definitions")
        println(io)
        println(io, "- `stable`: stable particle limit from `Models.solve_meson_density_from_meson_point`.")
        println(io, "- `strict_bw_stage1`: reduced strict-BW density from `Models.solve_strict_bw_meson_density_from_meson_point`.")
        println(io, "- `phase_shift_current`: phase-shift BU density with weight `delta`.")
        println(io, "- `phase_shift_gbu_reference`: generalized BU density with weight `delta - 0.5sin(2delta)`.")
        println(io)
        println(io, "## Policies")
        println(io)
        println(io, "- `real_axis_mode=$(opts.real_axis_mode)`")
        println(io, "- `phase_convention=$(opts.phase_convention)`")
        println(io, "- `phase_display=$(opts.phase_display)`")
        println(io, "- `density_policy=$(opts.density_policy)`")
        println(io, "- `noanom_policy=$(opts.noanom_policy)`")
        println(io, "- `gamma_zero_tol=$(opts.gamma_zero_tol)`")
        println(io, "- This entrypoint is a Bridge-style composition of scan path and density-regime strategy.")
        println(io)
        println(io, "## Status Counts")
        println(io)
        for key in sort(collect(keys(counts)))
            println(io, "- `$(key)`: `$(counts[key])`")
        end
    end
end

function run_combined_scan(opts::CombinedOptions)
    opts.path_strategy === :tmu || throw(ArgumentError("only path_strategy=:tmu is currently implemented"))

    rows = Dict{String, Any}[]
    flavor_profile = Models.FlavorChemicalProfiles.load_flavor_chemical_profile(profile=opts.flavor_profile)
    meson_profile = Models.MesonChemicalProfiles.load_meson_chemical_profile(profile=opts.meson_profile)

    point_index = 0
    for muq_MeV in opts.muq_values_MeV
        flavor_mev = Models.FlavorChemicalProfiles.flavor_mu_profile_MeV(flavor_profile, muq_MeV)
        chemical, common_density = _density_kwargs_for_profile(meson_profile, flavor_mev)
        flavor_fm = Models.FlavorChemicalProfiles.flavor_mu_profile_fm(flavor_profile, muq_MeV)
        flavor_override = flavor_profile.apply_to_equilibrium ? (
            flavor_fm.mu_u_fm,
            flavor_fm.mu_d_fm,
            flavor_fm.mu_s_fm,
        ) : nothing

        continuation_state = nothing
        muq_fm = muq_MeV / ħc_MeV_fm
        for T_MeV in _temperature_grid(opts)
            point_index += 1
            T_fm = T_MeV / ħc_MeV_fm
            meson_point = try
                Models.solve_gap_and_meson_point(
                    T_fm,
                    muq_fm;
                    xi=opts.xi,
                    mesons=(chemical.pi_channel, chemical.k_channel),
                    continuation_state=continuation_state,
                    mixed_branch_align=:strict_sign_binding,
                    flavor_mu_override=flavor_override,
                    p_num=opts.p_num,
                    t_num=opts.t_num,
                    solver_kwargs=(; iterations=opts.max_iter),
                    mass_kwargs=(; iterations=opts.max_iter),
                )
            catch err
                base = Dict{String, Any}(
                    "path_strategy" => opts.path_strategy,
                    "path_point_index" => point_index,
                    "T_MeV" => T_MeV,
                    "muq_MeV" => muq_MeV,
                    "muB_MeV" => 3.0 * muq_MeV,
                    "xi" => opts.xi,
                    "flavor_profile" => flavor_profile.profile_name,
                    "meson_profile" => meson_profile.profile_name,
                    "pi_channel" => chemical.pi_channel,
                    "k_channel" => chemical.k_channel,
                    "charge_resolved" => chemical.charge_resolved,
                    "mu_pi_MeV" => chemical.mu_pi_fm * ħc_MeV_fm,
                    "mu_K_MeV" => chemical.mu_K_fm * ħc_MeV_fm,
                    "d_pi" => chemical.d_pi,
                    "d_K" => chemical.d_K,
                )
                for regime in opts.regimes
                    push!(rows, _failure_row(base, regime, err))
                end
                continue
            end

            continuation_state = meson_point.continuation_state
            base = _base_row(opts, point_index, muq_MeV, T_MeV, flavor_profile, flavor_mev, meson_profile, chemical, meson_point)
            for regime in opts.regimes
                try
                    density = _solve_density_for_regime(regime, meson_point, common_density, opts)
                    push!(rows, _density_row(base, regime, density))
                catch err
                    push!(rows, _failure_row(base, regime, err))
                end
            end
        end
    end

    return rows
end

function main()
    opts = parse_args(ARGS)
    outdir = normpath(joinpath(PROJECT_ROOT, opts.output_dir))
    figdir = normpath(joinpath(PROJECT_ROOT, opts.figure_dir))
    mkpath(outdir)
    csv_path = joinpath(outdir, "combined_meson_density_scan.csv")
    summary_path = joinpath(outdir, "README.md")
    plot_path = opts.plot ? joinpath(figdir, "combined_meson_density_scan.svg") : nothing
    plot_manifest_path = opts.plot ? joinpath(figdir, "plot_manifest.json") : nothing

    if !opts.overwrite && (isfile(csv_path) || isfile(summary_path) || (plot_path !== nothing && isfile(plot_path)) || (plot_manifest_path !== nothing && isfile(plot_manifest_path)))
        error("output exists; pass --overwrite to replace $(relpath(outdir, PROJECT_ROOT))")
    end

    rows = run_combined_scan(opts)
    _write_csv(csv_path, opts, rows)
    if plot_path !== nothing
        _write_svg_plot(plot_path, rows)
        _write_plot_manifest(plot_manifest_path, opts, csv_path, [plot_path])
    end
    _write_summary(summary_path, opts, csv_path, plot_path, plot_manifest_path, rows)

    println("csv=$(_rel(csv_path))")
    if plot_path !== nothing
        println("plot=$(_rel(plot_path))")
        println("plot_manifest=$(_rel(plot_manifest_path))")
    end
    println("summary=$(_rel(summary_path))")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
