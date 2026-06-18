# Lightweight CLI and output contract for run_combined_meson_density_scan.jl.
# This file intentionally avoids loading Models or numerical kernels.

module CombinedMesonDensityScanContract

export CombinedOptions,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_TRHO_ASYMMETRIC_OUTPUT_DIR,
    DEFAULT_REGIMES,
    OUTPUT_COLUMNS,
    parse_args,
    print_usage

const DEFAULT_OUTPUT_DIR = joinpath(
    "data", "outputs", "results", "relaxtime", "meson_density",
    "combined_tmu_mu0_temperature_scan",
)

const DEFAULT_TRHO_ASYMMETRIC_OUTPUT_DIR = joinpath(
    "data", "outputs", "results", "relaxtime", "meson_density",
    "combined_trho_asymmetric_smoke_scan",
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
    "constraint_mode", "rho_target", "rho_norm",
    "rho_u_fm3", "rho_d_fm3", "rho_s_fm3", "rho_u_over_rho_d",
    "asym_ud_ratio_target", "asym_s_target", "constraint_residual_norm",
    "equilibrium_pressure_fm4", "trho_seed_candidate_count", "trho_branch_policy",
    "mu_u_MeV", "mu_d_MeV", "mu_s_MeV", "muQ_MeV", "muS_MeV",
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
    rho_values::Vector{Float64}
    trho_reverse_rho::Bool
    asym_ud_ratio_target::Float64
    asym_s_target::Float64
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
    bose_x_min::Float64
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
    println("  --path <tmu|trho_asymmetric> Scan path strategy (default tmu)")
    println("  --regimes <list>            Comma-separated regimes")
    println("                              default stable,strict_bw_stage1,phase_shift_current,phase_shift_gbu_reference")
    println("  --tmin/--tmax/--tstep <MeV> Temperature range")
    println("  --muq <MeV>                 Quark chemical potential (default 0)")
    println("  --mub <MeV>                 Baryon chemical potential; sets muq=mub/3")
    println("  --muq-values <list>         Comma-separated fixed-muq paths")
    println("  --mumin/--mumax/--mustep    Fixed-muq path grid")
    println("  --rho-values <list>         Comma-separated rho/rho0 targets for trho_asymmetric")
    println("  --rhomin/--rhomax/--rhostep FixedAsymmetricRho target grid")
    println("  --trho-reverse-rho          Scan rho values in reverse input order for trho_asymmetric (default)")
    println("  --no-trho-reverse-rho       Keep rho values in input order for trho_asymmetric")
    println("  --asym-ud-ratio-target <x>  FixedAsymmetricRho rho_u/rho_d target (default 0.876)")
    println("  --asym-s-target <x>         FixedAsymmetricRho rho_s target in fm^-3 (default 0)")
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
    println("  --bose-x-min <float>        Bose x cutoff for phase-shift density_policy=x_min_cut (default 0)")
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
    step > 0.0 || throw(ArgumentError("range step must be positive"))
    hi >= lo || throw(ArgumentError("range maximum must be >= range minimum"))
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
    value in (:trho_asymmetric, :fixed_asymmetric_rho, :fixedasymrho, :asymmetric_trho) && return :trho_asymmetric
    throw(ArgumentError("unsupported path strategy $(value); currently supported: tmu, trho_asymmetric"))
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
        :output_dir_set => false,
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
        :rho_values => nothing,
        :rhomin => nothing,
        :rhomax => nothing,
        :rhostep => nothing,
        :trho_reverse_rho => true,
        :asym_ud_ratio_target => 0.876,
        :asym_s_target => 0.0,
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
        :bose_x_min => 0.0,
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
            opts[:output_dir_set] = true
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
        elseif arg == "--rho-values"
            opts[:rho_values] = _split_floats(require_value())
        elseif arg == "--rhomin"
            opts[:rhomin] = parse(Float64, require_value())
        elseif arg == "--rhomax"
            opts[:rhomax] = parse(Float64, require_value())
        elseif arg == "--rhostep"
            opts[:rhostep] = parse(Float64, require_value())
        elseif arg == "--trho-reverse-rho"
            opts[:trho_reverse_rho] = true
        elseif arg == "--no-trho-reverse-rho"
            opts[:trho_reverse_rho] = false
        elseif arg == "--asym-ud-ratio-target"
            opts[:asym_ud_ratio_target] = parse(Float64, require_value())
        elseif arg == "--asym-s-target"
            opts[:asym_s_target] = parse(Float64, require_value())
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
        elseif arg == "--bose-x-min"
            opts[:bose_x_min] = parse(Float64, require_value())
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

    path_strategy = _path_strategy_symbol(Symbol(opts[:path]))
    if path_strategy === :trho_asymmetric && !Bool(opts[:output_dir_set])
        opts[:output_dir] = DEFAULT_TRHO_ASYMMETRIC_OUTPUT_DIR
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

    any_rho_range_key = opts[:rhomin] !== nothing || opts[:rhomax] !== nothing || opts[:rhostep] !== nothing
    opts[:rho_values] !== nothing && any_rho_range_key && throw(ArgumentError("use either --rho-values or --rhomin/--rhomax/--rhostep, not both"))
    rho_values = if opts[:rho_values] !== nothing
        Float64.(opts[:rho_values])
    elseif any_rho_range_key
        opts[:rhomin] === nothing && throw(ArgumentError("--rhomin is required with rho range"))
        opts[:rhomax] === nothing && throw(ArgumentError("--rhomax is required with rho range"))
        opts[:rhostep] === nothing && throw(ArgumentError("--rhostep is required with rho range"))
        _range_values(Float64(opts[:rhomin]), Float64(opts[:rhomax]), Float64(opts[:rhostep]))
    else
        Float64[0.05]
    end
    all(isfinite, rho_values) || throw(ArgumentError("rho values must be finite"))
    all(>=(0.0), rho_values) || throw(ArgumentError("rho values must be nonnegative"))
    Float64(opts[:asym_ud_ratio_target]) > 0.0 || throw(ArgumentError("asym-ud-ratio-target must be positive"))
    isfinite(Float64(opts[:asym_s_target])) || throw(ArgumentError("asym-s-target must be finite"))

    if path_strategy === :tmu
        opts[:rho_values] === nothing || throw(ArgumentError("--rho-values is only valid with --path trho_asymmetric"))
        any_rho_range_key && throw(ArgumentError("--rhomin/--rhomax/--rhostep are only valid with --path trho_asymmetric"))
    else
        opts[:muq_values] === nothing || throw(ArgumentError("--muq-values is only valid with --path tmu"))
        any_range_key && throw(ArgumentError("--mumin/--mumax/--mustep are only valid with --path tmu"))
    end

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
    Float64(opts[:bose_x_min]) >= 0.0 || throw(ArgumentError("bose-x-min must be nonnegative"))
    Float64(opts[:gamma_zero_tol]) >= 0.0 || throw(ArgumentError("gamma-zero-tol must be nonnegative"))

    return CombinedOptions(
        String(opts[:output_dir]),
        opts[:figure_dir] === nothing ? _default_figure_dir(String(opts[:output_dir])) : String(opts[:figure_dir]),
        path_strategy,
        regimes,
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        tstep,
        Float64(opts[:muq]),
        muq_values,
        rho_values,
        Bool(opts[:trho_reverse_rho]),
        Float64(opts[:asym_ud_ratio_target]),
        Float64(opts[:asym_s_target]),
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
        Float64(opts[:bose_x_min]),
        Symbol(opts[:noanom_policy]),
        _strict_bw_stage_symbol(Symbol(opts[:strict_bw_stage])),
        Float64(opts[:gamma_zero_tol]),
        Bool(opts[:overwrite]),
        Bool(opts[:plot]),
    )
end

end # module CombinedMesonDensityScanContract
