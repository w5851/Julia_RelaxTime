#!/usr/bin/env julia

using CSV
using JSON
using SHA

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

const CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
const DEFAULT_ANALYSIS_DIR = joinpath(
    PROJECT_ROOT,
    "docs",
    "analysis",
    "relaxtime",
    "phase_guided_transport",
    "phase_guided_transport_p128_xi001_analysis_v2",
)
const DEFAULT_CANDIDATE_CSV = joinpath(DEFAULT_ANALYSIS_DIR, "tables", "mechanism_window_candidates.csv")
const DEFAULT_OUT_DIR = joinpath(DEFAULT_ANALYSIS_DIR, "tables")

const BAND_EDGES = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, Inf]
const SAMPLE_DS_POINTS = [1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
const PRODUCTION_A_BUILDER_CONFIG = (p_nodes=16, p_max=20.0, cos_nodes=4, use_aniso=true)

Base.@kwdef struct MechanismOptions
    case_name::String = CASE
    candidate_csv::String = DEFAULT_CANDIDATE_CSV
    out_dir::String = DEFAULT_OUT_DIR
    window_limit::Int = 0
    channel_limit::Int = 4
    coverage_threshold::Float64 = 0.70
    ds_grid_n::Int = 60
    max_ds::Float64 = 2.0
    convergence_p_nodes::Int = 0
    convergence_angle_nodes::Int = 0
    convergence_phi_nodes::Int = 0
    convergence_n_sigma_points::Int = 0
    convergence_sigma_grid_n::Int = 0
    only_convergence_gate::Bool = false
    integration_mode::Symbol = :semi_infinite
end

function usage()
    println("""
    Usage:
      julia --project=. scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl [options]

    Options:
      --case-name <name>           formal result case name (default first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2)
      --candidate-csv <path>       mechanism_window_candidates.csv (default docs/analysis/.../tables)
      --out-dir <path>             output table directory (default docs/analysis/.../tables)
      --window-limit <int>         limit selected windows for smoke runs (0 means all selected windows)
      --channel-limit <int>        top channels per selected window (default 4)
      --coverage-threshold <float> cumulative contribution share target (default 0.70)
      --ds-grid-n <int>            denominator/sigma ds grid size (default 60)
      --max-ds <float>             maximum s-s_th sampled for denominator table (default 2.0)
      --convergence-p-nodes <int>  optional local high-precision rate gate p nodes
      --convergence-angle-nodes <int>
      --convergence-phi-nodes <int>
      --convergence-n-sigma-points <int>
      --convergence-sigma-grid-n <int>
      --only-convergence-gate      write only local_rate_convergence_gate.csv; leave denominator tables untouched
      --integration-mode <mode>    tau outer momentum integration: semi_infinite | finite_15 | finite_lambda (default semi_infinite)
    """)
end

function parse_args(argv=ARGS)
    opts = Dict{Symbol, Any}(
        :case_name => CASE,
        :candidate_csv => DEFAULT_CANDIDATE_CSV,
        :out_dir => DEFAULT_OUT_DIR,
        :window_limit => 0,
        :channel_limit => 4,
        :coverage_threshold => 0.70,
        :ds_grid_n => 60,
        :max_ds => 2.0,
        :convergence_p_nodes => 0,
        :convergence_angle_nodes => 0,
        :convergence_phi_nodes => 0,
        :convergence_n_sigma_points => 0,
        :convergence_sigma_grid_n => 0,
        :only_convergence_gate => false,
        :integration_mode => :semi_infinite,
    )
    i = 1
    while i <= length(argv)
        arg = argv[i]
        if arg in ("-h", "--help")
            usage()
            exit(0)
        end
        need() = (i += 1; i <= length(argv) || error("missing value for $arg"); argv[i])
        if arg == "--case-name"
            opts[:case_name] = need()
        elseif arg == "--candidate-csv"
            opts[:candidate_csv] = abspath(need())
        elseif arg == "--out-dir"
            opts[:out_dir] = abspath(need())
        elseif arg == "--window-limit"
            opts[:window_limit] = parse(Int, need())
        elseif arg == "--channel-limit"
            opts[:channel_limit] = parse(Int, need())
        elseif arg == "--coverage-threshold"
            opts[:coverage_threshold] = parse(Float64, need())
        elseif arg == "--ds-grid-n"
            opts[:ds_grid_n] = parse(Int, need())
        elseif arg == "--max-ds"
            opts[:max_ds] = parse(Float64, need())
        elseif arg == "--convergence-p-nodes"
            opts[:convergence_p_nodes] = parse(Int, need())
        elseif arg == "--convergence-angle-nodes"
            opts[:convergence_angle_nodes] = parse(Int, need())
        elseif arg == "--convergence-phi-nodes"
            opts[:convergence_phi_nodes] = parse(Int, need())
        elseif arg == "--convergence-n-sigma-points"
            opts[:convergence_n_sigma_points] = parse(Int, need())
        elseif arg == "--convergence-sigma-grid-n"
            opts[:convergence_sigma_grid_n] = parse(Int, need())
        elseif arg == "--only-convergence-gate"
            opts[:only_convergence_gate] = true
        elseif arg == "--integration-mode"
            mode = Symbol(need())
            mode in (:semi_infinite, :finite_15, :finite_lambda) || error("unknown integration mode: $(mode)")
            opts[:integration_mode] = mode
        else
            error("unknown argument: $arg")
        end
        i += 1
    end
    return MechanismOptions(;
        case_name=String(opts[:case_name]),
        candidate_csv=String(opts[:candidate_csv]),
        out_dir=String(opts[:out_dir]),
        window_limit=Int(opts[:window_limit]),
        channel_limit=Int(opts[:channel_limit]),
        coverage_threshold=Float64(opts[:coverage_threshold]),
        ds_grid_n=Int(opts[:ds_grid_n]),
        max_ds=Float64(opts[:max_ds]),
        convergence_p_nodes=Int(opts[:convergence_p_nodes]),
        convergence_angle_nodes=Int(opts[:convergence_angle_nodes]),
        convergence_phi_nodes=Int(opts[:convergence_phi_nodes]),
        convergence_n_sigma_points=Int(opts[:convergence_n_sigma_points]),
        convergence_sigma_grid_n=Int(opts[:convergence_sigma_grid_n]),
        only_convergence_gate=Bool(opts[:only_convergence_gate]),
        integration_mode=Symbol(opts[:integration_mode]),
    )
end

boolish(x) = x === true || (!(x === false || x === missing) && lowercase(string(x)) in ("true", "1", "yes"))
asstr(x) = x === missing ? "" : String(x)
asfloat(x) = x === missing || x == "" ? NaN : Float64(x)

function csv_cell(x)
    if x === missing
        return ""
    elseif x isa AbstractFloat
        return isfinite(x) ? string(x) : ""
    elseif x isa Bool
        return x ? "true" : "false"
    end
    text = string(x)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function write_csv_rows(path::AbstractString, fields::Vector{String}, rows::Vector{Dict{String, Any}})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(fields, ","))
        for row in rows
            println(io, join((csv_cell(get(row, field, "")) for field in fields), ","))
        end
    end
end

function trapz(x::Vector{Float64}, y::Vector{Float64})
    length(x) == length(y) || error("trapz length mismatch")
    length(x) <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(length(x) - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
end

function fval(row, name::Symbol)
    v = getproperty(row, name)
    return asfloat(v)
end

function sval(row, name::Symbol)
    return asstr(getproperty(row, name))
end

function closeval(a::Float64, b::Float64; atol::Float64=1.0e-6)
    return isfinite(a) && isfinite(b) && abs(a - b) <= atol
end

function selected_candidates(path::AbstractString, limit::Int)
    rows = [r for r in CSV.File(path) if boolish(getproperty(r, :selected_for_deep_scan))]
    sort!(rows, by=r -> begin
        order = getproperty(r, :selected_order)
        order === missing || order == "" ? 999 : (order isa Integer ? Int(order) : parse(Int, string(order)))
    end)
    return limit > 0 ? rows[1:min(limit, length(rows))] : rows
end

function mode_result_dirs(case_name::AbstractString)
    return Dict(
        "mode_a" => joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "transport", "phase_guided", "mode_a_fixed_muB_phase_scaled", case_name),
        "mode_b" => joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "transport", "phase_guided", "mode_b_fixed_T_sparse_muB", case_name),
    )
end

function mode_inputs(mode_key::String, opts::MechanismOptions)
    dir = mode_result_dirs(opts.case_name)[mode_key]
    return (
        main = collect(CSV.File(joinpath(dir, "phase_guided_transport_scan.csv"); comment="#")),
        channel = collect(CSV.File(joinpath(dir, "channel_diagnostics.csv"); comment="#")),
        config = JSON.parsefile(joinpath(dir, "effective_config.json")),
    )
end

function main_row_for_candidate(rows, cand)
    target_T = fval(cand, :T_MeV)
    target_muB = fval(cand, :muB_MeV)
    target_xi = fval(cand, :xi)
    matches = [r for r in rows if closeval(fval(r, :T_MeV), target_T; atol=1.0e-5) &&
        closeval(fval(r, :muB_MeV), target_muB; atol=1.0e-5) &&
        closeval(fval(r, :xi), target_xi; atol=1.0e-7)]
    isempty(matches) && error("main row not found for candidate $(sval(cand, :window_id))")
    return matches[1]
end

function neighbor_row(rows, target_row, xi::Float64)
    !isfinite(xi) && return nothing
    target_T = fval(target_row, :T_MeV)
    target_muB = fval(target_row, :muB_MeV)
    matches = [r for r in rows if closeval(fval(r, :T_MeV), target_T; atol=1.0e-5) &&
        closeval(fval(r, :muB_MeV), target_muB; atol=1.0e-5) &&
        closeval(fval(r, :xi), xi; atol=1.0e-7)]
    return isempty(matches) ? nothing : matches[1]
end

function channel_rows_for(channel_rows, main_row, species::String)
    target_T = fval(main_row, :T_MeV)
    target_muB = fval(main_row, :muB_MeV)
    target_xi = fval(main_row, :xi)
    return [r for r in channel_rows if sval(r, :species) == species &&
        closeval(fval(r, :T_MeV), target_T; atol=1.0e-5) &&
        closeval(fval(r, :muB_MeV), target_muB; atol=1.0e-5) &&
        closeval(fval(r, :xi), target_xi; atol=1.0e-7)]
end

function channel_map(rows)
    out = Dict{String, Any}()
    for r in rows
        out[sval(r, :channel)] = r
    end
    return out
end

function local_step(prev_map, target_map, next_map, channel::String)
    target = haskey(target_map, channel) ? fval(target_map[channel], :contribution) : NaN
    prev = haskey(prev_map, channel) ? fval(prev_map[channel], :contribution) : NaN
    nxt = haskey(next_map, channel) ? fval(next_map[channel], :contribution) : NaN
    vals = Float64[]
    isfinite(target) && isfinite(prev) && push!(vals, abs(target - prev))
    isfinite(target) && isfinite(nxt) && push!(vals, abs(nxt - target))
    return isempty(vals) ? 0.0 : maximum(vals)
end

function select_channels(channel_rows, main_rows, target_row, cand, opts::MechanismOptions)
    species = sval(cand, :primary_species)
    isempty(species) && return String[], 0.0, Dict{String, Any}()
    target_rows = channel_rows_for(channel_rows, target_row, species)
    prev_row = neighbor_row(main_rows, target_row, asfloat(getproperty(cand, :prev_xi)))
    next_row = neighbor_row(main_rows, target_row, asfloat(getproperty(cand, :next_xi)))
    prev_map = channel_map(prev_row === nothing ? [] : channel_rows_for(channel_rows, prev_row, species))
    target_map = channel_map(target_rows)
    next_map = channel_map(next_row === nothing ? [] : channel_rows_for(channel_rows, next_row, species))
    ranked = sort(collect(keys(target_map)); by=ch -> begin
        r = target_map[ch]
        abs(fval(r, :contribution)) + local_step(prev_map, target_map, next_map, ch)
    end, rev=true)

    selected = String[]
    coverage = 0.0
    for ch in ranked
        length(selected) >= opts.channel_limit && break
        push!(selected, ch)
        total = fval(target_map[ch], :total)
        contrib = fval(target_map[ch], :contribution)
        share = isfinite(total) && abs(total) > 0 ? contrib / total : 0.0
        coverage += max(share, 0.0)
        coverage >= opts.coverage_threshold && break
    end
    return selected, coverage, target_map
end

function p128_tau_config(config)
    opts = config["options"]
    return (
        p_nodes = Int(opts["tau_p_nodes"]),
        angle_nodes = Int(opts["tau_angle_nodes"]),
        phi_nodes = Int(opts["tau_phi_nodes"]),
        n_sigma_points = Int(opts["tau_n_sigma_points"]),
        sigma_grid_n = Int(opts["sigma_grid_n"]),
        propagator_xi_policy = Symbol(opts["propagator_xi_policy"]),
        sigma_cache_policy = Symbol(opts["sigma_cache_policy"]),
        threshold_subtraction = true,
        asym_window = 0.6,
        asym_fit_min_points = 8,
        asym_extra_points = 10,
        interpolation_mode = :linear,
        integration_mode = :semi_infinite,
    )
end

function convergence_enabled(opts::MechanismOptions)
    return opts.convergence_p_nodes > 0 ||
        opts.convergence_angle_nodes > 0 ||
        opts.convergence_phi_nodes > 0 ||
        opts.convergence_n_sigma_points > 0 ||
        opts.convergence_sigma_grid_n > 0
end

function convergence_tau_config(base, opts::MechanismOptions)
    return (
        p_nodes = opts.convergence_p_nodes > 0 ? opts.convergence_p_nodes : base.p_nodes,
        angle_nodes = opts.convergence_angle_nodes > 0 ? opts.convergence_angle_nodes : base.angle_nodes,
        phi_nodes = opts.convergence_phi_nodes > 0 ? opts.convergence_phi_nodes : base.phi_nodes,
        n_sigma_points = opts.convergence_n_sigma_points > 0 ? opts.convergence_n_sigma_points : base.n_sigma_points,
        sigma_grid_n = opts.convergence_sigma_grid_n > 0 ? opts.convergence_sigma_grid_n : base.sigma_grid_n,
        propagator_xi_policy = base.propagator_xi_policy,
        sigma_cache_policy = base.sigma_cache_policy,
        threshold_subtraction = base.threshold_subtraction,
        asym_window = base.asym_window,
        asym_fit_min_points = base.asym_fit_min_points,
        asym_extra_points = base.asym_extra_points,
        interpolation_mode = base.interpolation_mode,
        integration_mode = base.integration_mode,
    )
end

function with_integration_mode(base, integration_mode::Symbol)
    return (
        p_nodes = base.p_nodes,
        angle_nodes = base.angle_nodes,
        phi_nodes = base.phi_nodes,
        n_sigma_points = base.n_sigma_points,
        sigma_grid_n = base.sigma_grid_n,
        propagator_xi_policy = base.propagator_xi_policy,
        sigma_cache_policy = base.sigma_cache_policy,
        threshold_subtraction = base.threshold_subtraction,
        asym_window = base.asym_window,
        asym_fit_min_points = base.asym_fit_min_points,
        asym_extra_points = base.asym_extra_points,
        interpolation_mode = base.interpolation_mode,
        integration_mode = integration_mode,
    )
end

function tau_momentum_grid(tau_cfg)
    if tau_cfg.integration_mode === :finite_15
        return Main.GaussLegendre.gauleg(0.0, 15.0, tau_cfg.p_nodes)
    elseif tau_cfg.integration_mode === :finite_lambda
        return Main.GaussLegendre.gauleg(0.0, Main.Constants_PNJL.Λ_inv_fm, tau_cfg.p_nodes)
    elseif tau_cfg.integration_mode === :semi_infinite
        return nothing, nothing
    end
    error("unknown integration_mode=$(tau_cfg.integration_mode)")
end

@inline function meson_K(meson::Symbol, K_coeffs)
    if meson === :pi
        return K_coeffs.K123_plus
    elseif meson === :K
        return K_coeffs.K4567_plus
    elseif meson === :sigma_pi
        return K_coeffs.K123_minus
    elseif meson === :sigma_K
        return K_coeffs.K4567_minus
    end
    error("unsupported meson: $meson")
end

function build_state_from_main_row(row)
    T_fm = Float64(row.T_fm)
    muq_fm = Float64(row.muq_fm)
    xi = Float64(row.xi)
    masses = (u=Float64(row.m_u), d=Float64(row.m_d), s=Float64(row.m_s))
    Phi = Float64(row.Phi)
    Phibar = Float64(row.Phibar)
    qp_basic = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm))
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    A_vals = Main.AFieldBuilder.build_A_triplet(qp_basic, thermo_params;
        p_nodes=PRODUCTION_A_BUILDER_CONFIG.p_nodes,
        p_max=PRODUCTION_A_BUILDER_CONFIG.p_max,
        cos_nodes=PRODUCTION_A_BUILDER_CONFIG.cos_nodes,
        use_aniso=PRODUCTION_A_BUILDER_CONFIG.use_aniso,
    )

    A_u_iso = Main.A(masses.u, muq_fm, T_fm, Phi, Phibar, Main.DEFAULT_MOMENTUM_NODES_LIB, Main.DEFAULT_MOMENTUM_WEIGHTS_LIB)
    A_s_iso = Main.A(masses.s, muq_fm, T_fm, Phi, Phibar, Main.DEFAULT_MOMENTUM_NODES_LIB, Main.DEFAULT_MOMENTUM_WEIGHTS_LIB)
    G_u = Main.calculate_G_from_A(A_u_iso, masses.u)
    G_s = Main.calculate_G_from_A(A_s_iso, masses.s)
    K_coeffs = Main.calculate_effective_couplings(Main.G_fm2, Main.K_fm5, G_u, G_s)
    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=(u=A_vals.u, d=A_vals.d, s=A_vals.s))
    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function simple_p_branch_for_channel(process::Symbol, channel::Symbol, k0::Float64, k_norm::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    process_map = Main.Constants_PNJL.SCATTERING_MESON_MAP[process]
    channel_info = process_map[:channels][channel]
    simple_mesons = channel_info[:simple]
    T1, T2 = Main.TotalPropagator.get_flavor_factors_for_channel(process, channel)
    D_simple = 0.0 + 0.0im
    den_ref = NaN + NaN * im
    for meson in simple_mesons
        (meson === :pi || meson === :K) || continue
        pol = Main.TotalPropagator.get_polarization_params(meson, quark_params)
        Pi_re, Pi_im = Main.PolarizationCache.polarization_aniso_cached(
            pol.channel, k0, k_norm,
            pol.m1, pol.m2,
            pol.μ1, pol.μ2,
            thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
            pol.A1, pol.A2,
            pol.num_s_quark,
        )
        Pi = ComplexF64(Pi_re, Pi_im)
        den = 1.0 - 4.0 * meson_K(meson, K_coeffs) * Pi
        D_simple += T1 * Main.MesonPropagator.meson_propagator_simple(meson, K_coeffs, Pi) * T2
        isnan(real(den_ref)) && (den_ref = den)
    end
    return (D_simple=D_simple, den_ref=den_ref)
end

function mixed_p_branch_for_channel(process::Symbol, channel::Symbol, k0::Float64, k_norm::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    process_map = Main.Constants_PNJL.SCATTERING_MESON_MAP[process]
    channel_info = process_map[:channels][channel]
    if !Bool(channel_info[:mixed_P])
        return (has_mixed=false, D_mixed=0.0 + 0.0im, detM=NaN + NaN * im)
    end
    D_mixed = Main.TotalPropagator.total_propagator_mixed(
        process, channel, :P, k0, k_norm, quark_params, thermo_params, K_coeffs,
    )
    Pi_uu_re, Pi_uu_im = Main.PolarizationCache.polarization_aniso_cached(
        :P, k0, k_norm,
        quark_params.m.u, quark_params.m.u,
        quark_params.μ.u, quark_params.μ.u,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.u, quark_params.A.u,
        0,
    )
    Pi_ss_re, Pi_ss_im = Main.PolarizationCache.polarization_aniso_cached(
        :P, k0, k_norm,
        quark_params.m.s, quark_params.m.s,
        quark_params.μ.s, quark_params.μ.s,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.s, quark_params.A.s,
        2,
    )
    Pi_uu = ComplexF64(Pi_uu_re, Pi_uu_im)
    Pi_ss = ComplexF64(Pi_ss_re, Pi_ss_im)
    M00, M08, M88 = Main.MesonPropagator.calculate_coupling_elements(Pi_uu, Pi_ss, K_coeffs, :P)
    detM = M00 * M88 - M08 * M08
    return (has_mixed=true, D_mixed=D_mixed, detM=detM)
end

function ds_grid(opts::MechanismOptions)
    lower = collect(range(1.0e-4, min(0.5, opts.max_ds); length=max(20, opts.ds_grid_n ÷ 2)))
    upper = collect(range(min(0.5, opts.max_ds), opts.max_ds; length=max(20, opts.ds_grid_n ÷ 2)))
    return unique(sort(vcat(lower, upper, SAMPLE_DS_POINTS[SAMPLE_DS_POINTS .<= opts.max_ds])))
end

function band_index(ds::Float64)
    for i in 1:(length(BAND_EDGES) - 1)
        if ds > BAND_EDGES[i] && ds <= BAND_EDGES[i + 1]
            return i
        end
    end
    return length(BAND_EDGES) - 1
end

function process_branch_scan(window_id::String, process::Symbol, state, tau_cfg, opts::MechanismOptions)
    th = process_threshold_info(process, state.quark_params)
    grid = ds_grid(opts)
    sample_rows = Dict{String, Any}[]
    simple_inv = Float64[]
    mixed_inv = Float64[]
    d_simple = Float64[]
    d_mixed = Float64[]
    d_total = Float64[]
    sigma_vals = Float64[]
    for ds in grid
        s = th.s_th + ds
        tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
        t = 0.5 * (tb.t_min + tb.t_max)
        cms = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, state.quark_params)
        simple = simple_p_branch_for_channel(process, :s, cms.k0, cms.k, state.quark_params, state.thermo_params, state.K_coeffs)
        mixed = mixed_p_branch_for_channel(process, :s, cms.k0, cms.k, state.quark_params, state.thermo_params, state.K_coeffs)
        D_total = simple.D_simple + mixed.D_mixed
        sigma = Main.TotalCrossSection.total_cross_section(
            process, s, state.quark_params, state.thermo_params, state.K_coeffs;
            n_points=tau_cfg.n_sigma_points,
        )
        push!(simple_inv, isfinite(real(simple.den_ref)) ? 1.0 / max(abs2(simple.den_ref), 1.0e-30) : NaN)
        push!(mixed_inv, mixed.has_mixed ? 1.0 / max(abs2(mixed.detM), 1.0e-30) : NaN)
        push!(d_simple, abs2(simple.D_simple))
        push!(d_mixed, abs2(mixed.D_mixed))
        push!(d_total, abs2(D_total))
        push!(sigma_vals, sigma)
        if ds in SAMPLE_DS_POINTS || ds == grid[1] || ds == grid[end]
            push!(sample_rows, Dict{String, Any}(
                "window_id" => window_id,
                "channel" => String(process),
                "xi" => state.thermo_params.ξ,
                "ds" => ds,
                "s_value" => s,
                "den_simple_re" => real(simple.den_ref),
                "den_simple_im" => imag(simple.den_ref),
                "detM_re" => mixed.has_mixed ? real(mixed.detM) : NaN,
                "detM_im" => mixed.has_mixed ? imag(mixed.detM) : NaN,
                "invabs_den_simple" => simple_inv[end],
                "invabs_detM" => mixed_inv[end],
                "abs_D_simple_sq" => d_simple[end],
                "abs_D_mixed_sq" => d_mixed[end],
                "abs_D_total_sq" => d_total[end],
                "sigma" => sigma,
            ))
        end
    end
    sigma_idx = argmax(sigma_vals)
    simple_idx = argmax(replace(simple_inv, NaN => -Inf))
    mixed_idx = argmax(replace(mixed_inv, NaN => -Inf))
    simple_peak = simple_inv[simple_idx]
    mixed_peak = mixed_inv[mixed_idx]
    branch = (isfinite(mixed_peak) && mixed_peak > simple_peak) ? "mixed_detM" : "simple_1m4KPi"
    branch_idx = branch == "mixed_detM" ? mixed_idx : simple_idx
    aligned = band_index(grid[sigma_idx]) == band_index(grid[branch_idx])

    summary = Dict{String, Any}(
        "window_id" => window_id,
        "channel" => String(process),
        "xi" => state.thermo_params.ξ,
        "sigma_peak_ds" => grid[sigma_idx],
        "sigma_peak_value" => sigma_vals[sigma_idx],
        "simple_peak_ds" => grid[simple_idx],
        "simple_peak_invabs" => simple_peak,
        "mixed_peak_ds" => isfinite(mixed_peak) ? grid[mixed_idx] : NaN,
        "mixed_peak_invabs" => mixed_peak,
        "dominant_denominator_branch" => branch,
        "denominator_sigma_same_band" => aligned,
        "max_abs_D_simple_sq" => maximum(replace(d_simple, NaN => 0.0)),
        "max_abs_D_mixed_sq" => maximum(replace(d_mixed, NaN => 0.0)),
        "max_abs_D_total_sq" => maximum(replace(d_total, NaN => 0.0)),
    )

    band_rows = Dict{String, Any}[]
    for b in 1:(length(BAND_EDGES) - 1)
        idx = findall(ds -> ds > BAND_EDGES[b] && ds <= BAND_EDGES[b + 1], grid)
        length(idx) >= 2 || continue
        x = grid[idx]
        push!(band_rows, Dict{String, Any}(
            "window_id" => window_id,
            "channel" => String(process),
            "xi" => state.thermo_params.ξ,
            "band" => b,
            "ds_left" => BAND_EDGES[b],
            "ds_right" => BAND_EDGES[b + 1],
            "area_invabs_den_simple" => trapz(x, replace(simple_inv[idx], NaN => 0.0)),
            "area_invabs_detM" => trapz(x, replace(mixed_inv[idx], NaN => 0.0)),
            "area_abs_D_total_sq" => trapz(x, replace(d_total[idx], NaN => 0.0)),
            "area_sigma" => trapz(x, replace(sigma_vals[idx], NaN => 0.0)),
        ))
    end
    return summary, band_rows, sample_rows
end

function compute_rate_band_stats(process::Symbol, state, tau_cfg)
    p_grid, p_w = tau_momentum_grid(tau_cfg)
    cos_grid, cos_w = Main.GaussLegendre.gauleg(-1.0, 1.0, tau_cfg.angle_nodes)
    phi_grid, phi_w = Main.GaussLegendre.gauleg(0.0, 2 * pi, tau_cfg.phi_nodes)
    band_omega_ref = Ref(Vector{Float64}())
    band_omega_sigma_ref = Ref(Vector{Float64}())
    rate = Main.AverageScatteringRate.average_scattering_rate(
        process,
        state.quark_params,
        state.thermo_params,
        state.K_coeffs;
        p_nodes=tau_cfg.p_nodes,
        angle_nodes=tau_cfg.angle_nodes,
        phi_nodes=tau_cfg.phi_nodes,
        p_grid=p_grid,
        p_w=p_w,
        cos_grid=cos_grid,
        cos_w=cos_w,
        phi_grid=phi_grid,
        phi_w=phi_w,
        n_sigma_points=tau_cfg.n_sigma_points,
        sigma_grid_n=tau_cfg.sigma_grid_n,
        sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
        threshold_subtraction=tau_cfg.threshold_subtraction,
        asym_window=tau_cfg.asym_window,
        asym_fit_min_points=tau_cfg.asym_fit_min_points,
        asym_extra_points=tau_cfg.asym_extra_points,
        interpolation_mode=tau_cfg.interpolation_mode,
        band_edges=BAND_EDGES,
        band_omega_out=band_omega_ref,
        band_omega_sigma_out=band_omega_sigma_ref,
        propagator_xi_policy=tau_cfg.propagator_xi_policy,
        sigma_cache_policy=tau_cfg.sigma_cache_policy,
    )
    pi_sym, pj_sym, _, _ = Main.TotalCrossSection.parse_particles_from_process(process)
    mi = Main.AverageScatteringRate.get_mass(pi_sym, state.quark_params)
    mj = Main.AverageScatteringRate.get_mass(pj_sym, state.quark_params)
    mui = Main.AverageScatteringRate.get_mu(pi_sym, state.quark_params)
    muj = Main.AverageScatteringRate.get_mu(pj_sym, state.quark_params)
    rhoi = Main.AverageScatteringRate.number_density(pi_sym, mi, mui, state.thermo_params.T, state.thermo_params.Φ, state.thermo_params.Φbar, state.thermo_params.ξ; p_grid=p_grid, p_w=p_w, angle_nodes=tau_cfg.angle_nodes)
    rhoj = Main.AverageScatteringRate.number_density(pj_sym, mj, muj, state.thermo_params.T, state.thermo_params.Φ, state.thermo_params.Φbar, state.thermo_params.ξ; p_grid=p_grid, p_w=p_w, angle_nodes=tau_cfg.angle_nodes)
    prefactor = (Main.AverageScatteringRate.DQ^2) / (32.0 * pi^5 * rhoi * rhoj)
    return (rate=rate, prefactor=prefactor, omega=band_omega_ref[], omega_sigma=band_omega_sigma_ref[])
end

function direct_rate_for(channel_rows, main_row, species::String, channel::String)
    rows = channel_rows_for(channel_rows, main_row, species)
    matches = [r for r in rows if sval(r, :channel) == channel]
    return isempty(matches) ? NaN : fval(matches[1], :rate)
end

function append_convergence_rows!(conv_rows, window_id::String, mode_key::String, cand, selected::Vector{String},
    channel_rows, states::Dict{Float64, Any}, main_by_xi::Dict{Float64, Any}, high_cfg)
    species = sval(cand, :primary_species)
    target_xi = fval(cand, :xi)
    xi_values = sort(collect(keys(states)))
    for channel in selected
        process = Symbol(channel)
        prod_rate = Dict{Float64, Float64}()
        high_rate = Dict{Float64, Float64}()
        for xi in xi_values
            row = main_by_xi[xi]
            prod_rate[xi] = direct_rate_for(channel_rows, row, species, channel)
            stats = compute_rate_band_stats(process, states[xi], high_cfg)
            high_rate[xi] = stats.rate
        end
        neighbor_xis = [xi for xi in xi_values if !closeval(xi, target_xi; atol=1.0e-7)]
        prod_neighbors = [prod_rate[xi] for xi in neighbor_xis if isfinite(prod_rate[xi])]
        high_neighbors = [high_rate[xi] for xi in neighbor_xis if isfinite(high_rate[xi])]
        prod_target = get(prod_rate, target_xi, NaN)
        high_target = get(high_rate, target_xi, NaN)
        prod_neighbor_max = isempty(prod_neighbors) ? NaN : maximum(prod_neighbors)
        high_neighbor_max = isempty(high_neighbors) ? NaN : maximum(high_neighbors)
        prod_spike = isfinite(prod_target) && isfinite(prod_neighbor_max) && prod_neighbor_max > 0 ? prod_target / prod_neighbor_max : NaN
        high_spike = isfinite(high_target) && isfinite(high_neighbor_max) && high_neighbor_max > 0 ? high_target / high_neighbor_max : NaN
        rel_diff = isfinite(prod_target) && isfinite(high_target) && abs(prod_target) > 0 ? abs(high_target - prod_target) / abs(prod_target) : NaN
        ratio_retained = isfinite(prod_spike) && isfinite(high_spike) && prod_spike > 0 ? high_spike / prod_spike : NaN
        verdict = if isfinite(high_spike) && high_spike >= 3.0 && isfinite(ratio_retained) && ratio_retained >= 0.50
            "spike_persists_high_precision"
        elseif isfinite(high_spike) && high_spike >= 1.5
            "spike_weakened_but_visible"
        elseif isfinite(rel_diff) && rel_diff > 0.20
            "not_converged_or_spike_removed"
        else
            "inconclusive"
        end
        push!(conv_rows, Dict{String, Any}(
            "window_id" => window_id,
            "mode_key" => mode_key,
            "primary_species" => species,
            "channel" => channel,
            "target_xi" => target_xi,
            "neighbor_xis" => join(string.(neighbor_xis), ";"),
            "production_target_rate" => prod_target,
            "production_neighbor_max_rate" => prod_neighbor_max,
            "production_spike_ratio" => prod_spike,
            "high_target_rate" => high_target,
            "high_neighbor_max_rate" => high_neighbor_max,
            "high_spike_ratio" => high_spike,
            "target_rel_diff_high_vs_production" => rel_diff,
            "high_to_production_spike_ratio" => ratio_retained,
            "convergence_verdict" => verdict,
            "high_p_nodes" => high_cfg.p_nodes,
            "high_angle_nodes" => high_cfg.angle_nodes,
            "high_phi_nodes" => high_cfg.phi_nodes,
            "high_n_sigma_points" => high_cfg.n_sigma_points,
            "high_sigma_grid_n" => high_cfg.sigma_grid_n,
        ))
    end
end

function upstream_summary_row(window_id::String, rows)
    vals = Dict{String, Vector{Float64}}()
    names = ["T_MeV", "m_u", "m_s", "Phi", "Phibar", "K123_plus", "K4567_plus", "K0_plus", "K8_plus", "K08_plus"]
    for name in names
        vals[name] = Float64[]
    end
    for r in rows
        r === nothing && continue
        st = build_state_from_main_row(r)
        push!(vals["T_MeV"], fval(r, :T_MeV))
        push!(vals["m_u"], fval(r, :m_u))
        push!(vals["m_s"], fval(r, :m_s))
        push!(vals["Phi"], fval(r, :Phi))
        push!(vals["Phibar"], fval(r, :Phibar))
        push!(vals["K123_plus"], st.K_coeffs.K123_plus)
        push!(vals["K4567_plus"], st.K_coeffs.K4567_plus)
        push!(vals["K0_plus"], st.K_coeffs.K0_plus)
        push!(vals["K8_plus"], st.K_coeffs.K8_plus)
        push!(vals["K08_plus"], st.K_coeffs.K08_plus)
    end
    max_step = 0.0
    max_curv = 0.0
    driver = ""
    for name in names
        v = vals[name]
        length(v) >= 2 || continue
        denom = max(maximum(abs.(v)), 1.0e-12)
        step = maximum(abs.(diff(v))) / denom
        if step > max_step
            max_step = step
            driver = name
        end
        if length(v) >= 3
            curv = abs(v[3] - 2.0 * v[2] + v[1]) / denom
            max_curv = max(max_curv, curv)
        end
    end
    return Dict{String, Any}(
        "window_id" => window_id,
        "n_points" => length([r for r in rows if r !== nothing]),
        "max_rel_step" => max_step,
        "max_rel_curvature" => max_curv,
        "dominant_background_driver" => driver,
        "upstream_branch_flag" => max_step > 0.15 || max_curv > 0.15,
    )
end

function classify_window(coverage::Float64, rate_err::Float64, branch_aligned::Bool, upstream_flag::Bool)
    rate_ok = isfinite(rate_err) && rate_err <= 0.05
    if coverage >= 0.50 && rate_ok && branch_aligned
        return "small_denominator_supported", 0.95
    elseif upstream_flag && branch_aligned
        return "upstream_branch_candidate", 0.70
    elseif coverage >= 0.50 && rate_ok
        return "channel_competition_supported", 0.55
    elseif upstream_flag
        return "phase_space_or_density_supported", 0.45
    end
    return "inconclusive", 0.20
end

function main()
    opts = parse_args()
    opts.only_convergence_gate && !convergence_enabled(opts) && error("--only-convergence-gate requires convergence override options")
    mkpath(opts.out_dir)
    cands = selected_candidates(opts.candidate_csv, opts.window_limit)
    isempty(cands) && error("no selected candidates found in $(opts.candidate_csv)")

    inputs = Dict{String, Any}()
    for mode_key in keys(mode_result_dirs(opts.case_name))
        inputs[mode_key] = mode_inputs(mode_key, opts)
    end

    mechanism_rows = Dict{String, Any}[]
    denom_rows = Dict{String, Any}[]
    band_rows = Dict{String, Any}[]
    sample_rows = Dict{String, Any}[]
    upstream_rows = Dict{String, Any}[]
    global_rows = Dict{String, Any}[]
    convergence_rows = Dict{String, Any}[]

    for cand in cands
        window_id = sval(cand, :window_id)
        mode_key = sval(cand, :mode_key)
        data = inputs[mode_key]
        main_rows = data.main
        channel_rows = data.channel
        tau_cfg = with_integration_mode(p128_tau_config(data.config), opts.integration_mode)
        target_row = main_row_for_candidate(main_rows, cand)
        prev_row = neighbor_row(main_rows, target_row, asfloat(getproperty(cand, :prev_xi)))
        next_row = neighbor_row(main_rows, target_row, asfloat(getproperty(cand, :next_xi)))
        upstream = upstream_summary_row(window_id, [prev_row, target_row, next_row])
        push!(upstream_rows, upstream)

        selected, coverage, target_map = select_channels(channel_rows, main_rows, target_row, cand, opts)
        if isempty(selected)
            verdict, score = "inconclusive", 0.0
            push!(mechanism_rows, Dict{String, Any}(
                "window_id" => window_id,
                "window_label" => sval(cand, :window_label),
                "mode_key" => mode_key,
                "plot_panel" => sval(cand, :plot_panel),
                "plot_series" => sval(cand, :plot_series),
                "observable" => sval(cand, :observable),
                "primary_species" => sval(cand, :primary_species),
                "xi" => fval(cand, :xi),
                "mechanism_verdict" => verdict,
                "evidence_score" => score,
                "dominant_channels" => "",
                "channel_coverage_share" => coverage,
                "dominant_denominator_branch" => "",
                "max_rate_reproduction_rel_error" => "",
                "denominator_sigma_alignment" => false,
                "upstream_branch_flag" => upstream["upstream_branch_flag"],
                "mechanism_note" => "no primary species channel rows available",
            ))
            continue
        end

        states = Dict{Float64, Any}()
        main_by_xi = Dict{Float64, Any}()
        for r in (prev_row, target_row, next_row)
            r === nothing && continue
            xi = fval(r, :xi)
            states[xi] = build_state_from_main_row(r)
            main_by_xi[xi] = r
        end
        if convergence_enabled(opts)
            high_cfg = convergence_tau_config(tau_cfg, opts)
            append_convergence_rows!(
                convergence_rows,
                window_id,
                mode_key,
                cand,
                selected,
                channel_rows,
                states,
                main_by_xi,
                high_cfg,
            )
        end
        opts.only_convergence_gate && continue

        branch_aligned_any = false
        dominant_branch = ""
        dominant_branch_score = -Inf
        max_rate_err = -Inf
        dominant_notes = String[]
        for channel in selected
            process = Symbol(channel)
            direct_rate = haskey(target_map, channel) ? fval(target_map[channel], :rate) : NaN
            branch_score = if haskey(target_map, channel)
                total = fval(target_map[channel], :total)
                contrib = fval(target_map[channel], :contribution)
                isfinite(total) && abs(total) > 0 ? contrib / total : abs(contrib)
            else
                0.0
            end
            for (xi, st) in sort(collect(states); by=x -> x[1])
                summary, bands, samples = try
                    process_branch_scan(window_id, process, st, tau_cfg, opts)
                catch err
                    if closeval(xi, fval(cand, :xi); atol=1.0e-7)
                        push!(dominant_notes, "$(channel):denominator_scan_skipped=$(replace(sprint(showerror, err), ',' => ';'))")
                    end
                    continue
                end
                append!(denom_rows, [merge(summary, Dict{String, Any}("mode_key" => mode_key, "observable" => sval(cand, :observable)))])
                append!(band_rows, bands)
                append!(sample_rows, samples)

                if closeval(xi, fval(cand, :xi); atol=1.0e-7)
                    rate_stats = compute_rate_band_stats(process, st, tau_cfg)
                    for b in 1:(length(BAND_EDGES) - 1)
                        b <= length(rate_stats.omega) || continue
                        push!(band_rows, Dict{String, Any}(
                            "window_id" => window_id,
                            "channel" => channel,
                            "xi" => xi,
                            "band" => b,
                            "ds_left" => BAND_EDGES[b],
                            "ds_right" => BAND_EDGES[b + 1],
                            "prefactor" => rate_stats.prefactor,
                            "omega_bin" => rate_stats.omega[b],
                            "omega_sigma_bin" => rate_stats.omega_sigma[b],
                            "sigma_eff_bin" => rate_stats.omega[b] == 0.0 ? NaN : rate_stats.omega_sigma[b] / rate_stats.omega[b],
                            "rate_bin" => rate_stats.prefactor * rate_stats.omega_sigma[b],
                            "rate_func" => rate_stats.rate,
                        ))
                    end
                    err = isfinite(direct_rate) && abs(direct_rate) > 0 ? abs(rate_stats.rate - direct_rate) / abs(direct_rate) : 0.0
                    max_rate_err = max(max_rate_err, err)
                    branch_aligned_any |= Bool(summary["denominator_sigma_same_band"])
                    if Bool(summary["denominator_sigma_same_band"]) && branch_score > dominant_branch_score
                        dominant_branch = String(summary["dominant_denominator_branch"])
                        dominant_branch_score = branch_score
                    elseif isempty(dominant_branch)
                        dominant_branch = String(summary["dominant_denominator_branch"])
                    end
                    push!(dominant_notes, "$(channel):rate_rel_err=$(round(err; sigdigits=4)),branch=$(summary["dominant_denominator_branch"]),aligned=$(summary["denominator_sigma_same_band"])")
                end
            end
        end

        !isfinite(max_rate_err) && (max_rate_err = Inf)
        verdict, score = classify_window(coverage, max_rate_err, branch_aligned_any, Bool(upstream["upstream_branch_flag"]))
        applicability = sval(cand, :mechanism_applicability)
        extra_note = ""
        if applicability != "tau_species" && verdict == "small_denominator_supported"
            verdict = "small_denominator_supported_indirect"
            score = 0.75
            extra_note = " | composite observable: denominator evidence is through same-point primary-species channel proxy"
        end
        push!(mechanism_rows, Dict{String, Any}(
            "window_id" => window_id,
            "window_label" => sval(cand, :window_label),
            "mode_key" => mode_key,
            "plot_panel" => sval(cand, :plot_panel),
            "plot_series" => sval(cand, :plot_series),
            "observable" => sval(cand, :observable),
            "primary_species" => sval(cand, :primary_species),
            "xi" => fval(cand, :xi),
            "mechanism_verdict" => verdict,
            "evidence_score" => score,
            "dominant_channels" => join(selected, ";"),
            "channel_coverage_share" => coverage,
            "dominant_denominator_branch" => dominant_branch,
            "max_rate_reproduction_rel_error" => max_rate_err,
            "denominator_sigma_alignment" => branch_aligned_any,
            "upstream_branch_flag" => upstream["upstream_branch_flag"],
            "mechanism_note" => join(dominant_notes, " | ") * extra_note,
        ))
    end

    verdict_by_window = Dict(row["window_id"] => row for row in mechanism_rows)
    all_candidates = collect(CSV.File(opts.candidate_csv))
    for cand in all_candidates
        source = sval(cand, :candidate_source)
        occursin("nonmonotonic_turning_point", source) || continue
        window_id = sval(cand, :window_id)
        selected = boolish(getproperty(cand, :selected_for_deep_scan))
        classification = "inconclusive"
        evidence_window = ""
        note = "not selected for denominator-chain scan"
        if selected && haskey(verdict_by_window, window_id)
            classification = String(verdict_by_window[window_id]["mechanism_verdict"])
            evidence_window = window_id
            note = String(verdict_by_window[window_id]["mechanism_note"])
        else
            target_xi = fval(cand, :xi)
            same_point = [row for row in mechanism_rows if row["mode_key"] == sval(cand, :mode_key) &&
                row["plot_panel"] == sval(cand, :plot_panel) &&
                row["plot_series"] == sval(cand, :plot_series) &&
                abs(Float64(row["xi"]) - target_xi) <= 0.0500001 &&
                row["mechanism_verdict"] == "small_denominator_supported"]
            if !isempty(same_point)
                classification = "small_denominator_supported_indirect"
                evidence_window = String(same_point[1]["window_id"])
                note = "shares panel/xi neighborhood with a denominator-supported deep window"
            end
        end
        push!(global_rows, Dict{String, Any}(
            "trend_key" => join([sval(cand, :mode_key), sval(cand, :plot_panel), sval(cand, :plot_series), sval(cand, :observable)], "|"),
            "window_id" => window_id,
            "xi" => fval(cand, :xi),
            "observable" => sval(cand, :observable),
            "mechanism_classification" => classification,
            "evidence_window_id" => evidence_window,
            "classification_note" => note,
        ))
    end

    if opts.only_convergence_gate
        write_csv_rows(joinpath(opts.out_dir, "local_rate_convergence_gate.csv"),
            ["window_id", "mode_key", "primary_species", "channel", "target_xi", "neighbor_xis", "production_target_rate", "production_neighbor_max_rate", "production_spike_ratio", "high_target_rate", "high_neighbor_max_rate", "high_spike_ratio", "target_rel_diff_high_vs_production", "high_to_production_spike_ratio", "convergence_verdict", "high_p_nodes", "high_angle_nodes", "high_phi_nodes", "high_n_sigma_points", "high_sigma_grid_n"],
            convergence_rows)
        manifest = Dict(
            "schema" => "phase_guided_p128_mechanism_scan_manifest_v1",
            "case_name" => opts.case_name,
            "candidate_csv" => relpath(opts.candidate_csv, PROJECT_ROOT),
            "window_count" => length(cands),
            "channel_limit" => opts.channel_limit,
            "coverage_threshold" => opts.coverage_threshold,
            "only_convergence_gate" => opts.only_convergence_gate,
            "integration_mode" => String(opts.integration_mode),
            "production_a_builder_config" => Dict(
                "p_nodes" => PRODUCTION_A_BUILDER_CONFIG.p_nodes,
                "p_max" => PRODUCTION_A_BUILDER_CONFIG.p_max,
                "cos_nodes" => PRODUCTION_A_BUILDER_CONFIG.cos_nodes,
                "use_aniso" => PRODUCTION_A_BUILDER_CONFIG.use_aniso,
            ),
            "convergence_gate_enabled" => convergence_enabled(opts),
            "convergence_p_nodes" => opts.convergence_p_nodes,
            "convergence_angle_nodes" => opts.convergence_angle_nodes,
            "convergence_phi_nodes" => opts.convergence_phi_nodes,
            "convergence_n_sigma_points" => opts.convergence_n_sigma_points,
            "convergence_sigma_grid_n" => opts.convergence_sigma_grid_n,
            "repository_head" => readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`),
            "generator" => relpath(@__FILE__, PROJECT_ROOT),
            "generator_sha256" => bytes2hex(open(sha256, @__FILE__)),
        )
        open(joinpath(opts.out_dir, "mechanism_manifest.json"), "w") do io
            JSON.print(io, manifest, 2)
            println(io)
        end
        println("Wrote local convergence gate to $(relpath(opts.out_dir, PROJECT_ROOT))")
        println("  windows: $(length(cands))")
        println("  convergence rows: $(length(convergence_rows))")
        return
    end

    write_csv_rows(joinpath(opts.out_dir, "mechanism_window_summary.csv"),
        ["window_id", "window_label", "mode_key", "plot_panel", "plot_series", "observable", "primary_species", "xi", "mechanism_verdict", "evidence_score", "dominant_channels", "channel_coverage_share", "dominant_denominator_branch", "max_rate_reproduction_rel_error", "denominator_sigma_alignment", "upstream_branch_flag", "mechanism_note"],
        mechanism_rows)
    write_csv_rows(joinpath(opts.out_dir, "denominator_chain_summary.csv"),
        ["window_id", "mode_key", "observable", "channel", "xi", "sigma_peak_ds", "sigma_peak_value", "simple_peak_ds", "simple_peak_invabs", "mixed_peak_ds", "mixed_peak_invabs", "dominant_denominator_branch", "denominator_sigma_same_band", "max_abs_D_simple_sq", "max_abs_D_mixed_sq", "max_abs_D_total_sq"],
        denom_rows)
    write_csv_rows(joinpath(opts.out_dir, "denominator_chain_band_table.csv"),
        ["window_id", "channel", "xi", "band", "ds_left", "ds_right", "area_invabs_den_simple", "area_invabs_detM", "area_abs_D_total_sq", "area_sigma", "prefactor", "omega_bin", "omega_sigma_bin", "sigma_eff_bin", "rate_bin", "rate_func"],
        band_rows)
    write_csv_rows(joinpath(opts.out_dir, "denominator_ds_samples.csv"),
        ["window_id", "channel", "xi", "ds", "s_value", "den_simple_re", "den_simple_im", "detM_re", "detM_im", "invabs_den_simple", "invabs_detM", "abs_D_simple_sq", "abs_D_mixed_sq", "abs_D_total_sq", "sigma"],
        sample_rows)
    write_csv_rows(joinpath(opts.out_dir, "upstream_branch_smoothness_summary.csv"),
        ["window_id", "n_points", "max_rel_step", "max_rel_curvature", "dominant_background_driver", "upstream_branch_flag"],
        upstream_rows)
    write_csv_rows(joinpath(opts.out_dir, "global_nonmonotonic_mechanism_summary.csv"),
        ["trend_key", "window_id", "xi", "observable", "mechanism_classification", "evidence_window_id", "classification_note"],
        global_rows)
    write_csv_rows(joinpath(opts.out_dir, "local_rate_convergence_gate.csv"),
        ["window_id", "mode_key", "primary_species", "channel", "target_xi", "neighbor_xis", "production_target_rate", "production_neighbor_max_rate", "production_spike_ratio", "high_target_rate", "high_neighbor_max_rate", "high_spike_ratio", "target_rel_diff_high_vs_production", "high_to_production_spike_ratio", "convergence_verdict", "high_p_nodes", "high_angle_nodes", "high_phi_nodes", "high_n_sigma_points", "high_sigma_grid_n"],
        convergence_rows)

    manifest = Dict(
        "schema" => "phase_guided_p128_mechanism_scan_manifest_v1",
        "case_name" => opts.case_name,
        "candidate_csv" => relpath(opts.candidate_csv, PROJECT_ROOT),
        "window_count" => length(cands),
        "channel_limit" => opts.channel_limit,
        "coverage_threshold" => opts.coverage_threshold,
        "ds_grid_n" => opts.ds_grid_n,
        "max_ds" => opts.max_ds,
        "integration_mode" => String(opts.integration_mode),
        "production_a_builder_config" => Dict(
            "p_nodes" => PRODUCTION_A_BUILDER_CONFIG.p_nodes,
            "p_max" => PRODUCTION_A_BUILDER_CONFIG.p_max,
            "cos_nodes" => PRODUCTION_A_BUILDER_CONFIG.cos_nodes,
            "use_aniso" => PRODUCTION_A_BUILDER_CONFIG.use_aniso,
        ),
        "convergence_gate_enabled" => convergence_enabled(opts),
        "convergence_p_nodes" => opts.convergence_p_nodes,
        "convergence_angle_nodes" => opts.convergence_angle_nodes,
        "convergence_phi_nodes" => opts.convergence_phi_nodes,
        "convergence_n_sigma_points" => opts.convergence_n_sigma_points,
        "convergence_sigma_grid_n" => opts.convergence_sigma_grid_n,
        "only_convergence_gate" => opts.only_convergence_gate,
        "repository_head" => readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`),
        "generator" => relpath(@__FILE__, PROJECT_ROOT),
        "generator_sha256" => bytes2hex(open(sha256, @__FILE__)),
    )
    open(joinpath(opts.out_dir, "mechanism_manifest.json"), "w") do io
        JSON.print(io, manifest, 2)
        println(io)
    end
    println("Wrote p128 mechanism tables to $(relpath(opts.out_dir, PROJECT_ROOT))")
    println("  windows: $(length(cands))")
    println("  mechanism rows: $(length(mechanism_rows))")
    println("  convergence rows: $(length(convergence_rows))")
end

main()
