#!/usr/bin/env julia

# Scan |M|^2 with respect to s and t for a given scattering process.
# Supports 2D grid (s×t) and 1D scans (fix s, vary t with physical t-bounds; or fix t, vary s).
# Outputs CSV plus a heatmap/line plot under data/processed.

using CSV
using Plots
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src", "relaxtime"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalPropagator.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScatteringAmplitude: scattering_amplitude_squared
using .TotalPropagator: get_quark_masses_for_process
using .TotalCrossSection: calculate_t_bounds

struct CLIOptions
    mode::Symbol
    process::Symbol
    s_min::Float64
    s_max::Float64
    s_steps::Int
    s_fixed::Float64
    t_min::Float64
    t_max::Float64
    t_steps::Int
    t_fixed::Float64
    use_t_bounds::Bool
    T_MeV::Float64
    mu_u_MeV::Float64
    mu_d_MeV::Float64
    mu_s_MeV::Float64
    m_u_MeV::Float64
    m_d_MeV::Float64
    m_s_MeV::Float64
    phi::Float64
    phibar::Float64
    xi::Float64
    output_csv::String
    output_fig::String
end

function print_usage()
    println("Usage: julia scripts/relaxtime/scan_scattering_amplitude.jl [options]\n")
    println("Options:")
    println("  --mode <grid|scan-t|scan-s> Mode: 2D grid (default), scan t at fixed s, or scan s at fixed t")
    println("  --process <symbol>         Scattering process (default uu_to_uu)")
    println("  --s-min <float>            Minimum s in fm^-2 (default 10.0)")
    println("  --s-max <float>            Maximum s in fm^-2 (default 20.0)")
    println("  --s-steps <int>            Number of s samples (default 40)")
    println("  --s-fixed <float>          s value for scan-t mode (default 15.0 fm^-2)")
    println("  --t-min <float>            Minimum t in fm^-2 (default -1.0)")
    println("  --t-max <float>            Maximum t in fm^-2 (default -0.05)")
    println("  --t-steps <int>            Number of t samples (default 40)")
    println("  --t-fixed <float>          t value for scan-s mode (default -0.3 fm^-2)")
    println("  --use-t-bounds             In scan-t mode, derive t range from physical bounds (default true)")
    println("  --no-use-t-bounds          Disable physical t bounds and use provided t-min/max")
    println("  --T-MeV <float>            Temperature in MeV (default 150)")
    println("  --mu-u <float>             Mu_u in MeV (default 0)")
    println("  --mu-d <float>             Mu_d in MeV (default 0)")
    println("  --mu-s <float>             Mu_s in MeV (default 0)")
    println("  --m-u <float>              Constituent m_u in MeV (default 300)")
    println("  --m-d <float>              Constituent m_d in MeV (default 300)")
    println("  --m-s <float>              Constituent m_s in MeV (default 500)")
    println("  --phi <float>             Polyakov loop Phi (default 0.5)")
    println("  --phibar <float>          Conjugate Polyakov loop Phi_bar (default 0.5)")
    println("  --xi <float>              Momentum anisotropy xi (default 0.0)")
    println("  --output-csv <path>       CSV output path (default data/processed/results/relaxtime/scattering_amplitude_surface.csv)")
    println("  --output-fig <path>       Figure output path (default data/processed/figures/relaxtime/scattering_amplitude_surface.png)")
    println("  -h, --help                Show this help message")
end

function parse_args(args)::CLIOptions
    opts = Dict{Symbol, Any}(
        :mode => "grid",
        :process => "uu_to_uu",
        :s_min => 10.0,
        :s_max => 20.0,
        :s_steps => 40,
        :s_fixed => 15.0,
        :t_min => -1.0,
        :t_max => -0.05,
        :t_steps => 40,
        :t_fixed => -0.3,
        :use_t_bounds => true,
        :T_MeV => 150.0,
        :mu_u_MeV => 0.0,
        :mu_d_MeV => 0.0,
        :mu_s_MeV => 0.0,
        :m_u_MeV => 300.0,
        :m_d_MeV => 300.0,
        :m_s_MeV => 500.0,
        :phi => 0.5,
        :phibar => 0.5,
        :xi => 0.0,
        :output_csv => joinpath("data", "processed", "results", "relaxtime", "scattering_amplitude_surface.csv"),
        :output_fig => joinpath("data", "processed", "figures", "relaxtime", "scattering_amplitude_surface.png"),
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

        if arg == "--process"
            opts[:process] = require_value()
        elseif arg == "--mode"
            opts[:mode] = require_value()
        elseif arg == "--s-min"
            opts[:s_min] = parse(Float64, require_value())
        elseif arg == "--s-max"
            opts[:s_max] = parse(Float64, require_value())
        elseif arg == "--s-steps"
            opts[:s_steps] = parse(Int, require_value())
        elseif arg == "--s-fixed"
            opts[:s_fixed] = parse(Float64, require_value())
        elseif arg == "--t-min"
            opts[:t_min] = parse(Float64, require_value())
        elseif arg == "--t-max"
            opts[:t_max] = parse(Float64, require_value())
        elseif arg == "--t-steps"
            opts[:t_steps] = parse(Int, require_value())
        elseif arg == "--t-fixed"
            opts[:t_fixed] = parse(Float64, require_value())
        elseif arg == "--use-t-bounds"
            opts[:use_t_bounds] = true
        elseif arg == "--no-use-t-bounds"
            opts[:use_t_bounds] = false
        elseif arg == "--T-MeV"
            opts[:T_MeV] = parse(Float64, require_value())
        elseif arg == "--mu-u"
            opts[:mu_u_MeV] = parse(Float64, require_value())
        elseif arg == "--mu-d"
            opts[:mu_d_MeV] = parse(Float64, require_value())
        elseif arg == "--mu-s"
            opts[:mu_s_MeV] = parse(Float64, require_value())
        elseif arg == "--m-u"
            opts[:m_u_MeV] = parse(Float64, require_value())
        elseif arg == "--m-d"
            opts[:m_d_MeV] = parse(Float64, require_value())
        elseif arg == "--m-s"
            opts[:m_s_MeV] = parse(Float64, require_value())
        elseif arg == "--phi"
            opts[:phi] = parse(Float64, require_value())
        elseif arg == "--phibar"
            opts[:phibar] = parse(Float64, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--output-csv"
            opts[:output_csv] = require_value()
        elseif arg == "--output-fig"
            opts[:output_fig] = require_value()
        elseif arg == "-h" || arg == "--help"
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    mode_sym = Symbol(replace(String(opts[:mode]), "-" => "_"))

    return CLIOptions(
        mode_sym,
        Symbol(opts[:process]),
        Float64(opts[:s_min]),
        Float64(opts[:s_max]),
        Int(opts[:s_steps]),
        Float64(opts[:s_fixed]),
        Float64(opts[:t_min]),
        Float64(opts[:t_max]),
        Int(opts[:t_steps]),
        Float64(opts[:t_fixed]),
        Bool(opts[:use_t_bounds]),
        Float64(opts[:T_MeV]),
        Float64(opts[:mu_u_MeV]),
        Float64(opts[:mu_d_MeV]),
        Float64(opts[:mu_s_MeV]),
        Float64(opts[:m_u_MeV]),
        Float64(opts[:m_d_MeV]),
        Float64(opts[:m_s_MeV]),
        Float64(opts[:phi]),
        Float64(opts[:phibar]),
        Float64(opts[:xi]),
        String(opts[:output_csv]),
        String(opts[:output_fig]),
    )
end

function build_physical_parameters(opts::CLIOptions)
    T = opts.T_MeV / ħc_MeV_fm
    μ_u = opts.mu_u_MeV / ħc_MeV_fm
    μ_d = opts.mu_d_MeV / ħc_MeV_fm
    μ_s = opts.mu_s_MeV / ħc_MeV_fm
    m_u = opts.m_u_MeV / ħc_MeV_fm
    m_d = opts.m_d_MeV / ħc_MeV_fm
    m_s = opts.m_s_MeV / ħc_MeV_fm

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS

    A_u = A(m_u, μ_u, T, opts.phi, opts.phibar, nodes_p, weights_p)
    A_d = A(m_d, μ_d, T, opts.phi, opts.phibar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, opts.phi, opts.phibar, nodes_p, weights_p)

    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)

    quark_params = (
        m = (u=m_u, d=m_d, s=m_s),
        μ = (u=μ_u, d=μ_d, s=μ_s),
        A = (u=A_u, d=A_d, s=A_s),
    )
    thermo_params = (T=T, Φ=opts.phi, Φbar=opts.phibar, ξ=opts.xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function determine_t_range(opts::CLIOptions, params; s_val::Float64)
    if opts.use_t_bounds
        m1, m2, m3, m4 = get_quark_masses_for_process(opts.process, params.quark_params)
        bounds = calculate_t_bounds(s_val, m1, m2, m3, m4)
        return collect(range(bounds.t_min, bounds.t_max; length=opts.t_steps))
    else
        return collect(range(opts.t_min, opts.t_max; length=opts.t_steps))
    end
end

function compute_surface(opts::CLIOptions, params)
    s_values = collect(range(opts.s_min, opts.s_max; length=opts.s_steps))
    t_values = collect(range(opts.t_min, opts.t_max; length=opts.t_steps))

    grid = fill(NaN, length(t_values), length(s_values))
    rows = Vector{NamedTuple{(:process, :s, :t, :M_squared), Tuple{String, Float64, Float64, Float64}}}()

    for (i_s, s) in enumerate(s_values)
        for (i_t, t) in enumerate(t_values)
            M_val = try
                scattering_amplitude_squared(opts.process, s, t, params.quark_params, params.thermo_params, params.K_coeffs)
            catch err
                @warn "Skip point" process=opts.process s=s t=t error=err
                NaN
            end
            grid[i_t, i_s] = M_val
            push!(rows, (process=String(opts.process), s=s, t=t, M_squared=M_val))
        end
    end

    return (s_values=s_values, t_values=t_values, grid=grid, rows=rows)
end

function compute_scan_t(opts::CLIOptions, params)
    s_val = opts.s_fixed
    t_values = determine_t_range(opts, params; s_val=s_val)
    rows = Vector{NamedTuple{(:process, :s, :t, :M_squared), Tuple{String, Float64, Float64, Float64}}}()
    m_vals = Float64[]

    for t in t_values
        M_val = try
            scattering_amplitude_squared(opts.process, s_val, t, params.quark_params, params.thermo_params, params.K_coeffs)
        catch err
            @warn "Skip point" process=opts.process s=s_val t=t error=err
            NaN
        end
        push!(m_vals, M_val)
        push!(rows, (process=String(opts.process), s=s_val, t=t, M_squared=M_val))
    end

    return (mode=:scan_t, s_val=s_val, t_values=t_values, m_vals=m_vals, rows=rows)
end

function compute_scan_s(opts::CLIOptions, params)
    s_values = collect(range(opts.s_min, opts.s_max; length=opts.s_steps))
    t_val = opts.t_fixed
    rows = Vector{NamedTuple{(:process, :s, :t, :M_squared), Tuple{String, Float64, Float64, Float64}}}()
    m_vals = Float64[]

    for s in s_values
        M_val = try
            scattering_amplitude_squared(opts.process, s, t_val, params.quark_params, params.thermo_params, params.K_coeffs)
        catch err
            @warn "Skip point" process=opts.process s=s t=t_val error=err
            NaN
        end
        push!(m_vals, M_val)
        push!(rows, (process=String(opts.process), s=s, t=t_val, M_squared=M_val))
    end

    return (mode=:scan_s, s_values=s_values, t_val=t_val, m_vals=m_vals, rows=rows)
end

function save_outputs(opts::CLIOptions, surface)
    mkpath(dirname(opts.output_csv))
    mkpath(dirname(opts.output_fig))
    CSV.write(opts.output_csv, surface.rows)

    if haskey(surface, :grid)
        if all(isnan, surface.grid)
            @warn "All grid values are NaN; skip plotting" output=opts.output_fig
            return
        end
        plt = heatmap(
            surface.s_values,
            surface.t_values,
            surface.grid;
            xlabel="s [fm^-2]",
            ylabel="t [fm^-2]",
            title=string("|M|^2 for ", opts.process),
            colorbar_title="fm^-4",
        )
        savefig(plt, opts.output_fig)
    elseif surface.mode == :scan_t
        plt = plot(
            surface.t_values,
            surface.m_vals;
            xlabel="t [fm^-2]",
            ylabel="|M|^2 [fm^-4]",
            title=@sprintf("|M|^2 vs t @ s=%.3f (%s)", surface.s_val, opts.process),
        )
        savefig(plt, opts.output_fig)
    elseif surface.mode == :scan_s
        plt = plot(
            surface.s_values,
            surface.m_vals;
            xlabel="s [fm^-2]",
            ylabel="|M|^2 [fm^-4]",
            title=@sprintf("|M|^2 vs s @ t=%.3f (%s)", surface.t_val, opts.process),
        )
        savefig(plt, opts.output_fig)
    end
end

function main(args)
    opts = parse_args(args)
    params = build_physical_parameters(opts)
    surface = if opts.mode == :scan_t
        compute_scan_t(opts, params)
    elseif opts.mode == :scan_s
        compute_scan_s(opts, params)
    else
        compute_surface(opts, params)
    end
    save_outputs(opts, surface)
    println(@sprintf("Done. CSV: %s | Figure: %s", opts.output_csv, opts.output_fig))
end

main(ARGS)
