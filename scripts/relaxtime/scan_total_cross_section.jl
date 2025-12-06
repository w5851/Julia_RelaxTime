#!/usr/bin/env julia

# Scan total cross section σ(s) for a given scattering process.
# Uses TotalCrossSection.total_cross_section (t integral included) and outputs CSV + line plot.

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
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section, DEFAULT_T_INTEGRAL_POINTS

struct CLIOptions
    process::Symbol
    s_min::Float64
    s_max::Float64
    s_steps::Int
    n_points::Int
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
    println("Usage: julia scripts/relaxtime/scan_total_cross_section.jl [options]\n")
    println("Options:")
    println("  --process <symbol>   Scattering process (default uu_to_uu)")
    println("  --s-min <float>      Minimum s in fm^-2 (default 10.0)")
    println("  --s-max <float>      Maximum s in fm^-2 (default 20.0)")
    println("  --s-steps <int>      Number of s samples (default 40)")
    println("  --n-points <int>     Gauss-Legendre t-integration points (default 32)")
    println("  --T-MeV <float>      Temperature in MeV (default 150)")
    println("  --mu-u <float>       Mu_u in MeV (default 0)")
    println("  --mu-d <float>       Mu_d in MeV (default 0)")
    println("  --mu-s <float>       Mu_s in MeV (default 0)")
    println("  --m-u <float>        Constituent m_u in MeV (default 300)")
    println("  --m-d <float>        Constituent m_d in MeV (default 300)")
    println("  --m-s <float>        Constituent m_s in MeV (default 500)")
    println("  --phi <float>       Polyakov loop Phi (default 0.5)")
    println("  --phibar <float>    Conjugate Polyakov loop Phi_bar (default 0.5)")
    println("  --xi <float>        Momentum anisotropy xi (default 0.0)")
    println("  --output-csv <path> CSV output (default data/processed/results/relaxtime/total_cross_section_scan.csv)")
    println("  --output-fig <path> Figure output (default data/processed/figures/relaxtime/total_cross_section_scan.png)")
    println("  -h, --help          Show this help message")
end

function parse_args(args)::CLIOptions
    opts = Dict{Symbol, Any}(
        :process => "uu_to_uu",
        :s_min => 10.0,
        :s_max => 20.0,
        :s_steps => 40,
        :n_points => DEFAULT_T_INTEGRAL_POINTS,
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
        :output_csv => joinpath("data", "processed", "results", "relaxtime", "total_cross_section_scan.csv"),
        :output_fig => joinpath("data", "processed", "figures", "relaxtime", "total_cross_section_scan.png"),
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
        elseif arg == "--s-min"
            opts[:s_min] = parse(Float64, require_value())
        elseif arg == "--s-max"
            opts[:s_max] = parse(Float64, require_value())
        elseif arg == "--s-steps"
            opts[:s_steps] = parse(Int, require_value())
        elseif arg == "--n-points"
            opts[:n_points] = parse(Int, require_value())
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

    return CLIOptions(
        Symbol(opts[:process]),
        Float64(opts[:s_min]),
        Float64(opts[:s_max]),
        Int(opts[:s_steps]),
        Int(opts[:n_points]),
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

    G_u = calculate_G_from_A(A_u)
    G_s = calculate_G_from_A(A_s)

    quark_params = (
        m = (u=m_u, d=m_d, s=m_s),
        μ = (u=μ_u, d=μ_d, s=μ_s),
        A = (u=A_u, d=A_d, s=A_s),
    )
    thermo_params = (T=T, Φ=opts.phi, Φbar=opts.phibar, ξ=opts.xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function compute_scan(opts::CLIOptions, params)
    s_values = collect(range(opts.s_min, opts.s_max; length=opts.s_steps))
    rows = Vector{NamedTuple{(:process, :s, :sigma), Tuple{String, Float64, Float64}}}()
    sigma_vals = Float64[]

    for s in s_values
        σ = try
            total_cross_section(opts.process, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=opts.n_points)
        catch err
            @warn "Skip point" process=opts.process s=s error=err
            NaN
        end
        push!(sigma_vals, σ)
        push!(rows, (process=String(opts.process), s=s, sigma=σ))
    end

    return (s_values=s_values, sigma_vals=sigma_vals, rows=rows)
end

function save_outputs(opts::CLIOptions, data)
    mkpath(dirname(opts.output_csv))
    mkpath(dirname(opts.output_fig))
    CSV.write(opts.output_csv, data.rows)

    if all(isnan, data.sigma_vals)
        @warn "All σ values are NaN; skip plotting" output=opts.output_fig
        return
    end

    plt = plot(
        data.s_values,
        data.sigma_vals;
        xlabel="s [fm^-2]",
        ylabel="σ [fm^2]",
        title=@sprintf("σ(s) for %s", opts.process),
    )
    savefig(plt, opts.output_fig)
end

function main(args)
    opts = parse_args(args)
    params = build_physical_parameters(opts)
    data = compute_scan(opts, params)
    save_outputs(opts, data)
    println(@sprintf("Done. CSV: %s | Figure: %s", opts.output_csv, opts.output_fig))
end

main(ARGS)
