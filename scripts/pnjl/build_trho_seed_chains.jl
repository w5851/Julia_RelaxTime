#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.TrhoSeedChain

struct Options
    T_min::Float64
    T_max::Float64
    T_step::Float64
    xi_values::Vector{Float64}
    output_dir::String
    rho_max::Float64
    coarse_step::Float64
    medium_switch::Float64
    medium_step::Float64
    fine_switch::Float64
    fine_step::Float64
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :T_min => 50.0,
        :T_max => 200.0,
        :T_step => 10.0,
        :xi_values => "0.0",
        :output_dir => TrhoSeedChain.DEFAULT_CHAIN_DIR,
        :rho_max => 3.0,
        :coarse_step => 0.05,
        :medium_switch => 0.5,
        :medium_step => 0.02,
        :fine_switch => 0.3,
        :fine_step => 0.01,
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
        if arg == "--T-min"
            opts[:T_min] = parse(Float64, require_value())
        elseif arg == "--T-max"
            opts[:T_max] = parse(Float64, require_value())
        elseif arg == "--T-step"
            opts[:T_step] = parse(Float64, require_value())
        elseif arg == "--xi-values"
            opts[:xi_values] = require_value()
        elseif arg == "--output-dir"
            opts[:output_dir] = require_value()
        elseif arg == "--rho-max"
            opts[:rho_max] = parse(Float64, require_value())
        elseif arg == "--coarse-step"
            opts[:coarse_step] = parse(Float64, require_value())
        elseif arg == "--medium-switch"
            opts[:medium_switch] = parse(Float64, require_value())
        elseif arg == "--medium-step"
            opts[:medium_step] = parse(Float64, require_value())
        elseif arg == "--fine-switch"
            opts[:fine_switch] = parse(Float64, require_value())
        elseif arg == "--fine-step"
            opts[:fine_step] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    T_min = Float64(opts[:T_min])
    T_max = Float64(opts[:T_max])
    T_step = Float64(opts[:T_step])
    T_min <= T_max || error("T_min must be <= T_max")
    T_step > 0 || error("T_step must be positive")
    xi_values = _parse_float_list(String(opts[:xi_values]))
    return Options(
        T_min,
        T_max,
        T_step,
        xi_values,
        String(opts[:output_dir]),
        Float64(opts[:rho_max]),
        Float64(opts[:coarse_step]),
        Float64(opts[:medium_switch]),
        Float64(opts[:medium_step]),
        Float64(opts[:fine_switch]),
        Float64(opts[:fine_step]),
    )
end

function print_usage()
    println("Usage: julia scripts/pnjl/build_trho_seed_chains.jl [options]\n")
    println("Options:")
    println("  --T-min <value>          Minimum temperature in MeV (default 50)")
    println("  --T-max <value>          Maximum temperature in MeV (default 200)")
    println("  --T-step <value>         Temperature step in MeV (default 10)")
    println("  --xi-values <list>       Comma separated ξ list (default 0.0)")
    println("  --output-dir <path>      Destination directory for chain .jld2 files")
    println("  --rho-max <value>        Maximum ρ/ρ₀ for the dense grid (default 3.0)")
    println("  --coarse-step <value>    Step size for high-density region (default 0.05)")
    println("  --medium-switch <value>  ρ threshold to switch to medium step (default 0.5)")
    println("  --medium-step <value>    Step size between switch and fine region (default 0.02)")
    println("  --fine-switch <value>    ρ threshold to switch to fine step (default 0.3)")
    println("  --fine-step <value>      Step size below fine_switch (default 0.01)")
    println("  -h, --help               Show this message")
end

function _parse_float_list(str::AbstractString)
    parts = split(str, ',')
    values = Float64[]
    for part in parts
        token = strip(part)
        isempty(token) && continue
        push!(values, parse(Float64, token))
    end
    isempty(values) && error("xi-values list cannot be empty")
    return values
end

function build_grid(opts::Options)
    return TrhoSeedChain.build_dense_rho_grid(
        ; rho_max = opts.rho_max,
          coarse_step = opts.coarse_step,
          medium_switch = opts.medium_switch,
          medium_step = opts.medium_step,
          fine_switch = opts.fine_switch,
          fine_step = opts.fine_step,
    )
end

function main()
    opts = parse_args(copy(ARGS))
    T_values = collect(opts.T_min:opts.T_step:opts.T_max)
    rho_grid = build_grid(opts)
    chains = TrhoSeedChain.build_seed_chains(
        ; T_values = T_values,
          xi_values = opts.xi_values,
          rho_grid = rho_grid,
          output_dir = opts.output_dir,
    )
    total_samples = sum(length(chain.rho_values) for chain in chains)
    @info "Seed chains built" chains = length(chains) samples = total_samples dir = opts.output_dir
end

main()
