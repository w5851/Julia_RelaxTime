#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.TrhoSeedChain
const ħc_MeV_fm = Main.Constants_PNJL.ħc_MeV_fm

struct Options
    chain_dir::String
    output::String
    rho_min::Float64
    rho_max::Float64
    stride::Int
    xi_filter::Union{Nothing, Float64}
    xi_tol::Float64
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :chain_dir => TrhoSeedChain.DEFAULT_CHAIN_DIR,
        :output => joinpath(PROJECT_ROOT, "data", "raw", "pnjl", "seeds", "trho_seed_table.csv"),
        :rho_min => 0.0,
        :rho_max => 3.0,
        :stride => 1,
        :xi_filter => nothing,
        :xi_tol => 1e-6,
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
        if arg == "--chain-dir"
            opts[:chain_dir] = require_value()
        elseif arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--rho-min"
            opts[:rho_min] = parse(Float64, require_value())
        elseif arg == "--rho-max"
            opts[:rho_max] = parse(Float64, require_value())
        elseif arg == "--stride"
            value = parse(Int, require_value())
            value > 0 || error("--stride must be positive")
            opts[:stride] = value
        elseif arg == "--xi"
            opts[:xi_filter] = parse(Float64, require_value())
        elseif arg == "--xi-tol"
            opts[:xi_tol] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    opts[:rho_min] <= opts[:rho_max] || error("rho_min must be <= rho_max")
    return Options(
        String(opts[:chain_dir]),
        String(opts[:output]),
        Float64(opts[:rho_min]),
        Float64(opts[:rho_max]),
        Int(opts[:stride]),
        opts[:xi_filter],
        Float64(opts[:xi_tol]),
    )
end

function print_usage()
    println("Usage: julia scripts/pnjl/export_trho_seeds.jl [options]\n")
    println("Options:")
    println("  --chain-dir <path>      Directory containing Trho chain cache (default data/raw/pnjl/seeds/trho_continuation)")
    println("  --output <path>         Destination CSV path (default data/raw/pnjl/seeds/trho_seed_table.csv)")
    println("  --rho-min <value>       Minimum ρ/ρ₀ to include (default 0.0)")
    println("  --rho-max <value>       Maximum ρ/ρ₀ to include (default 3.0)")
    println("  --stride <int>          Keep every N-th sample along each chain (default 1)")
    println("  --xi <value>            Only export chains whose ξ matches the value (optional)")
    println("  --xi-tol <value>        Allowed deviation when filtering by ξ (default 1e-6)")
    println("  -h, --help              Show this help message")
end

function collect_rows(opts::Options)
    rows = Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}[]
    files = TrhoSeedChain.list_chain_files(opts.chain_dir)
    for path in files
        chain = TrhoSeedChain.load_chain(path)
        if opts.xi_filter !== nothing && abs(chain.xi - opts.xi_filter) > opts.xi_tol
            continue
        end
        for (idx, rho) in enumerate(chain.rho_values)
            rho < opts.rho_min && continue
            rho > opts.rho_max && continue
            if opts.stride > 1 && ((idx - 1) % opts.stride != 0)
                continue
            end
            state = chain.states[idx]
            phi_u, phi_d, phi_s = state[1:3]
            Phi1, Phi2 = state[4], state[5]
            mu_u = state[6] * ħc_MeV_fm
            mu_d = state[7] * ħc_MeV_fm
            mu_s = state[8] * ħc_MeV_fm
            mu_avg = (mu_u + mu_d + mu_s) / 3.0
            push!(rows, (chain.T, mu_avg, rho, chain.xi, phi_u, phi_d, phi_s, Phi1, Phi2, mu_u, mu_d, mu_s))
        end
    end
    sort!(rows; by = x -> (x[1], x[3]))
    return rows
end

const HEADER = (
    "T_MeV",
    "mu_MeV",
    "rho",
    "xi",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi1",
    "Phi2",
    "mu_u_MeV",
    "mu_d_MeV",
    "mu_s_MeV",
)

function write_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(HEADER, ','))
        for row in rows
            println(io, join(row, ','))
        end
    end
end

function main()
    opts = parse_args(copy(ARGS))
    rows = collect_rows(opts)
    @info "Exporting seeds" count=length(rows) output=opts.output
    write_csv(opts.output, rows)
end

main()
