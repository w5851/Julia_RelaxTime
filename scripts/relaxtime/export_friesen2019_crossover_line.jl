"""
Export crossover line under a selected PNJL/physics profile pair.

Primary use:
- produce the Friesen-2019-like chiral crossover line as an explicit artifact
- separate path construction from meson-density follow-up scans
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :pnjl_profile => "friesen2019_poly",
        :physics_profile => "friesen2019_base",
        :mu_max_MeV => 600.0,
        :T_min_MeV => 1.0,
        :T_max_MeV => 220.0,
        :n_mu => 24,
        :xi => 0.0,
        :method => :peak,
        :variable => :phi_u,
        :output_dir => joinpath(
            PROJECT_ROOT,
            "data", "outputs", "results", "relaxtime", "meson_density", "crossover_validation",
            "friesen2019_crossover_line",
        ),
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

        if arg == "--pnjl-profile"
            opts[:pnjl_profile] = require_value()
        elseif arg == "--physics-profile"
            opts[:physics_profile] = require_value()
        elseif arg == "--mu-max"
            opts[:mu_max_MeV] = parse(Float64, require_value())
        elseif arg == "--T-min"
            opts[:T_min_MeV] = parse(Float64, require_value())
        elseif arg == "--T-max"
            opts[:T_max_MeV] = parse(Float64, require_value())
        elseif arg == "--n-mu"
            opts[:n_mu] = parse(Int, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--method"
            opts[:method] = Symbol(require_value())
        elseif arg == "--variable"
            opts[:variable] = Symbol(require_value())
        elseif arg == "--output-dir"
            opts[:output_dir] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = lowercase(require_value()) == "true"
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/relaxtime/export_friesen2019_crossover_line.jl [options]")
            println("  --pnjl-profile <name>       model profile (default friesen2019_poly)")
            println("  --physics-profile <name>    physics profile (default friesen2019_base)")
            println("  --mu-max <MeV>              max quark chemical potential (default 600)")
            println("  --T-min <MeV>               min search temperature (default 1)")
            println("  --T-max <MeV>               max search temperature (default 220)")
            println("  --n-mu <int>                number of mu samples (default 24)")
            println("  --xi <float>                anisotropy xi (default 0)")
            println("  --method <peak|inflection>  crossover locator (default peak)")
            println("  --variable <symbol>         crossover variable (default phi_u)")
            println("  --output-dir <path>         output directory")
            println("  --overwrite <true|false>    overwrite output dir (default false)")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return opts
end

function write_crossover_csv(path::String, rows)
    open(path, "w") do io
        println(io, "mu_MeV,muB_MeV,T_crossover_MeV,T_over_muB,rho,method,converged,derivative,variable")
        for row in rows
            muq = Float64(row.mu_MeV)
            muB = 3.0 * muq
            T_mev = Float64(row.T_crossover_MeV)
            x = muB > 0.0 ? T_mev / muB : NaN
            println(io, "$(muq),$(muB),$(T_mev),$(x),$(row.rho),$(row.method),$(row.converged),$(row.derivative),$(row.variable)")
        end
    end
end

function write_summary(path::String, rows, opts)
    valid = [
        row for row in rows
        if Bool(row.converged) && isfinite(Float64(row.mu_MeV)) && isfinite(Float64(row.T_crossover_MeV))
    ]
    open(path, "w") do io
        println(io, "# Friesen 2019 crossover line export")
        println(io)
        println(io, "- pnjl_profile: $(opts[:pnjl_profile])")
        println(io, "- physics_profile: $(opts[:physics_profile])")
        println(io, "- method: $(opts[:method])")
        println(io, "- variable: $(opts[:variable])")
        println(io, "- xi: $(opts[:xi])")
        println(io, "- total_points: $(length(rows))")
        println(io, "- converged_points: $(length(valid))")
        if !isempty(valid)
            first_row = first(valid)
            last_row = last(valid)
            println(io, "- first_valid_point: mu_B=$(3 * first_row.mu_MeV) MeV, T_c=$(first_row.T_crossover_MeV) MeV")
            println(io, "- last_valid_point: mu_B=$(3 * last_row.mu_MeV) MeV, T_c=$(last_row.T_crossover_MeV) MeV")
            println(io, "- T_c(mu_B=0): $(rows[1].T_crossover_MeV) MeV")
        end
    end
end

function main(args::Vector{String})
    opts = parse_args(args)

    ENV["PNJL_PARAM_PROFILE"] = String(opts[:pnjl_profile])
    ENV["PHYSICS_PARAM_PROFILE"] = String(opts[:physics_profile])

    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

    outdir = abspath(String(opts[:output_dir]))
    if isdir(outdir)
        Bool(opts[:overwrite]) || error("output_dir already exists: $outdir")
    end
    mkpath(outdir)

    rows = Base.invokelatest(
        Main.Models.build_crossover_line;
        mu_max_MeV=Float64(opts[:mu_max_MeV]),
        T_min_MeV=Float64(opts[:T_min_MeV]),
        T_max_MeV=Float64(opts[:T_max_MeV]),
        xi=Float64(opts[:xi]),
        n_mu=Int(opts[:n_mu]),
        method=Symbol(opts[:method]),
        variable=Symbol(opts[:variable]),
    )

    csv_path = joinpath(outdir, "crossover_line.csv")
    write_crossover_csv(csv_path, rows)
    write_summary(joinpath(outdir, "README.md"), rows, opts)

    println("friesen crossover line export completed")
    println("output_dir=$(outdir)")
    println("csv_path=$(csv_path)")
    println("points=$(length(rows))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
