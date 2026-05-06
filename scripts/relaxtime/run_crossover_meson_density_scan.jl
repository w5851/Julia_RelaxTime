"""
沿 crossover line 运行 charged/aggregate meson-density workflow 扫描。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_crossover_meson_density_scan

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_crossover_meson_density_scan.jl [options]\n")
    println("Options:")
    println("  --mu-min <MeV>                      quark chemical potential lower bound (default 10)")
    println("  --mu-max <MeV>                      quark chemical potential upper bound (required)")
    println("  --T-min <MeV>                       crossover search T lower bound (required)")
    println("  --T-max <MeV>                       crossover search T upper bound (required)")
    println("  --n-mu <int>                        mu sampling count for crossover line (default 12)")
    println("  --xi <float>                        anisotropy xi (default 0)")
    println("  --method <peak|inflection>          crossover locator (default peak)")
    println("  --variable <symbol>                 crossover variable (default phi_u)")
    println("  --flavor-chemical-profile <name>    flavor chemical profile (default default)")
    println("  --meson-chemical-profile <name>     meson chemical profile (default default)")
    println("  --regime <stable|phase_shift_current|phase_shift_gbu|strict_bw_stage1|strict_bw_stage2>")
    println("  --output <path>                     output csv path")
    println("  --overwrite <true|false>            overwrite existing output (default true)")
    println("  --p-num <int>                       equilibrium momentum nodes (default 24)")
    println("  --t-num <int>                       equilibrium theta nodes (default 8)")
    println("  --max-iter <int>                    solver / meson iterations (default 40)")
    println("  --qmax <float>                      q upper bound for BW/BU kernels (default 12)")
    println("  --q-nodes <int>                     q nodes for BW/BU kernels (default 48)")
    println("  --omega-min <float>                 phase-shift omega lower bound (default 0.05)")
    println("  --omega-max <float>                 omega upper bound (default 10)")
    println("  --omega-nodes <int>                 omega nodes (default 48)")
    println("  --eta <float>                       phase-shift eta (default 1e-6)")
    println("  -h, --help                          show help")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :mu_min_MeV => 10.0,
        :mu_max_MeV => nothing,
        :T_min_MeV => nothing,
        :T_max_MeV => nothing,
        :n_mu => 12,
        :xi => 0.0,
        :method => :peak,
        :variable => :phi_u,
        :flavor_chemical_profile_name => "default",
        :meson_chemical_profile_name => "default",
        :regime => :stable,
        :output_path => nothing,
        :overwrite => true,
        :p_num => 24,
        :t_num => 8,
        :max_iter => 40,
        :qmax => 12.0,
        :q_nodes => 48,
        :omega_min => 0.05,
        :omega_max => 10.0,
        :omega_nodes => 48,
        :eta => 1e-6,
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
        if arg == "--mu-min"
            opts[:mu_min_MeV] = parse(Float64, require_value())
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
        elseif arg == "--flavor-chemical-profile"
            opts[:flavor_chemical_profile_name] = require_value()
        elseif arg == "--meson-chemical-profile"
            opts[:meson_chemical_profile_name] = require_value()
        elseif arg == "--regime"
            opts[:regime] = Symbol(require_value())
        elseif arg == "--output"
            opts[:output_path] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = lowercase(require_value()) == "true"
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
        elseif arg == "--eta"
            opts[:eta] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    for key in (:mu_max_MeV, :T_min_MeV, :T_max_MeV)
        isnothing(opts[key]) && error("missing required option: $(key)")
    end
    return opts
end

function main(args::Vector{String})
    opts = parse_args(args)
    result = run_crossover_meson_density_scan(
        mu_min_MeV=Float64(opts[:mu_min_MeV]),
        mu_max_MeV=Float64(opts[:mu_max_MeV]),
        T_min_MeV=Float64(opts[:T_min_MeV]),
        T_max_MeV=Float64(opts[:T_max_MeV]),
        n_mu=Int(opts[:n_mu]),
        xi=Float64(opts[:xi]),
        method=Symbol(opts[:method]),
        variable=Symbol(opts[:variable]),
        flavor_chemical_profile_name=String(opts[:flavor_chemical_profile_name]),
        meson_chemical_profile_name=String(opts[:meson_chemical_profile_name]),
        regime=Symbol(opts[:regime]),
        output_path=isnothing(opts[:output_path]) ? Models.CrossoverMesonDensityScan.DEFAULT_CROSSOVER_MESON_DENSITY_OUTPUT_PATH : String(opts[:output_path]),
        overwrite=Bool(opts[:overwrite]),
        p_num=Int(opts[:p_num]),
        t_num=Int(opts[:t_num]),
        max_iter=Int(opts[:max_iter]),
        strict_bw_qmax=Float64(opts[:qmax]),
        strict_bw_q_nodes=Int(opts[:q_nodes]),
        strict_bw_omega_max=Float64(opts[:omega_max]),
        strict_bw_omega_nodes=Int(opts[:omega_nodes]),
        phase_shift_qmax=Float64(opts[:qmax]),
        phase_shift_q_nodes=Int(opts[:q_nodes]),
        phase_shift_omega_min=Float64(opts[:omega_min]),
        phase_shift_omega_max=Float64(opts[:omega_max]),
        phase_shift_omega_nodes=Int(opts[:omega_nodes]),
        phase_shift_eta=Float64(opts[:eta]),
    )

    println("crossover meson-density scan completed")
    println("output_path=$(result.output_path)")
    println("points=$(result.points)")
    println("workflow_entry=$(result.workflow_entry)")
    return result
end

abspath(PROGRAM_FILE) == @__FILE__ && main(ARGS)
