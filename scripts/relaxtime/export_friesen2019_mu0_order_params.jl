"""
Export mu_B = 0 order-parameter curves for the Friesen-2019-aligned profile pair.

Artifacts:
- mu0_order_params.csv
- README.md
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :pnjl_profile => "friesen2019_poly",
        :physics_profile => "friesen2019_base",
        :T_min_MeV => 1.0,
        :T_max_MeV => 350.0,
        :T_step_MeV => 5.0,
        :xi => 0.0,
        :p_num => 12,
        :t_num => 4,
        :output_dir => joinpath(
            PROJECT_ROOT,
            "data", "outputs", "results", "relaxtime", "meson_density", "crossover_validation",
            "friesen2019_mu0_order_params",
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
        elseif arg == "--T-min"
            opts[:T_min_MeV] = parse(Float64, require_value())
        elseif arg == "--T-max"
            opts[:T_max_MeV] = parse(Float64, require_value())
        elseif arg == "--T-step"
            opts[:T_step_MeV] = parse(Float64, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--output-dir"
            opts[:output_dir] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = lowercase(require_value()) == "true"
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/relaxtime/export_friesen2019_mu0_order_params.jl [options]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return opts
end

const CLI_OPTS = parse_args(ARGS)
ENV["PNJL_PARAM_PROFILE"] = String(CLI_OPTS[:pnjl_profile])
ENV["PHYSICS_PARAM_PROFILE"] = String(CLI_OPTS[:physics_profile])
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

function write_csv(path::String, rows)
    open(path, "w") do io
        println(io, "T_MeV,phi_u,phi_d,phi_s,phi_u_norm,phi_s_norm,Phi,Phibar,converged,iterations,residual_norm")
        for row in rows
            println(io,
                "$(row.T_MeV),$(row.phi_u),$(row.phi_d),$(row.phi_s),$(row.phi_u_norm),$(row.phi_s_norm),$(row.Phi),$(row.Phibar),$(row.converged),$(row.iterations),$(row.residual_norm)"
            )
        end
    end
end

function main(args::Vector{String})
    opts = parse_args(args)
    outdir = abspath(String(opts[:output_dir]))
    isdir(outdir) && !Bool(opts[:overwrite]) && error("output_dir already exists: $outdir")
    mkpath(outdir)

    model = Main.Models.create_model(:PNJL)
    mode = Main.Models.FixedMu()
    hbarc = Float64(getfield(Main.Models.Constants_PNJL, Symbol("ħc_MeV_fm")))

    T_values = collect(Float64(opts[:T_min_MeV]):Float64(opts[:T_step_MeV]):Float64(opts[:T_max_MeV]))
    seed = Float64.(Main.Models.HADRON_SEED_5)
    rows = NamedTuple[]
    phi0_u = nothing
    phi0_s = nothing

    for T_MeV in T_values
        T_fm = T_MeV / hbarc
        result = Main.Models.solve_constraint(
            model,
            mode,
            T_fm;
            μ_fm=0.0,
            seed_guess=seed,
            xi=Float64(opts[:xi]),
            p_num=Int(opts[:p_num]),
            t_num=Int(opts[:t_num]),
        )

        xv = result.x_state
        phi_u = Float64(xv[1])
        phi_d = Float64(xv[2])
        phi_s = Float64(xv[3])
        Phi = Float64(xv[4])
        Phibar = Float64(xv[5])
        if Bool(result.converged)
            seed = [phi_u, phi_d, phi_s, Phi, Phibar]
        end
        if phi0_u === nothing && Bool(result.converged)
            phi0_u = phi_u
        end
        if phi0_s === nothing && Bool(result.converged)
            phi0_s = phi_s
        end
        push!(rows, (
            T_MeV=T_MeV,
            phi_u=phi_u,
            phi_d=phi_d,
            phi_s=phi_s,
            phi_u_norm=(phi0_u === nothing || iszero(phi0_u)) ? NaN : phi_u / phi0_u,
            phi_s_norm=(phi0_s === nothing || iszero(phi0_s)) ? NaN : phi_s / phi0_s,
            Phi=Phi,
            Phibar=Phibar,
            converged=Bool(result.converged),
            iterations=Int(getproperty(result, :iterations)),
            residual_norm=Float64(getproperty(result, :residual_norm)),
        ))
    end

    csv_path = joinpath(outdir, "mu0_order_params.csv")
    write_csv(csv_path, rows)

    open(joinpath(outdir, "README.md"), "w") do io
        conv = [row for row in rows if row.converged]
        println(io, "# Friesen 2019 muB=0 order-parameter export")
        println(io)
        println(io, "- pnjl_profile: $(opts[:pnjl_profile])")
        println(io, "- physics_profile: $(opts[:physics_profile])")
        println(io, "- xi: $(opts[:xi])")
        println(io, "- T window: $(opts[:T_min_MeV]) to $(opts[:T_max_MeV]) MeV, step $(opts[:T_step_MeV]) MeV")
        println(io, "- converged points: $(length(conv))/$(length(rows))")
        if !isempty(conv)
            println(io, "- reference T≈0 point: $(conv[1].T_MeV) MeV")
            println(io, "- reference phi_u: $(conv[1].phi_u)")
            println(io, "- reference phi_s: $(conv[1].phi_s)")
        end
    end

    println("friesen mu0 order-parameter export completed")
    println("output_dir=$(outdir)")
    println("csv_path=$(csv_path)")
    println("points=$(length(rows))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
