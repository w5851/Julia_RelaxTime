#!/usr/bin/env julia

"""
Benchmark the runtime overhead introduced by adaptive sigma-cache refinement.

The public relaxation_times API does not currently expose the adaptive toggle,
so this benchmark compares the actual changed stage directly: sigma-cache
construction for all required scattering processes at representative thermal
points.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using Statistics
using .Constants_PNJL: G_fm2, K_fm5, Λ_inv_fm, ħc_MeV_fm
using .GaussLegendre: gauleg
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings

const PNJL = Models.pnjl_module()
const Integrals = getproperty(PNJL, :Integrals)
const DEFAULT_MOMENTUM_NODES = getproperty(Integrals, :DEFAULT_MOMENTUM_NODES)
const DEFAULT_MOMENTUM_WEIGHTS = getproperty(Integrals, :DEFAULT_MOMENTUM_WEIGHTS)
using .OneLoopIntegrals: A

const TransportWorkflow = Models.transport_workflow_module()
const RT_ASR = Main.AverageScatteringRate
const RT_TCS = Main.TotalCrossSection
const RelaxationTime = Main.RelaxationTime
const REQUIRED_PROCESSES = RelaxationTime.REQUIRED_PROCESSES

struct Options
    out::String
    temperatures_mev::Vector{Float64}
    muB_mev::Float64
    xi::Float64
    repeats::Int
    integration_mode::Symbol
    sigma_grid_n::Int
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/benchmark_adaptive_refinement_overhead.jl [options]\n")
    println("Options:")
    println("  --out <path>           output CSV (default data/outputs/results/relaxtime/benchmark/adaptive_refinement_overhead.csv)")
    println("  --temperatures <csv>   benchmark temperatures in MeV (default 150,190,250)")
    println("  --muB <MeV>            baryon chemical potential (default 0)")
    println("  --xi <value>           anisotropy parameter (default 0)")
    println("  --repeats <int>        timed repeats per variant (default 3)")
    println("  --mode <mode>          semi_infinite|finite_15|finite_lambda (default semi_infinite)")
    println("  --sigma-grid-n <int>   sigma cache grid size (default $(RT_ASR.DEFAULT_SIGMA_GRID_N))")
    println("  -h, --help             show help")
end

function parse_float_list(raw::AbstractString)::Vector{Float64}
    vals = Float64[]
    for part in split(raw, ',')
        s = strip(part)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && error("empty float list")
    return vals
end

function parse_args(args::Vector{String})::Options
    opts = Dict{Symbol,Any}(
        :out => joinpath("data", "outputs", "results", "relaxtime", "benchmark", "adaptive_refinement_overhead.csv"),
        :temperatures_mev => Float64[150.0, 190.0, 250.0],
        :muB_mev => 0.0,
        :xi => 0.0,
        :repeats => 3,
        :integration_mode => :semi_infinite,
        :sigma_grid_n => RT_ASR.DEFAULT_SIGMA_GRID_N,
        :tau_p_nodes => RT_ASR.DEFAULT_P_NODES,
        :tau_angle_nodes => RT_ASR.DEFAULT_ANGLE_NODES,
        :tau_phi_nodes => RT_ASR.DEFAULT_PHI_NODES,
        :tau_n_sigma_points => RT_TCS.DEFAULT_T_INTEGRAL_POINTS,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            value = args[i + 1]
            i += 1
            return value
        end

        if arg == "--out"
            opts[:out] = require_value()
        elseif arg == "--temperatures"
            opts[:temperatures_mev] = parse_float_list(require_value())
        elseif arg == "--muB"
            opts[:muB_mev] = parse(Float64, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--repeats"
            opts[:repeats] = parse(Int, require_value())
        elseif arg == "--mode"
            mode = Symbol(require_value())
            mode in (:semi_infinite, :finite_15, :finite_lambda) || error("unknown mode: $mode")
            opts[:integration_mode] = mode
        elseif arg == "--sigma-grid-n"
            opts[:sigma_grid_n] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return Options(
        String(opts[:out]),
        Float64.(opts[:temperatures_mev]),
        Float64(opts[:muB_mev]),
        Float64(opts[:xi]),
        Int(opts[:repeats]),
        Symbol(opts[:integration_mode]),
        Int(opts[:sigma_grid_n]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma_points]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function build_K_data(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Phi::Float64, Phibar::Float64)
    A_u = A(masses.u, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (
        K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s),
        A_vals=(u=A_u, d=A_u, s=A_s),
    )
end

function solve_base_state(T_mev::Float64, muB_mev::Float64, xi::Float64)
    T_fm = T_mev / ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / ħc_MeV_fm
    base = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
        T_fm,
        muq_fm;
        xi=xi,
        solver_backend=:models,
        p_num=12,
        t_num=6,
        models_residual_norm_max=1e-4,
    )
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))
    Phi = Float64(base.x_state[4])
    Phibar = Float64(base.x_state[5])
    kdata = build_K_data(T_fm, muq_fm, masses, Phi, Phibar)
    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=kdata.A_vals)
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    return (
        T_fm=T_fm,
        muq_fm=muq_fm,
        quark_params=quark_params,
        thermo_params=thermo_params,
        K_coeffs=kdata.K_coeffs,
    )
end

function benchmark_variant(state, opts::Options, adaptive_refinement::Bool)
    function run_once()
        caches = Dict{Symbol,RT_ASR.CrossSectionCache}()
        for process in REQUIRED_PROCESSES
            caches[process] = RT_ASR.build_w0cdf_pchip_cache(
                process,
                state.quark_params,
                state.thermo_params,
                state.K_coeffs;
                N=opts.sigma_grid_n,
                p_cutoff=Λ_inv_fm,
                n_sigma_points=opts.tau_n_sigma_points,
                adaptive_refinement=adaptive_refinement,
            )
        end
        return caches
    end

    run_once()
    GC.gc()

    times = Float64[]
    for _ in 1:opts.repeats
        GC.gc()
        elapsed = @elapsed begin
            caches = run_once()
            length(caches) == length(REQUIRED_PROCESSES) || error("incomplete cache set during benchmark")
        end
        push!(times, elapsed)
    end
    return times
end

function main()
    opts = parse_args(copy(ARGS))
    ensure_parent_dir(opts.out)

    rows = NamedTuple[]
    open(opts.out, "w") do io
        println(io, "T_MeV,muB_MeV,xi,variant,repeat,elapsed_s")
        for T_mev in opts.temperatures_mev
            state = solve_base_state(T_mev, opts.muB_mev, opts.xi)
            baseline_times = benchmark_variant(state, opts, false)
            refined_times = benchmark_variant(state, opts, true)

            for (index, value) in enumerate(baseline_times)
                println(io, join((T_mev, opts.muB_mev, opts.xi, "adaptive_off", index, value), ','))
                push!(rows, (T_MeV=T_mev, variant="adaptive_off", elapsed_s=value))
            end
            for (index, value) in enumerate(refined_times)
                println(io, join((T_mev, opts.muB_mev, opts.xi, "adaptive_on", index, value), ','))
                push!(rows, (T_MeV=T_mev, variant="adaptive_on", elapsed_s=value))
            end

            baseline_mean = mean(baseline_times)
            refined_mean = mean(refined_times)
            slowdown = (refined_mean / baseline_mean - 1.0) * 100.0
            println("T=$(T_mev) MeV | adaptive_off=$(round(baseline_mean; digits=3)) s | adaptive_on=$(round(refined_mean; digits=3)) s | slowdown=$(round(slowdown; digits=1))%")
        end
    end

    baseline_all = [row.elapsed_s for row in rows if row.variant == "adaptive_off"]
    refined_all = [row.elapsed_s for row in rows if row.variant == "adaptive_on"]
    total_slowdown = (mean(refined_all) / mean(baseline_all) - 1.0) * 100.0
    println("Overall mean slowdown: $(round(total_slowdown; digits=1))%")
    println("Wrote benchmark CSV: $(opts.out)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end