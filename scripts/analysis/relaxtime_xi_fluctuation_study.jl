using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, gauleg
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.OneLoopIntegrals: A

const TransportWorkflow = Models.transport_workflow_module()

struct Config
    label::String
    mode::Symbol
    p_nodes::Int
    angle_nodes::Int
    phi_nodes::Int
    n_sigma_points::Int
    sigma_grid_n::Int
end

struct SolveResult
    result
    backend::Symbol
end

struct Options
    T_mev::Float64
    muB_mev::Float64
    xi_values::Vector{Float64}
    top_n::Int
end

function print_usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime_xi_fluctuation_study.jl [options]\n")
    println("Options:")
    println("  --T <MeV>             temperature (required)")
    println("  --muB <MeV>           baryon chemical potential (default 0)")
    println("  --xi-list <csv>       xi list such as -0.32,-0.30,-0.28 (required)")
    println("  --top-n <int>         top differing channels per species (default 5)")
    println("  -h, --help            show help")
end

function parse_args(args::Vector{String})
    T_mev = nothing
    muB_mev = 0.0
    xi_values = Float64[]
    top_n = 5

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            value = args[i + 1]
            i += 1
            return value
        end

        if arg == "--T"
            T_mev = parse(Float64, require_value())
        elseif arg == "--muB"
            muB_mev = parse(Float64, require_value())
        elseif arg == "--xi-list"
            xi_values = Float64[parse(Float64, strip(v)) for v in split(require_value(), ',') if !isempty(strip(v))]
        elseif arg == "--top-n"
            top_n = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    T_mev === nothing && error("--T is required")
    isempty(xi_values) && error("--xi-list is required")
    return Options(T_mev, muB_mev, sort(unique(xi_values)), top_n)
end

function build_K_data(T_fm::Float64, mu_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    A_u = A(masses.u, mu_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, mu_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function integration_grids(cfg::Config)
    if cfg.mode == :finite_15
        p_grid, p_w = gauleg(0.0, 15.0, cfg.p_nodes)
        return (p_grid, p_w, Λ_inv_fm)
    elseif cfg.mode == :finite_lambda
        p_grid, p_w = gauleg(0.0, Λ_inv_fm, cfg.p_nodes)
        return (p_grid, p_w, Λ_inv_fm)
    end
    return (nothing, nothing, nothing)
end

function solve_equilibrium_with_fallback(T_fm::Float64, muq_fm::Float64, xi::Float64)
    try
        base = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:models,
            p_num=12,
            t_num=6,
            seed_state=TransportWorkflow.PNJL.HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )
        return (base, :models)
    catch err
        @warn "models equilibrium solver failed, fallback to legacy" xi=xi err=err
        base = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            muq_fm;
            xi=xi,
            solver_backend=:legacy,
            p_num=12,
            t_num=6,
            seed_state=nothing,
            solver_kwargs=(iterations=40,),
        )
        return (base, :legacy)
    end
end

function solve_tau_result(T_mev::Float64, muB_mev::Float64, xi::Float64, cfg::Config)
    T_fm = T_mev / ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / ħc_MeV_fm

    base, backend = solve_equilibrium_with_fallback(T_fm, muq_fm, xi)

    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))
    ktmp = build_K_data(T_fm, muq_fm, masses, Φ, Φbar)

    p_grid, p_w, sigma_cutoff = integration_grids(cfg)
    cos_grid, cos_w = gauleg(-1.0, 1.0, cfg.angle_nodes)
    phi_grid, phi_w = gauleg(0.0, 2 * pi, cfg.phi_nodes)

    result = TransportWorkflow.solve_gap_and_transport(
        T_fm,
        muq_fm;
        xi=xi,
        equilibrium=base,
        compute_tau=true,
        K_coeffs=ktmp.K_coeffs,
        tau=nothing,
        compute_bulk=false,
        p_num=12,
        t_num=6,
        seed_state=Vector(base.x_state),
        solver_kwargs=(iterations=40,),
        tau_kwargs=(
            p_nodes=cfg.p_nodes,
            angle_nodes=cfg.angle_nodes,
            phi_nodes=cfg.phi_nodes,
            n_sigma_points=cfg.n_sigma_points,
            sigma_grid_n=cfg.sigma_grid_n,
            p_grid=p_grid,
            p_w=p_w,
            cos_grid=cos_grid,
            cos_w=cos_w,
            phi_grid=phi_grid,
            phi_w=phi_w,
            sigma_cutoff=sigma_cutoff,
        ),
        transport_config=TransportWorkflow.TransportIntegrationConfig(p_nodes=24, p_max=8.0),
    )
    return SolveResult(result, backend)
end

function expand_isospin_densities(densities)
    return (
        u=densities.u,
        d=densities.u,
        s=densities.s,
        ubar=densities.ubar,
        dbar=densities.ubar,
        sbar=densities.sbar,
    )
end

function top_rate_deltas(base_rates, ref_rates; top_n::Int)
    rows = NamedTuple[]
    for key in propertynames(ref_rates)
        ref_val = getproperty(ref_rates, key)
        base_val = getproperty(base_rates, key)
        absdelta = abs(base_val - ref_val)
        relerr = absdelta / max(abs(ref_val), 1e-12)
        push!(rows, (channel=String(key), base=base_val, ref=ref_val, absdelta=absdelta, relerr=relerr))
    end
    sort!(rows, by=row -> row.absdelta, rev=true)
    return first(rows, min(top_n, length(rows)))
end

function top_contribution_deltas(base_res, ref_res, species::Symbol; top_n::Int)
    base_rows = Main.RelaxationTime.relaxation_rate_contribution_rows(expand_isospin_densities(base_res.densities), base_res.rates)
    ref_rows = Main.RelaxationTime.relaxation_rate_contribution_rows(expand_isospin_densities(ref_res.densities), ref_res.rates)
    ref_map = Dict((row.species, row.channel) => row for row in ref_rows)
    rows = NamedTuple[]
    for row in base_rows
        row.species == species || continue
        ref_row = ref_map[(row.species, row.channel)]
        delta = row.contribution - ref_row.contribution
        relerr = abs(delta) / max(abs(ref_row.contribution), 1e-12)
        push!(rows, (channel=String(row.channel), delta=delta, relerr=relerr, base=row.contribution, ref=ref_row.contribution))
    end
    sort!(rows, by=row -> abs(row.delta), rev=true)
    return first(rows, min(top_n, length(rows)))
end

function report_species(species::Symbol, base_res, ref_res, top_n::Int)
    tau_base = getproperty(base_res.tau, species)
    tau_ref = getproperty(ref_res.tau, species)
    rel = abs(tau_base - tau_ref) / max(abs(tau_ref), 1e-12)
    @printf("  %-5s tau_base=%10.6f tau_ref=%10.6f relerr=%8.4f\n", String(species), tau_base, tau_ref, rel)
    for row in top_contribution_deltas(base_res, ref_res, species; top_n=top_n)
        @printf("    contrib %-18s delta=% .6e relerr=%8.4f base=% .6e ref=% .6e\n",
            row.channel, row.delta, row.relerr, row.base, row.ref)
    end
end

function main()
    opts = parse_args(copy(ARGS))
    current_cfg = Config("current", :semi_infinite, 28, 8, 8, 6, 128)
    ref_cfg = Config("reference", :semi_infinite, 40, 8, 12, 24, 512)

    println("[xi-fluctuation-study] current vs reference")
    @printf("  T=%6.2f MeV  muB=%6.2f MeV\n", opts.T_mev, opts.muB_mev)
    println("  current   : semi_infinite p28 a8 phi8 sigma_n6 grid128")
    println("  reference : semi_infinite p40 a8 phi12 sigma_n24 grid512")

    previous_tau = Dict{Symbol,Float64}()
    for xi in opts.xi_values
        println()
        @printf("xi = % .4f\n", xi)
        try
            current = solve_tau_result(opts.T_mev, opts.muB_mev, xi, current_cfg)
            ref = solve_tau_result(opts.T_mev, opts.muB_mev, xi, ref_cfg)

            @printf("  equilibrium backend: current=%s reference=%s\n", String(current.backend), String(ref.backend))

            for species in (:u, :s, :ubar, :sbar)
                if haskey(previous_tau, species)
                    jump = getproperty(current.result.tau, species) - previous_tau[species]
                    @printf("  current Δtau %-5s = % .6f\n", String(species), jump)
                end
                previous_tau[species] = getproperty(current.result.tau, species)
            end

            println("  top rate deltas")
            for row in top_rate_deltas(current.result.rates, ref.result.rates; top_n=opts.top_n)
                @printf("    rate %-21s delta=% .6e relerr=%8.4f base=% .6e ref=% .6e\n",
                    row.channel, row.base - row.ref, row.relerr, row.base, row.ref)
            end

            report_species(:u, current.result, ref.result, opts.top_n)
            report_species(:s, current.result, ref.result, opts.top_n)
            report_species(:ubar, current.result, ref.result, opts.top_n)
            report_species(:sbar, current.result, ref.result, opts.top_n)
        catch err
            @printf("  ERROR at xi=% .4f: %s\n", xi, sprint(showerror, err))
        end
    end
end

main()