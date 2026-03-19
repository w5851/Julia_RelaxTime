using Printf
using Statistics

if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    include(joinpath(@__DIR__, "..", "..", "tests", "validation", "relaxtime", "literature_validation_helpers.jl"))
end

const HBARC = Main.Constants_PNJL.ħc_MeV_fm
const MB_PER_FM2 = 10.0
const LEGACY_SIGMA_TARGETS = joinpath(
    @__DIR__, "..", "..", "tests", "validation", "data", "targets", "relaxtime", "legacy", "sigma",
    "relaxtime_sigma_usbar_mu0_legacy_consensus_targets_v1.csv",
)

function load_legacy_sigma_targets(path::String)
    rows = NamedTuple[]
    for (idx, line) in enumerate(readlines(path))
        idx == 1 && continue
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        cols = split(s, ',')
        push!(rows, (
            record_id=strip(cols[1]),
            series=strip(cols[2]),
            process=Symbol(strip(cols[3])),
            T_MeV=parse(Float64, strip(cols[4])),
            muB_MeV=parse(Float64, strip(cols[5])),
            sqrt_s_GeV=parse(Float64, strip(cols[6])),
            fortran_sigma_mb=parse(Float64, strip(cols[7])),
            cpp_sigma_mb=parse(Float64, strip(cols[8])),
            consensus_sigma_mb=parse(Float64, strip(cols[9])),
        ))
    end
    return rows
end

function sigma_mb_at_params(params, process::Symbol, sqrt_s_GeV::Float64; n_points::Int)
    sqrt_s_fm = (sqrt_s_GeV * 1000.0) / HBARC
    s_fm2 = sqrt_s_fm^2
    return Main.TotalCrossSection.total_cross_section(
        process,
        s_fm2,
        params.quark_params,
        params.thermo_params,
        params.K_coeffs;
        n_points=n_points,
    ) * MB_PER_FM2
end

function threshold_gap_mev(params, sqrt_s_GeV::Float64)
    masses = params.quark_params.m
    threshold_mev = (masses.u + masses.s) * HBARC
    return sqrt_s_GeV * 1000.0 - threshold_mev
end

function summarize_error_block(name::String, values::Vector{Float64})
    @printf("  %-24s mean=%10.6f  median=%10.6f  max=%10.6f\n",
        name, mean(values), median(values), maximum(values))
end

function analyze_sigma_accuracy()
    rows = load_legacy_sigma_targets(LEGACY_SIGMA_TARGETS)
    old_errors = Float64[]
    old_errors_nonthreshold = Float64[]
    prod6_errors = Float64[]
    prod6_errors_nonthreshold = Float64[]
    val64_errors = Float64[]
    val64_errors_nonthreshold = Float64[]

    println("[sigma] legacy consensus vs Julia high-accuracy reference")
    println("  reference: current Julia equilibrium/K, total_cross_section n_points=256")
    println("  prod6    : current production-like sigma quadrature n_points=6")
    println("  val64    : validation helper sigma quadrature n_points=64")

    for row in rows
        params = _build_relaxtime_literature_sigma_params(row.T_MeV, row.muB_MeV)
        ref256 = sigma_mb_at_params(params, row.process, row.sqrt_s_GeV; n_points=256)
        prod6 = sigma_mb_at_params(params, row.process, row.sqrt_s_GeV; n_points=6)
        val64 = sigma_mb_at_params(params, row.process, row.sqrt_s_GeV; n_points=64)
        gap_mev = threshold_gap_mev(params, row.sqrt_s_GeV)
        is_threshold_sensitive = gap_mev < 5.0

        old_err = abs(row.consensus_sigma_mb - ref256) / max(abs(ref256), 1e-12)
        prod6_err = abs(prod6 - ref256) / max(abs(ref256), 1e-12)
        val64_err = abs(val64 - ref256) / max(abs(ref256), 1e-12)
        push!(old_errors, old_err)
        push!(prod6_errors, prod6_err)
        push!(val64_errors, val64_err)
        if !is_threshold_sensitive
            push!(old_errors_nonthreshold, old_err)
            push!(prod6_errors_nonthreshold, prod6_err)
            push!(val64_errors_nonthreshold, val64_err)
        end

        if row.sqrt_s_GeV < 0.55 || row.sqrt_s_GeV > 1.15
            @printf("    T=%5.1f sqrt(s)=%.6f GeV gap=%8.3f MeV | old=%9.5f ref=%9.5f prod6=%9.5f\n",
                row.T_MeV, row.sqrt_s_GeV, gap_mev, row.consensus_sigma_mb, ref256, prod6)
        end
    end

    println("  all target points:")
    summarize_error_block("legacy consensus", old_errors)
    summarize_error_block("julia prod n=6", prod6_errors)
    summarize_error_block("julia val n=64", val64_errors)

    println("  excluding near-threshold points (gap < 5 MeV):")
    summarize_error_block("legacy consensus", old_errors_nonthreshold)
    summarize_error_block("julia prod n=6", prod6_errors_nonthreshold)
    summarize_error_block("julia val n=64", val64_errors_nonthreshold)
end

function compute_tau_point(T_MeV::Float64, muB_MeV::Float64; mode::Symbol, p_nodes::Int, angle_nodes::Int,
    phi_nodes::Int, n_sigma_points::Int, sigma_grid_n::Int)
    return compute_tau_result(T_MeV, muB_MeV;
        mode=mode,
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        n_sigma_points=n_sigma_points,
        sigma_grid_n=sigma_grid_n,
    ).tau
end

function compute_tau_result(T_MeV::Float64, muB_MeV::Float64; mode::Symbol, p_nodes::Int, angle_nodes::Int,
    phi_nodes::Int, n_sigma_points::Int, sigma_grid_n::Int)
    tp = _compute_relaxtime_literature_transport_point(T_MeV, muB_MeV)
    equilibrium = tp.equilibrium
    T_fm = T_MeV / HBARC
    muq_fm = (muB_MeV / 3.0) / HBARC

    p_grid = nothing
    p_w = nothing
    sigma_cutoff = nothing
    if mode == :finite_15
        p_grid, p_w = Main.GaussLegendre.gauleg(0.0, 15.0, p_nodes)
        sigma_cutoff = Main.Constants_PNJL.Λ_inv_fm
    elseif mode == :finite_lambda
        p_grid, p_w = Main.GaussLegendre.gauleg(0.0, Main.Constants_PNJL.Λ_inv_fm, p_nodes)
        sigma_cutoff = Main.Constants_PNJL.Λ_inv_fm
    end

    cos_grid, cos_w = Main.GaussLegendre.gauleg(-1.0, 1.0, angle_nodes)
    phi_grid, phi_w = Main.GaussLegendre.gauleg(0.0, 2π, phi_nodes)

    return RELAXTIME_LITERATURE_VALIDATION_TRANSPORT_WORKFLOW.solve_gap_and_transport(
        T_fm,
        muq_fm;
        xi=0.0,
        equilibrium=equilibrium,
        compute_tau=true,
        K_coeffs=tp.K_coeffs,
        compute_bulk=false,
        p_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_P_NUM,
        t_num=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_T_NUM,
        seed_state=Vector(equilibrium.x_state),
        solver_kwargs=(iterations=RELAXTIME_LITERATURE_VALIDATION_VALIDATION_MAX_ITER,),
        tau_kwargs=(
            p_nodes=p_nodes,
            angle_nodes=angle_nodes,
            phi_nodes=phi_nodes,
            p_grid=p_grid,
            p_w=p_w,
            cos_grid=cos_grid,
            cos_w=cos_w,
            phi_grid=phi_grid,
            phi_w=phi_w,
            n_sigma_points=n_sigma_points,
            sigma_grid_n=sigma_grid_n,
            sigma_cutoff=sigma_cutoff,
        ),
    )
end

function report_tau_error(label::String, tau, ref)
    err_u = abs(tau.u - ref.u) / max(abs(ref.u), 1e-12)
    err_s = abs(tau.s - ref.s) / max(abs(ref.s), 1e-12)
    err_ubar = abs(tau.ubar - ref.ubar) / max(abs(ref.ubar), 1e-12)
    err_sbar = abs(tau.sbar - ref.sbar) / max(abs(ref.sbar), 1e-12)
    @printf("  %-16s err_u=%9.5f err_s=%9.5f err_ubar=%9.5f err_sbar=%9.5f\n",
        label, err_u, err_s, err_ubar, err_sbar)
end

function analyze_tau_accuracy()
    current_cfg = (
        mode=:semi_infinite,
        p_nodes=Main.AverageScatteringRate.DEFAULT_P_NODES,
        angle_nodes=Main.AverageScatteringRate.DEFAULT_ANGLE_NODES,
        phi_nodes=Main.AverageScatteringRate.DEFAULT_PHI_NODES,
        n_sigma_points=6,
        sigma_grid_n=128,
    )
    points = [
        (T=210.0, muB=0.0, tag="muB0_T210"),
        (T=250.0, muB=0.0, tag="muB0_T250"),
        (T=260.0, muB=0.0, tag="muB0_T260"),
        (T=150.0, muB=800.0, tag="muB800_T150"),
    ]

    println("[tau] same-equilibrium numerical convergence study")
    println("  ref        : semi_infinite, p40 angle8 phi12, sigma_n=24, sigma_grid=512")
    println("  current    : semi_infinite, p$(current_cfg.p_nodes) angle$(current_cfg.angle_nodes) phi$(current_cfg.phi_nodes), sigma_n=6,  sigma_grid=128")
    println("  legacy-like: finite_15,     p20 angle4 phi8, sigma_n=6,  sigma_grid=128")

    for point in points
        println("  point $(point.tag)")
        ref = compute_tau_point(point.T, point.muB;
            mode=:semi_infinite,
            p_nodes=40,
            angle_nodes=8,
            phi_nodes=12,
            n_sigma_points=24,
            sigma_grid_n=512,
        )
        current = compute_tau_point(point.T, point.muB; current_cfg...)
        legacy_like = compute_tau_point(point.T, point.muB;
            mode=:finite_15,
            p_nodes=20,
            angle_nodes=4,
            phi_nodes=8,
            n_sigma_points=6,
            sigma_grid_n=128,
        )
        report_tau_error("current", current, ref)
        report_tau_error("legacy-like", legacy_like, ref)
        @printf("    ref tau_u=%9.5f tau_s=%9.5f tau_ubar=%9.5f tau_sbar=%9.5f\n",
            ref.u, ref.s, ref.ubar, ref.sbar)
    end
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

function summarize_top_rate_deltas(base_rates, ref_rates; top_n::Int=8)
    rows = NamedTuple[]
    for key in propertynames(ref_rates)
        ref_val = getproperty(ref_rates, key)
        base_val = getproperty(base_rates, key)
        rel = abs(base_val - ref_val) / max(abs(ref_val), 1e-12)
        absdelta = abs(base_val - ref_val)
        push!(rows, (channel=String(key), base=base_val, ref=ref_val, relerr=rel, absdelta=absdelta))
    end
    sort!(rows, by=row -> row.relerr, rev=true)
    return first(rows, min(top_n, length(rows)))
end

function summarize_top_contribution_deltas(base_res, ref_res, species::Symbol; top_n::Int=8)
    base_rows = Main.RelaxationTime.relaxation_rate_contribution_rows(expand_isospin_densities(base_res.densities), base_res.rates)
    ref_rows = Main.RelaxationTime.relaxation_rate_contribution_rows(expand_isospin_densities(ref_res.densities), ref_res.rates)
    ref_map = Dict((row.species, row.channel) => row for row in ref_rows)
    rows = NamedTuple[]
    for row in base_rows
        row.species == species || continue
        ref_row = ref_map[(row.species, row.channel)]
        delta = row.contribution - ref_row.contribution
        rel = abs(delta) / max(abs(ref_row.contribution), 1e-12)
        push!(rows, (
            channel=row.channel,
            base=row.contribution,
            ref=ref_row.contribution,
            delta=delta,
            relerr=rel,
        ))
    end
    sort!(rows, by=row -> abs(row.delta), rev=true)
    return first(rows, min(top_n, length(rows)))
end

function analyze_tau_root_cause()
    T = 260.0
    muB = 0.0
    current_cfg = (
        label="current",
        mode=:semi_infinite,
        p_nodes=Main.AverageScatteringRate.DEFAULT_P_NODES,
        angle_nodes=Main.AverageScatteringRate.DEFAULT_ANGLE_NODES,
        phi_nodes=Main.AverageScatteringRate.DEFAULT_PHI_NODES,
        n_sigma_points=6,
        sigma_grid_n=128,
    )
    println("[tau-root-cause] point muB0_T260")
    println("  reference config: semi_infinite p40 a8 phi12 sigma_n24 grid512")

    ref_res = compute_tau_result(T, muB;
        mode=:semi_infinite,
        p_nodes=40,
        angle_nodes=8,
        phi_nodes=12,
        n_sigma_points=24,
        sigma_grid_n=512,
    )

    configs = [
        current_cfg,
        (label="grid512", mode=:semi_infinite, p_nodes=20, angle_nodes=4, phi_nodes=8, n_sigma_points=6, sigma_grid_n=512),
        (label="sigma24", mode=:semi_infinite, p_nodes=20, angle_nodes=4, phi_nodes=8, n_sigma_points=24, sigma_grid_n=128),
        (label="phase40812", mode=:semi_infinite, p_nodes=40, angle_nodes=8, phi_nodes=12, n_sigma_points=6, sigma_grid_n=128),
        (label="finite15", mode=:finite_15, p_nodes=20, angle_nodes=4, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
    ]

    for cfg in configs
        res = compute_tau_result(T, muB;
            mode=cfg.mode,
            p_nodes=cfg.p_nodes,
            angle_nodes=cfg.angle_nodes,
            phi_nodes=cfg.phi_nodes,
            n_sigma_points=cfg.n_sigma_points,
            sigma_grid_n=cfg.sigma_grid_n,
        )
        tau_s_err = abs(res.tau.s - ref_res.tau.s) / max(abs(ref_res.tau.s), 1e-12)
        tau_u_err = abs(res.tau.u - ref_res.tau.u) / max(abs(ref_res.tau.u), 1e-12)
        @printf("  %-10s tau_u_err=%9.5f tau_s_err=%9.5f tau_s=%9.5f\n", cfg.label, tau_u_err, tau_s_err, res.tau.s)
    end

    baseline = compute_tau_result(T, muB;
        mode=current_cfg.mode,
        p_nodes=current_cfg.p_nodes,
        angle_nodes=current_cfg.angle_nodes,
        phi_nodes=current_cfg.phi_nodes,
        n_sigma_points=current_cfg.n_sigma_points,
        sigma_grid_n=current_cfg.sigma_grid_n,
    )

    println("  top rate relative errors: baseline vs reference")
    for row in summarize_top_rate_deltas(baseline.rates, ref_res.rates)
        @printf("    %-18s relerr=%9.5f base=%10.6f ref=%10.6f\n", row.channel, row.relerr, row.base, row.ref)
    end

    println("  top s-channel contribution deltas: baseline vs reference")
    for row in summarize_top_contribution_deltas(baseline, ref_res, :s)
        @printf("    %-18s delta=%10.6f relerr=%9.5f base=%10.6f ref=%10.6f\n", String(row.channel), row.delta, row.relerr, row.base, row.ref)
    end

    println("  top sbar-channel contribution deltas: baseline vs reference")
    for row in summarize_top_contribution_deltas(baseline, ref_res, :sbar)
        @printf("    %-18s delta=%10.6f relerr=%9.5f base=%10.6f ref=%10.6f\n", String(row.channel), row.delta, row.relerr, row.base, row.ref)
    end

    println("  secondary sensitivity points")
    secondary_points = [
        (T=210.0, muB=0.0, tag="muB0_T210"),
        (T=150.0, muB=800.0, tag="muB800_T150"),
    ]
    configs = [
        current_cfg,
        (label="grid512", mode=:semi_infinite, p_nodes=20, angle_nodes=4, phi_nodes=8, n_sigma_points=6, sigma_grid_n=512),
        (label="sigma24", mode=:semi_infinite, p_nodes=20, angle_nodes=4, phi_nodes=8, n_sigma_points=24, sigma_grid_n=128),
        (label="phase40812", mode=:semi_infinite, p_nodes=40, angle_nodes=8, phi_nodes=12, n_sigma_points=6, sigma_grid_n=128),
        (label="finite15", mode=:finite_15, p_nodes=20, angle_nodes=4, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
    ]
    for point in secondary_points
        ref = compute_tau_result(point.T, point.muB;
            mode=:semi_infinite,
            p_nodes=40,
            angle_nodes=8,
            phi_nodes=12,
            n_sigma_points=24,
            sigma_grid_n=512,
        )
        println("    point $(point.tag)")
        for cfg in configs
            res = compute_tau_result(point.T, point.muB;
                mode=cfg.mode,
                p_nodes=cfg.p_nodes,
                angle_nodes=cfg.angle_nodes,
                phi_nodes=cfg.phi_nodes,
                n_sigma_points=cfg.n_sigma_points,
                sigma_grid_n=cfg.sigma_grid_n,
            )
            err_u = abs(res.tau.u - ref.tau.u) / max(abs(ref.tau.u), 1e-12)
            err_s = abs(res.tau.s - ref.tau.s) / max(abs(ref.tau.s), 1e-12)
            err_ubar = abs(res.tau.ubar - ref.tau.ubar) / max(abs(ref.tau.ubar), 1e-12)
            err_sbar = abs(res.tau.sbar - ref.tau.sbar) / max(abs(ref.tau.sbar), 1e-12)
            @printf("      %-10s err_u=%8.5f err_s=%8.5f err_ubar=%8.5f err_sbar=%8.5f\n",
                cfg.label, err_u, err_s, err_ubar, err_sbar)
        end
    end
end

function analyze_phase_space_candidates()
    current_cfg = (
        label="current",
        mode=:semi_infinite,
        p_nodes=Main.AverageScatteringRate.DEFAULT_P_NODES,
        angle_nodes=Main.AverageScatteringRate.DEFAULT_ANGLE_NODES,
        phi_nodes=Main.AverageScatteringRate.DEFAULT_PHI_NODES,
        n_sigma_points=6,
        sigma_grid_n=128,
    )
    points = [
        (T=210.0, muB=0.0, tag="muB0_T210"),
        (T=250.0, muB=0.0, tag="muB0_T250"),
        (T=260.0, muB=0.0, tag="muB0_T260"),
        (T=150.0, muB=800.0, tag="muB800_T150"),
    ]
    candidates = [
        current_cfg,
        (label="p24_a6_f8", mode=:semi_infinite, p_nodes=24, angle_nodes=6, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
        (label="p24_a6_f10", mode=:semi_infinite, p_nodes=24, angle_nodes=6, phi_nodes=10, n_sigma_points=6, sigma_grid_n=128),
        (label="p28_a6_f8", mode=:semi_infinite, p_nodes=28, angle_nodes=6, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
        (label="p28_a8_f8", mode=:semi_infinite, p_nodes=28, angle_nodes=8, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
        (label="p28_a6_f10", mode=:semi_infinite, p_nodes=28, angle_nodes=6, phi_nodes=10, n_sigma_points=6, sigma_grid_n=128),
        (label="p32_a6_f8", mode=:semi_infinite, p_nodes=32, angle_nodes=6, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
        (label="p32_a8_f8", mode=:semi_infinite, p_nodes=32, angle_nodes=8, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
        (label="p32_a8_f10", mode=:semi_infinite, p_nodes=32, angle_nodes=8, phi_nodes=10, n_sigma_points=6, sigma_grid_n=128),
        (label="p36_a8_f8", mode=:semi_infinite, p_nodes=36, angle_nodes=8, phi_nodes=8, n_sigma_points=6, sigma_grid_n=128),
    ]

    println("[tau-phase-space-candidates] balanced semi_infinite candidates")
    println("  reference: semi_infinite, p40 angle8 phi12, sigma_n=24, sigma_grid=512")

    for point in points
        ref = compute_tau_result(point.T, point.muB;
            mode=:semi_infinite,
            p_nodes=40,
            angle_nodes=8,
            phi_nodes=12,
            n_sigma_points=24,
            sigma_grid_n=512,
        )
        println("  point $(point.tag)")
        for cfg in candidates
            res = compute_tau_result(point.T, point.muB;
                mode=cfg.mode,
                p_nodes=cfg.p_nodes,
                angle_nodes=cfg.angle_nodes,
                phi_nodes=cfg.phi_nodes,
                n_sigma_points=cfg.n_sigma_points,
                sigma_grid_n=cfg.sigma_grid_n,
            )
            err_u = abs(res.tau.u - ref.tau.u) / max(abs(ref.tau.u), 1e-12)
            err_s = abs(res.tau.s - ref.tau.s) / max(abs(ref.tau.s), 1e-12)
            err_ubar = abs(res.tau.ubar - ref.tau.ubar) / max(abs(ref.tau.ubar), 1e-12)
            err_sbar = abs(res.tau.sbar - ref.tau.sbar) / max(abs(ref.tau.sbar), 1e-12)
            @printf("    %-10s err_u=%8.5f err_s=%8.5f err_ubar=%8.5f err_sbar=%8.5f\n",
                cfg.label, err_u, err_s, err_ubar, err_sbar)
        end
    end
end

function main()
    analyze_sigma_accuracy()
    println()
    analyze_tau_accuracy()
    println()
    analyze_tau_root_cause()
    println()
    analyze_phase_space_candidates()
end

main()