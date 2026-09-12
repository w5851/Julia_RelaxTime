"""Local physical-background ``muB`` sensitivity scan for charged GBU.

The scan keeps the baseline freeze-out temperature fixed at each selected
sqrt(s_NN) and perturbs the *upstream* baryon chemical potential.  Every
nonzero perturbation re-solves ``FixedMuBConservedCharges`` with
rho_Q/rho_B=0.4 and rho_S=0; it is not a free meson-chemical-potential scan.

This is an analysis-only probe.  It reuses the infinite-thermal GBU kernel and
phase/root bookkeeping, but uses a sparse q/omega mesh and does not establish a
new production baseline.  Failed or negative rows are retained explicitly.
"""
module ChargedGBUMuBSensitivity

using CSV
using Dates
using JSON3
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
isdefined(Main, :Models) || Base.include(Main, joinpath(ROOT, "src", "models", "Models.jl"))
isdefined(Main, :ChargedGBUResearchWorkflow) || Base.include(
    Main, joinpath(ROOT, "src", "models", "workflow_apps", "ChargedGBUResearchWorkflow.jl"))

const W = Main.ChargedGBUResearchWorkflow
const Models = Main.Models
const R = W.R
const I = W.I
const P = W.P
const Y = W.Y
const DEFAULT_CONFIG = joinpath(ROOT, "config", "models", "pnjl", "charged_gbu_infinite_v1.toml")
const DEFAULT_BASELINE = joinpath(ROOT, "data", "outputs", "results", "relaxtime", "meson_density",
    "charged_gbu_infinite", "freezeout_20260907_v2")
const DEFAULT_OUTPUT = joinpath(ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_phase_backend", "muB_sensitivity_20260911_v1")
const DEFAULT_ENERGIES = [3.0, 5.0, 7.7]
const DEFAULT_OFFSETS = [-50.0, 0.0, 50.0]
const DEFAULT_CHANNELS = (:pi_plus, :K_plus)

hashfile(path) = bytes2hex(sha256(read(path)))

function parse_float_list(value::AbstractString, name::AbstractString)
    tokens = split(strip(value), ',')
    isempty(tokens) && throw(ArgumentError("$(name) must not be empty"))
    values = try
        parse.(Float64, tokens)
    catch
        throw(ArgumentError("$(name) must be a comma-separated Float64 list"))
    end
    all(isfinite, values) && length(unique(values)) == length(values) ||
        throw(ArgumentError("$(name) must contain finite unique values"))
    return values
end

function validate_grids(energies, offsets)
    es = Float64.(collect(energies))
    ds = Float64.(collect(offsets))
    !isempty(es) && all(x -> isfinite(x) && x > 0.0, es) && length(unique(es)) == length(es) ||
        throw(ArgumentError("energies must be finite, positive and unique"))
    !isempty(ds) && all(isfinite, ds) && length(unique(ds)) == length(ds) ||
        throw(ArgumentError("muB offsets must be finite and unique"))
    0.0 in ds || throw(ArgumentError("muB offsets must include 0.0 for a central sensitivity estimate"))
    any(<(0.0), ds) && any(>(0.0), ds) ||
        throw(ArgumentError("muB offsets must include both negative and positive perturbations"))
    return (energies=sort(es), offsets=sort(ds))
end

"""Return deterministic points with the baseline (delta_muB=0) first per energy."""
function probe_points(energies, offsets; profile=Models.load_freezeout_profile(profile="default"))
    grids = validate_grids(energies, offsets)
    base = Models.build_freezeout_scan_points(grids.energies; profile=profile, traversal=:sqrts_ascending)
    rows = NamedTuple[]
    ordered_offsets = vcat([0.0], filter(!=(0.0), grids.offsets))
    for point in base
        for delta in ordered_offsets
            muB = point.muB_MeV + delta
            muB > 0.0 || throw(ArgumentError("muB perturbation leaves nonpositive muB at $(point.sqrt_s_NN_GeV) GeV"))
            push!(rows, (sqrt_s_NN_GeV=Float64(point.sqrt_s_NN_GeV),
                T_MeV=Float64(point.T_MeV), baseline_muB_MeV=Float64(point.muB_MeV),
                delta_muB_MeV=Float64(delta), muB_MeV=Float64(muB),
                on_baseline_freezeout=delta == 0.0,
                probe_type="fixed_T_muB_perturbation"))
        end
    end
    return rows
end

"""Sparse settings for a mechanism probe, not a replacement for production gates."""
Base.@kwdef struct SparseSettings
    mesh::Int = 96
    cut_nodes::Int = 48
    tail_nodes::Int = 48
    omega_nodes::Int = 64
    q_nodes::Int = 8
    qmax_inv_fm::Float64 = 8.0
    upper_inv_fm::Float64 = 48.0
end

function validate_settings(s::SparseSettings)
    all(x -> x >= 8, (s.mesh, s.cut_nodes, s.tail_nodes, s.omega_nodes, s.q_nodes)) ||
        throw(ArgumentError("sparse node counts must be >= 8"))
    all(isfinite, (s.qmax_inv_fm, s.upper_inv_fm)) && 0.0 < s.qmax_inv_fm < s.upper_inv_fm ||
        throw(ArgumentError("sparse q/omega limits must be finite and ordered"))
    return s
end

function sparse_channel_density(bg, channel::Symbol, settings::SparseSettings)
    validate_settings(settings)
    qs, weights = R.gauleg(0.0, settings.qmax_inv_fm, settings.q_nodes)
    rows = NamedTuple[]
    for (q, weight) in zip(qs, weights)
        kernel = I.kernel(bg, channel, q;
            cut_nodes=settings.cut_nodes,
            split_inv_fm=max(36.0, q + 28bg.T + 2.0))
        profile = P.profile(kernel; mesh=settings.mesh, tail_nodes=settings.tail_nodes)
        shell = Y.shell(profile; nodes=settings.omega_nodes, upper=settings.upper_inv_fm)
        push!(rows, (channel=String(channel), q_inv_fm=Float64(q), quadrature_weight=Float64(weight),
            density_inv_fm3=Float64(shell.density), bound_inv_fm3=Float64(shell.bound),
            landau_inv_fm3=Float64(shell.landau), pair_inv_fm3=Float64(shell.pair),
            root_count=Int(shell.root_count), negative_root_count=Int(shell.negative_root_count),
            omega_tail_bound=Float64(shell.omega_tail_conditional_bound),
            shell_passed=isfinite(shell.density), production_authorized=false))
    end
    density = sum(r.quadrature_weight * r.density_inv_fm3 for r in rows)
    bound = sum(r.quadrature_weight * r.bound_inv_fm3 for r in rows)
    landau = sum(r.quadrature_weight * r.landau_inv_fm3 for r in rows)
    pair = sum(r.quadrature_weight * r.pair_inv_fm3 for r in rows)
    return (channel=String(channel), density_inv_fm3=density, bound_inv_fm3=bound,
        landau_inv_fm3=landau, pair_inv_fm3=pair, passed=isfinite(density),
        negative_density=density < 0.0, max_root_count=maximum(r.root_count for r in rows),
        negative_root_nodes=count(r -> r.negative_root_count > 0, rows),
        rows=rows, status="sparse_diagnostic", production_authorized=false)
end

function ratio_from_channels(channel_results::AbstractDict)
    haskey(channel_results, "pi_plus") && haskey(channel_results, "K_plus") ||
        return (ratio=NaN, passed=false, reason="missing_plus_channels")
    pion = channel_results["pi_plus"]
    kaon = channel_results["K_plus"]
    ok = pion.passed && kaon.passed && pion.density_inv_fm3 > 0.0 && kaon.density_inv_fm3 >= 0.0
    return (ratio=ok ? kaon.density_inv_fm3 / pion.density_inv_fm3 : NaN,
        passed=ok, reason=ok ? "" : "nonpositive_or_failed_density")
end

function sensitivity_rows(ratio_rows; delta_muB_MeV::Real)
    delta = Float64(delta_muB_MeV)
    delta > 0.0 || throw(ArgumentError("delta_muB_MeV must be positive"))
    groups = Dict{Float64,Vector{Any}}()
    for row in ratio_rows
        push!(get!(groups, Float64(row.sqrt_s_NN_GeV), Any[]), row)
    end
    out = NamedTuple[]
    for energy in sort(collect(keys(groups)))
        group = groups[energy]
        at(d) = findfirst(r -> isapprox(r.delta_muB_MeV, d; atol=0.0, rtol=0.0), group)
        im, i0, ip = at(-delta), at(0.0), at(delta)
        if im === nothing || i0 === nothing || ip === nothing
            push!(out, (sqrt_s_NN_GeV=energy, delta_muB_MeV=delta,
                ratio_minus=NaN, ratio_0=NaN, ratio_plus=NaN,
                d_ratio_dmuB=NaN, d_logratio_dmuB=NaN, curvature=NaN,
                passed=false, reason="missing_central_triplet"))
            continue
        end
        rm, r0, rp = group[im].ratio, group[i0].ratio, group[ip].ratio
        ok = all(isfinite, (rm, r0, rp)) && rm > 0.0 && r0 > 0.0 && rp > 0.0
        push!(out, (sqrt_s_NN_GeV=energy, delta_muB_MeV=delta,
            ratio_minus=rm, ratio_0=r0, ratio_plus=rp,
            d_ratio_dmuB=ok ? (rp - rm) / (2delta) : NaN,
            d_logratio_dmuB=ok ? (log(rp) - log(rm)) / (2delta) : NaN,
            curvature=ok ? (rp - 2r0 + rm) / delta^2 : NaN,
            passed=ok, reason=ok ? "" : "nonpositive_or_failed_triplet"))
    end
    return out
end

function background_row(point, bg; source::AbstractString, seed_reused::Bool)
    h = Main.Constants_PNJL.ħc_MeV_fm
    mu_pi = (bg.mu.u - bg.mu.d) * h
    mu_k = (bg.mu.u - bg.mu.s) * h
    return (sqrt_s_NN_GeV=point.sqrt_s_NN_GeV, T_MeV=point.T_MeV,
        baseline_muB_MeV=point.baseline_muB_MeV, delta_muB_MeV=point.delta_muB_MeV,
        muB_MeV=point.muB_MeV, on_baseline_freezeout=point.on_baseline_freezeout,
        source=source, seed_reused=seed_reused, residual=bg.residual,
        m_u_inv_fm=bg.m.u, m_d_inv_fm=bg.m.d, m_s_inv_fm=bg.m.s,
        mu_u_inv_fm=bg.mu.u, mu_d_inv_fm=bg.mu.d, mu_s_inv_fm=bg.mu.s,
        mu_pi_plus_MeV=mu_pi, mu_K_plus_MeV=mu_k, Phi=bg.Phi, PhiBar=bg.PhiBar,
        production_authorized=false)
end

function read_baseline_point(baseline_output::AbstractString, point)
    dir = joinpath(baseline_output, "energy_$(point.sqrt_s_NN_GeV)")
    bg_path = joinpath(dir, "background.json")
    isfile(bg_path) || return nothing
    data = JSON3.read(read(bg_path, String))
    return (bg=W.restore_background(data), seed=Float64.(data.seed),
        source="freezeout_20260907_v2", path=bg_path)
end

function read_baseline_ratio(baseline_output::AbstractString, point)
    dir = joinpath(baseline_output, "energy_$(point.sqrt_s_NN_GeV)")
    kp = joinpath(dir, "K_plus.json")
    pp = joinpath(dir, "pi_plus.json")
    isfile(kp) && isfile(pp) || return NaN
    k = JSON3.read(read(kp, String)); p = JSON3.read(read(pp, String))
    k.passed && p.passed && isfinite(Float64(k.density)) && isfinite(Float64(p.density)) &&
        Float64(k.density) >= 0.0 && Float64(p.density) > 0.0 || return NaN
    return Float64(k.density) / Float64(p.density)
end

function source_hashes(config, baseline_output, energies)
    paths = [
        @__FILE__,
        joinpath(ROOT, "src", "models", "workflow_apps", "ChargedGBUResearchWorkflow.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_qgate.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_acceptance.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_yield.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_profile.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_thermal.jl"),
        joinpath(ROOT, "config", "physics", "freezeout", "default.toml"),
        abspath(config),
    ]
    baseline_manifest = joinpath(baseline_output, "manifest.json")
    isfile(baseline_manifest) && push!(paths, baseline_manifest)
    for energy in energies
        dir = joinpath(baseline_output, "energy_$(energy)")
        for name in ("background.json", "pi_plus.json", "K_plus.json")
            path = joinpath(dir, name)
            isfile(path) && push!(paths, path)
        end
    end
    return Dict(replace(relpath(path, ROOT), '\\' => '/') => hashfile(path) for path in unique(paths))
end

function write_json(path, value)
    open(path, "w") do io
        JSON3.write(io, value)
    end
end

function render_plot(ratio_rows, output)
    try
        @eval using Plots
    catch err
        return (status="plot_skipped", reason=sprint(showerror, err))
    end
    return Base.invokelatest(_render_plot_impl, ratio_rows, output)
end

function _render_plot_impl(ratio_rows, output)
    plt = plot(;xlabel="muB [MeV]", ylabel="K+/pi+ partial-yield ratio",
        title="Charged GBU muB sensitivity at fixed baseline T",
        size=(1100, 680), dpi=220, legend=:best, grid=true)
    energies = sort(unique(r.sqrt_s_NN_GeV for r in ratio_rows))
    for energy in energies
        group = sort(filter(r -> r.sqrt_s_NN_GeV == energy && r.ratio_passed, ratio_rows), by=r -> r.muB_MeV)
        isempty(group) && continue
        plot!(plt, [r.muB_MeV for r in group], [r.ratio for r in group];
            marker=:circle, linewidth=2, label="sqrt(s)=$(energy) GeV")
    end
    png = joinpath(output, "muB_sensitivity.png")
    pdf = joinpath(output, "muB_sensitivity.pdf")
    savefig(plt, png)
    savefig(plt, pdf)
    return (status="plot_written", png=png, pdf=pdf)
end

function replace_setting(settings::SparseSettings, field::Symbol, value)
    return SparseSettings(
        mesh=field === :mesh ? Int(value) : settings.mesh,
        cut_nodes=field === :cut_nodes ? Int(value) : settings.cut_nodes,
        tail_nodes=settings.tail_nodes,
        omega_nodes=field === :omega_nodes ? Int(value) : settings.omega_nodes,
        q_nodes=field === :q_nodes ? Int(value) : settings.q_nodes,
        qmax_inv_fm=field === :qmax_inv_fm ? Float64(value) : settings.qmax_inv_fm,
        upper_inv_fm=field === :upper_inv_fm ? Float64(value) : settings.upper_inv_fm,
    )
end

function run_scan(;output=DEFAULT_OUTPUT, baseline_output=DEFAULT_BASELINE, config=DEFAULT_CONFIG,
    energies=DEFAULT_ENERGIES, offsets=DEFAULT_OFFSETS, settings=SparseSettings(), make_plot=true)
    validate_settings(settings)
    grids = validate_grids(energies, offsets)
    output = abspath(output)
    ispath(output) && error("refusing to overwrite existing output: $(output)")
    mkpath(output)
    c = TOML.parsefile(config)
    profile = Models.load_freezeout_profile(profile=String(c["freezeout_profile"]))
    points = probe_points(grids.energies, grids.offsets; profile=profile)
    hashes = source_hashes(config, baseline_output, grids.energies)
    manifest = Dict{String,Any}(
        "schema" => "charged_gbu_mub_sensitivity_v1",
        "status" => "running",
        "background" => "FixedMuBConservedCharges quark-only BQS; rhoQ/rhoB=0.4; rhoS=0",
        "scan_axis" => "muB at fixed baseline freezeout T",
        "probe_type" => "physical upstream background perturbation; BQS re-solved",
        "energies_GeV" => grids.energies,
        "muB_offsets_MeV" => grids.offsets,
        "settings" => settings,
        "baseline_output" => replace(relpath(baseline_output, ROOT), '\\' => '/'),
        "baseline_route" => "charged_gbu_infinite_v1 accepted smoke result",
        "source_hashes" => hashes,
        "production_authorized" => false,
        "diagnostic_only" => true,
        "full_q_convergence" => false,
        "meson_feedback" => false,
    )
    write_json(joinpath(output, "progress.json"), manifest)

    model = Models.create_model(:PNJL)
    backgrounds = NamedTuple[]
    channel_totals = NamedTuple[]
    shells = NamedTuple[]
    ratios = NamedTuple[]
    failures = NamedTuple[]
    baseline_cache = Dict{Float64,Any}()
    for point in points
        base = get!(baseline_cache, point.sqrt_s_NN_GeV) do
            read_baseline_point(baseline_output, point)
        end
        solved = nothing
        source = "new_bqs_solve"
        seed_reused = false
        if point.delta_muB_MeV == 0.0 && base !== nothing
            solved = base
            source = base.source
            seed_reused = true
        else
            seed = base === nothing ? nothing : base.seed
            try
                solved = W.background(model,
                    (T_MeV=point.T_MeV, muB_MeV=point.muB_MeV), c; seed=seed)
                solved = (bg=solved.bg, seed=solved.seed, source=source, path=nothing)
            catch err
                push!(failures, (sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
                    delta_muB_MeV=point.delta_muB_MeV, stage="background",
                    reason=sprint(showerror, err), production_authorized=false))
                continue
            end
        end
        push!(backgrounds, background_row(point, solved.bg; source=source, seed_reused=seed_reused))

        result_map = Dict{String,Any}()
        for channel in DEFAULT_CHANNELS
            result = try
                sparse_channel_density(solved.bg, channel, settings)
            catch err
                push!(failures, (sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
                    delta_muB_MeV=point.delta_muB_MeV, stage=String(channel),
                    reason=sprint(showerror, err), production_authorized=false))
                nothing
            end
            result === nothing && continue
            result_map[String(channel)] = result
            push!(channel_totals, (sqrt_s_NN_GeV=point.sqrt_s_NN_GeV, T_MeV=point.T_MeV,
                baseline_muB_MeV=point.baseline_muB_MeV, delta_muB_MeV=point.delta_muB_MeV,
                muB_MeV=point.muB_MeV, channel=result.channel,
                density_inv_fm3=result.density_inv_fm3, bound_inv_fm3=result.bound_inv_fm3,
                landau_inv_fm3=result.landau_inv_fm3, pair_inv_fm3=result.pair_inv_fm3,
                passed=result.passed, negative_density=result.negative_density,
                max_root_count=result.max_root_count, negative_root_nodes=result.negative_root_nodes,
                q_nodes=settings.q_nodes, qmax_inv_fm=settings.qmax_inv_fm,
                omega_nodes=settings.omega_nodes, status=result.status, production_authorized=false))
            for shell_row in result.rows
                push!(shells, (sqrt_s_NN_GeV=point.sqrt_s_NN_GeV, T_MeV=point.T_MeV,
                    baseline_muB_MeV=point.baseline_muB_MeV, delta_muB_MeV=point.delta_muB_MeV,
                    muB_MeV=point.muB_MeV, channel=shell_row.channel, density_inv_fm3=shell_row.density_inv_fm3,
                    bound_inv_fm3=shell_row.bound_inv_fm3, landau_inv_fm3=shell_row.landau_inv_fm3,
                    pair_inv_fm3=shell_row.pair_inv_fm3, passed=shell_row.shell_passed,
                    negative_density=shell_row.density_inv_fm3 < 0.0, root_count=shell_row.root_count,
                    negative_root_count=shell_row.negative_root_count, q_nodes=settings.q_nodes,
                    qmax_inv_fm=settings.qmax_inv_fm, omega_nodes=settings.omega_nodes,
                    status="q_shell", production_authorized=false))
            end
        end
        ratio = ratio_from_channels(result_map)
        baseline_ratio = point.delta_muB_MeV == 0.0 ? read_baseline_ratio(baseline_output, point) : NaN
        push!(ratios, (sqrt_s_NN_GeV=point.sqrt_s_NN_GeV, T_MeV=point.T_MeV,
            baseline_muB_MeV=point.baseline_muB_MeV, delta_muB_MeV=point.delta_muB_MeV,
            muB_MeV=point.muB_MeV, ratio=ratio.ratio, ratio_passed=ratio.passed,
            reason=ratio.reason, baseline_production_ratio=baseline_ratio,
            source=source, sparse_route=true, production_authorized=false))
        CSV.write(joinpath(output, "backgrounds.csv"), backgrounds)
        CSV.write(joinpath(output, "channel_totals.csv"), channel_totals)
        CSV.write(joinpath(output, "shells.csv"), shells)
        CSV.write(joinpath(output, "ratios.csv"), ratios)
        CSV.write(joinpath(output, "failures.csv"), failures)
        println("[mub-sensitivity] sqrt_s=$(point.sqrt_s_NN_GeV) delta=$(point.delta_muB_MeV) ratio=$(ratio.ratio) pass=$(ratio.passed)")
        flush(stdout)
    end

    positive_offsets = filter(>(0.0), grids.offsets)
    derivative_delta = minimum(positive_offsets)
    derivatives = sensitivity_rows(ratios; delta_muB_MeV=derivative_delta)
    CSV.write(joinpath(output, "derivatives.csv"), derivatives)
    plot_record = make_plot ? render_plot(ratios, output) : (status="plot_disabled",)
    manifest["status"] = isempty(failures) ? "diagnostic_completed" : "diagnostic_completed_with_failures"
    manifest["point_count"] = length(points)
    manifest["background_count"] = length(backgrounds)
    manifest["channel_total_count"] = length(channel_totals)
    manifest["shell_count"] = length(shells)
    manifest["ratio_count"] = length(ratios)
    manifest["failed_count"] = length(failures)
    if plot_record.status == "plot_written"
        manifest["plot"] = merge(plot_record, (
            png=replace(relpath(plot_record.png, ROOT), '\\' => '/'),
            pdf=replace(relpath(plot_record.pdf, ROOT), '\\' => '/'),
            png_sha256=hashfile(plot_record.png), pdf_sha256=hashfile(plot_record.pdf)))
    else
        manifest["plot"] = plot_record
    end
    manifest["derivative_count"] = length(derivatives)
    manifest["production_default"] = false
    manifest["source_git_head"] = readchomp(`git -C $ROOT rev-parse HEAD`)
    manifest["generated_at_utc"] = string(Dates.now(Dates.UTC))
    write_json(joinpath(output, "manifest.json"), manifest)
    open(joinpath(output, "manifest.json.sha256"), "w") do io
        write(io, hashfile(joinpath(output, "manifest.json")))
    end
    return (output=output, status=manifest["status"], points=length(points), ratios=length(ratios), failures=length(failures))
end

function parse_args(args)
    options = Dict{Symbol,Any}(
        :output => DEFAULT_OUTPUT,
        :baseline_output => DEFAULT_BASELINE,
        :config => DEFAULT_CONFIG,
        :energies => DEFAULT_ENERGIES,
        :offsets => DEFAULT_OFFSETS,
        :settings => SparseSettings(),
        :make_plot => true,
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--help", "-h")
            return nothing
        elseif arg == "--no-plot"
            options[:make_plot] = false
        elseif arg in ("--output", "--baseline-output", "--config", "--sqrts-list", "--mub-offsets",
                       "--q-nodes", "--qmax", "--mesh", "--cut-nodes", "--omega-nodes", "--upper")
            i < length(args) || throw(ArgumentError("missing value for $(arg)"))
            i += 1
            value = args[i]
            if arg == "--output"
                options[:output] = value
            elseif arg == "--baseline-output"
                options[:baseline_output] = value
            elseif arg == "--config"
                options[:config] = value
            elseif arg == "--sqrts-list"
                options[:energies] = parse_float_list(value, "sqrts-list")
            elseif arg == "--mub-offsets"
                options[:offsets] = parse_float_list(value, "mub-offsets")
            else
                old = options[:settings]
                field = arg == "--q-nodes" ? :q_nodes : arg == "--qmax" ? :qmax_inv_fm :
                    arg == "--mesh" ? :mesh : arg == "--cut-nodes" ? :cut_nodes :
                    arg == "--omega-nodes" ? :omega_nodes : :upper_inv_fm
                parsed = field in (:q_nodes, :mesh, :cut_nodes, :omega_nodes) ? parse(Int, value) : parse(Float64, value)
                options[:settings] = replace_setting(old, field, parsed)
            end
        else
            throw(ArgumentError("unknown option: $(arg)"))
        end
        i += 1
    end
    return options
end

function main(args=ARGS)
    options = parse_args(args)
    options === nothing && return println("Usage: julia --project=. scripts/analysis/relaxtime/scan_charged_gbu_mub_sensitivity.jl [--output DIR] [--baseline-output DIR] [--config TOML] [--sqrts-list 3,5,7.7] [--mub-offsets -50,0,50] [--q-nodes 8] [--qmax 8] [--mesh 96] [--cut-nodes 48] [--omega-nodes 64] [--upper 48] [--no-plot]")
    result = run_scan(;options...)
    println(result)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end # module ChargedGBUMuBSensitivity
