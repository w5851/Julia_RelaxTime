"""Compare direct finite-q and q=0 lambda-reference GBU freezeout diagnostics.

This is an analysis-only report.  It reconstructs component densities from the
retained q-shell rows, writes a current-route comparison, and summarizes the
historical FixedAsymmetricRho artifact without treating it as a matched
baseline.  No production default or formal data product is changed.
"""
module CausalGBUFreezeoutComparison

include("causal_gbu_research_utils.jl")
using .CausalGBUResearch
using CSV, JSON3, Statistics
using Plots
const R = CausalGBUResearch

const COMPONENTS = (:bound, :unitary, :landau)
const ROUTE_LABELS = Dict(
    "direct_finite_q" => "direct finite-q",
    "q0_lambda_reference" => "q=0 lambda reference",
)

function _rows(path)
    isfile(path) || error("missing CSV input: $(path)")
    return collect(CSV.File(path; comment="#"))
end

"""Integrate component shell rows with the same q quadrature as the producer."""
function aggregate_components(rows; qmax=8.0, nq=24)
    qnodes, weights = R.gauleg(0.0, qmax, nq)
    grouped = Dict{Tuple{Float64,String,String},Vector{Any}}()
    for row in rows
        key = (Float64(row.sqrt_s_NN_GeV), String(row.channel), String(row.route))
        push!(get!(grouped, key, Any[]), row)
    end
    out = NamedTuple[]
    for (key, group) in sort!(collect(grouped); by=first)
        length(group) == nq || error("q-shell count mismatch for $(key): $(length(group)) != $(nq)")
        sort!(group; by=r -> Float64(r.q_inv_fm))
        for (i, q) in enumerate(qnodes)
            abs(Float64(group[i].q_inv_fm) - q) < 1e-10 ||
                error("q-node mismatch for $(key) at index $(i)")
        end
        # The retained component columns are pre-measure shell terms.  Only
        # shell_inv_fm2 has already combined them with q^2/(2*pi^2).
        values = Dict(c => sum(weights[i] * (Float64(qnodes[i])^2 / (2pi^2)) *
                              Float64(getproperty(group[i], c)) for i in eachindex(group))
                      for c in COMPONENTS)
        total = sum(values[c] for c in COMPONENTS)
        e, channel, route = key
        push!(out, (sqrt_s_NN_GeV=e, channel=channel, route=route,
            bound_density_inv_fm3=values[:bound], unitary_density_inv_fm3=values[:unitary],
            landau_density_inv_fm3=values[:landau], total_density_inv_fm3=total,
            bound_fraction=total == 0 ? NaN : values[:bound] / total,
            unitary_fraction=total == 0 ? NaN : values[:unitary] / total,
            landau_fraction=total == 0 ? NaN : values[:landau] / total,
            q_nodes=nq, qmax_inv_fm=qmax, passed=all(Bool(r.passed) for r in group),
            production_authorized=false))
    end
    return out
end

function _density_map(rows)
    Dict((Float64(r.sqrt_s_NN_GeV), String(r.channel)) => Float64(r.density_inv_fm3)
         for r in rows)
end

function ratio_comparison(direct_density, reference_density)
    d = _density_map(direct_density)
    r = _density_map(reference_density)
    out = NamedTuple[]
    for energy in (200.0, 7.7, 3.0)
        dp = d[(energy, "pi_plus")]; dm = d[(energy, "pi_minus")]
        dk = d[(energy, "K_plus")]; dkm = d[(energy, "K_minus")]
        rp = r[(energy, "pi_plus")]; rm = r[(energy, "pi_minus")]
        rk = r[(energy, "K_plus")]; rkm = r[(energy, "K_minus")]
        push!(out, (sqrt_s_NN_GeV=energy,
            direct_Kplus_over_pi_plus=dk/dp, reference_Kplus_over_pi_plus=rk/rp,
            direct_Kminus_over_pi_minus=dkm/dm, reference_Kminus_over_pi_minus=rkm/rm,
            Kplus_density_factor_reference_over_direct=rk/dk,
            Kminus_density_factor_reference_over_direct=rkm/dkm,
            pi_plus_density_factor_reference_over_direct=rp/dp,
            pi_minus_density_factor_reference_over_direct=rm/dm,
            production_authorized=false))
    end
    return out
end

function _historical_context(path)
    rows = _rows(path)
    selected = [r for r in rows if String(r.regime) == "phase_shift_gbu_reference" &&
        isfinite(Float64(r.kpi_ratio)) && String(r.status) == "ok"]
    isempty(selected) && error("no usable historical GBU rows")
    out = NamedTuple[]
    for T in sort(unique(Float64(r.T_MeV) for r in selected))
        g = [r for r in selected if Float64(r.T_MeV) == T]
        vals = Float64.(getproperty.(g, :kpi_ratio))
        push!(out, (T_MeV=T, rho_min=minimum(Float64.(getproperty.(g, :rho_target))),
            rho_max=maximum(Float64.(getproperty.(g, :rho_target))),
            ratio_min=minimum(vals), ratio_max=maximum(vals), ratio_median=median(vals),
            muB_min_MeV=minimum(Float64.(getproperty.(g, :muB_MeV))),
            muB_max_MeV=maximum(Float64.(getproperty.(g, :muB_MeV))),
            mu_u_min_MeV=minimum(Float64.(getproperty.(g, :mu_u_MeV))),
            mu_u_max_MeV=maximum(Float64.(getproperty.(g, :mu_u_MeV))),
            mu_d_min_MeV=minimum(Float64.(getproperty.(g, :mu_d_MeV))),
            mu_d_max_MeV=maximum(Float64.(getproperty.(g, :mu_d_MeV))),
            mu_s_min_MeV=minimum(Float64.(getproperty.(g, :mu_s_MeV))),
            mu_s_max_MeV=maximum(Float64.(getproperty.(g, :mu_s_MeV))),
            constraint_mode=String(first(g).constraint_mode), path_strategy=String(first(g).path_strategy),
            density_policy=String(first(g).density_policy), bose_x_min=Float64(first(g).bose_x_min),
            route=String(first(g).regime), matched_to_current_bqs=false,
            comparison_status="historical_context_unmatched_background", production_authorized=false))
    end
    return out
end

function _render(current, historical, path)
    energies = Float64.(getproperty.(current, :sqrt_s_NN_GeV))
    p1 = plot(; xlabel="sqrt(s_NN) [GeV]", ylabel="K/pi ratio", xscale=:log10,
        title="Current FixedMuBConservedCharges BQS", legend=:outertopright, grid=true, linewidth=2)
    plot!(p1, energies, Float64.(getproperty.(current, :direct_Kplus_over_pi_plus));
        marker=:circle, label="direct finite-q K+/pi+", color=:blue)
    plot!(p1, energies, Float64.(getproperty.(current, :reference_Kplus_over_pi_plus));
        marker=:diamond, label="q=0 lambda reference K+/pi+", color=:navy)
    plot!(p1, energies, Float64.(getproperty.(current, :direct_Kminus_over_pi_minus));
        marker=:circle, linestyle=:dash, label="direct finite-q K-/pi-", color=:red)
    plot!(p1, energies, Float64.(getproperty.(current, :reference_Kminus_over_pi_minus));
        marker=:diamond, linestyle=:dash, label="q=0 lambda reference K-/pi-", color=:darkred)

    p2 = plot(; xlabel="sqrt(s_NN) [GeV]", ylabel="reference / direct density", xscale=:log10,
        title="Route density factor", legend=:outertopright, grid=true, linewidth=2)
    plot!(p2, energies, Float64.(getproperty.(current, :pi_plus_density_factor_reference_over_direct));
        marker=:circle, label="pi+", color=:blue)
    plot!(p2, energies, Float64.(getproperty.(current, :Kplus_density_factor_reference_over_direct));
        marker=:circle, label="K+", color=:green)
    plot!(p2, energies, Float64.(getproperty.(current, :pi_minus_density_factor_reference_over_direct));
        marker=:diamond, label="pi-", color=:red)
    plot!(p2, energies, Float64.(getproperty.(current, :Kminus_density_factor_reference_over_direct));
        marker=:diamond, label="K-", color=:orange)

    ht = Float64.(getproperty.(historical, :T_MeV))
    p3 = plot(; xlabel="T [MeV]", ylabel="historical K+/pi+", yscale=:log10,
        title="FixedAsymmetricRho context (unmatched)", legend=:outertopright, grid=true, linewidth=2)
    plot!(p3, ht, Float64.(getproperty.(historical, :ratio_median));
        marker=:diamond, label="median over rho", color=:purple)
    plot!(p3, ht, Float64.(getproperty.(historical, :ratio_min));
        linestyle=:dash, label="min/max envelope", color=:purple)
    plot!(p3, ht, Float64.(getproperty.(historical, :ratio_max)); label="", color=:purple)
    fig = plot(p1, p2, p3; layout=(1, 3), size=(1800, 520), margin=5Plots.mm,
        plot_title="GBU freezeout route audit (diagnostic; no production authorization)")
    savefig(fig, path)
end

function main()
    base = joinpath(R.ROOT, "data", "outputs", "results", "relaxtime", "analysis", "charged_rpa_phase_backend")
    direct_dir = get(ENV, "GBU_DIRECT_DIR", joinpath(base, "method_v1_saved_sparse"))
    reference_dir = get(ENV, "GBU_REFERENCE_DIR", joinpath(base, "method_v1_freezeout_q0_reference_v3"))
    historical = get(ENV, "GBU_HISTORICAL_CSV", joinpath(R.ROOT, "data", "outputs", "results", "relaxtime", "meson_density", "trho_asymmetric_kplus_piplus_scan_v1", "combined_meson_density_scan.csv"))
    output = get(ENV, "GBU_COMPARISON_OUTPUT", joinpath(base, "method_v1_freezeout_comparison_v2"))
    ispath(output) && error("refusing to overwrite $(output)")
    mkpath(output)
    direct_shell = _rows(joinpath(direct_dir, "refined_shells.csv"))
    reference_shell = _rows(joinpath(reference_dir, "reference_shells.csv"))
    direct_density = _rows(joinpath(direct_dir, "refined_densities.csv"))
    reference_density = _rows(joinpath(reference_dir, "reference_densities.csv"))
    # Direct density rows use no route column; bind the route explicitly for the
    # ratio helper while retaining source CSVs unchanged.
    direct_density_bound = [(merge((route="direct_finite_q",), NamedTuple(r))) for r in direct_density]
    reference_density_bound = [(merge((route="q0_lambda_reference",), NamedTuple(r))) for r in reference_density]
    direct_components = aggregate_components(direct_shell; qmax=8.0, nq=24)
    reference_components = aggregate_components(reference_shell; qmax=8.0, nq=24)
    current = ratio_comparison(direct_density_bound, reference_density_bound)
    hist = _historical_context(historical)
    CSV.write(joinpath(output, "component_totals.csv"), vcat(direct_components, reference_components))
    CSV.write(joinpath(output, "ratio_comparison.csv"), current)
    CSV.write(joinpath(output, "historical_context.csv"), hist)
    _render(current, hist, joinpath(output, "freezeout_route_audit.png"))
    manifest = Dict(
        "schema" => "charged_gbu_freezeout_comparison_v1",
        "status" => "diagnostic_only_unmatched_historical_context",
        "direct_input" => abspath(direct_dir), "reference_input" => abspath(reference_dir),
        "historical_input" => abspath(historical),
        "input_hashes" => Dict(
            "direct_shells" => R.hashfile(joinpath(direct_dir, "refined_shells.csv")),
            "direct_densities" => R.hashfile(joinpath(direct_dir, "refined_densities.csv")),
            "reference_shells" => R.hashfile(joinpath(reference_dir, "reference_shells.csv")),
            "reference_densities" => R.hashfile(joinpath(reference_dir, "reference_densities.csv")),
            "historical_csv" => R.hashfile(historical)),
        "routes" => ["direct_finite_q", "q0_lambda_reference"],
        "background" => "FixedMuBConservedCharges quark-only BQS; rhoQ/rhoB=0.4; rhoS=0",
        "observable" => "fixed_quark_only_gbu_partial_yield",
        "components" => ["bound", "unitary", "landau"],
        "component_integration" => "same 24-node Gauss-Legendre q weights as retained runs; shell already contains q^2/(2pi^2)",
        "historical_comparison" => "FixedAsymmetricRho, x_min_cut, bose_x_min=0.01; context only, no matched overlay",
        "git_head" => readchomp(`git -C $(R.ROOT) rev-parse HEAD`),
        "production_authorized" => false,
        "output_hashes" => Dict(f => R.hashfile(joinpath(output, f)) for f in ("component_totals.csv", "ratio_comparison.csv", "historical_context.csv", "freezeout_route_audit.png")))
    open(joinpath(output, "manifest.json"), "w") do io
        JSON3.write(io, manifest)
    end
    println("[gbu-freezeout-comparison] output: $(output)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
end
