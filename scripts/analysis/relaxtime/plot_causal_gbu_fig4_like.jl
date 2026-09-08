"""Render FIG.4-like ratio diagnostics for current and historical routes.

The current panel uses explicit FixedMuBConservedCharges freezeout points.
The historical panel is parameterized by T because the retained
FixedAsymmetricRho artifact is a T-rho scan, not a sqrt(s_NN) freezeout curve.
"""
module CausalGBUFig4Like

using CSV, JSON3, Plots, SHA

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

function rows(path)
    isfile(path) || error("missing input $(path)")
    collect(CSV.File(path; comment="#"))
end

function current_plot(path)
    r = rows(path)
    sort!(r; by=x -> Float64(x.sqrt_s_NN_GeV))
    x = Float64.(getproperty.(r, :sqrt_s_NN_GeV))
    finite_values = filter(isfinite, vcat(
        Float64.(getproperty.(r, :direct_Kplus_over_pi_plus)),
        Float64.(getproperty.(r, :reference_Kplus_over_pi_plus)),
        Float64.(getproperty.(r, :direct_Kminus_over_pi_minus)),
        Float64.(getproperty.(r, :reference_Kminus_over_pi_minus))))
    isempty(finite_values) && error("no finite ratios to display")
    ymax = max(0.62, 1.15maximum(finite_values))
    ymin = min(0.0, 1.15minimum(finite_values))
    p = plot(; xlabel="sqrt(s_NN) [GeV]", ylabel="K/pi ratio", xscale=:log10,
        ylims=(ymin, ymax), title="Current quark-only BQS GBU partial yield ($(length(r)) energies)",
        grid=true, legend=:outertopright, linewidth=2.2, size=(1200, 700),
        bottom_margin=8Plots.mm, top_margin=5Plots.mm)
    for (field,gate,marker,color,style,label) in (
        (:direct_Kplus_over_pi_plus,:direct_plus_passed,:circle,:black,:solid,"K+/pi+ direct finite-q"),
        (:reference_Kplus_over_pi_plus,:reference_plus_passed,:diamond,:black,:dash,"K+/pi+ q=0 reference"),
        (:direct_Kminus_over_pi_minus,:direct_minus_passed,:circle,:red3,:solid,"K-/pi- direct finite-q"),
        (:reference_Kminus_over_pi_minus,:reference_minus_passed,:diamond,:red3,:dash,"K-/pi- q=0 reference"))
        raw=Float64.(getproperty.(r,field))
        accepted=[!hasproperty(row,gate) || Bool(getproperty(row,gate)) for row in r]
        plot!(p,x,[ok ? v : NaN for (ok,v) in zip(accepted,raw)];marker=marker,color=color,linestyle=style,label=label)
        failed=findall(i->!accepted[i] && isfinite(raw[i]),eachindex(raw))
        isempty(failed) || scatter!(p,x[failed],raw[failed];marker=:xcross,color=:orange,label="failed: $(label)")
    end
    annotate!(p, minimum(x)*1.04, 0.96ymax,
        text("rhoQ/rhoB=0.4, rhoS=0; no meson feedback\nDiagnostic: line segments only guide the eye", 9, :left))
    return p
end

function historical_plot(path)
    r = rows(path)
    sort!(r; by=x -> Float64(x.T_MeV))
    x = Float64.(getproperty.(r, :T_MeV))
    lo = Float64.(getproperty.(r, :ratio_min))
    med = Float64.(getproperty.(r, :ratio_median))
    hi = Float64.(getproperty.(r, :ratio_max))
    p = plot(; xlabel="T [MeV]", ylabel="K+/pi+ ratio", yscale=:log10,
        title="Historical FixedAsymmetricRho context", grid=true,
        legend=:outertopright, linewidth=2.2, size=(900, 560))
    plot!(p, x, med; color=:purple, marker=:diamond, label="median over rho")
    plot!(p, x, lo; color=:purple, linestyle=:dash, label="rho envelope")
    plot!(p, x, hi; color=:purple, linestyle=:dash, label="")
    plot!(p, x, lo; fillrange=hi, fillalpha=0.12, fillcolor=:purple, linealpha=0,
        label="rho=0.05..1.00 envelope")
    annotate!(p, 121, maximum(hi) / 1.5,
        text("Unmatched: FixedAsymmetricRho, rho_u/rho_d=0.876, x_min_cut, bose_x_min=0.01", 8, :left))
    return p
end

function main()
    base = joinpath(ROOT, "data", "outputs", "results", "relaxtime", "analysis",
        "charged_rpa_phase_backend", "method_v1_freezeout_comparison_v2")
    current = get(ENV, "GBU_FIG4_CURRENT", joinpath(base, "ratio_comparison.csv"))
    historical = get(ENV, "GBU_FIG4_HISTORICAL", joinpath(base, "historical_context.csv"))
    out = get(ENV, "GBU_FIG4_OUTPUT",
        joinpath(ROOT, "data", "outputs", "results", "relaxtime", "analysis",
            "charged_rpa_phase_backend", "fig4_like_freezeout_ratio_20260905"))
    ispath(out) && error("refusing to overwrite $(out)")
    mkpath(out)
    p1 = current_plot(current)
    p2 = historical_plot(historical)
    combined = plot(p1, p2; layout=(1, 2), size=(1840, 610), margin=5Plots.mm,
        plot_title="FIG.4-like charged GBU ratio audit - diagnostic only")
    savefig(p1, joinpath(out, "current_bqs_fig4_like.png"))
    savefig(p2, joinpath(out, "historical_fixedasymrho_context.png"))
    savefig(combined, joinpath(out, "fig4_like_current_vs_historical.png"))
    report = Dict(
        "schema" => "charged_gbu_fig4_like_v1",
        "paper_reference" => "Blaschke et al., Particles 3 (2020) 169-181, DOI 10.3390/particles3010014",
        "current_input" => abspath(current), "historical_input" => abspath(historical),
        "current_plot" => "current_bqs_fig4_like.png",
        "historical_plot" => "historical_fixedasymrho_context.png",
        "combined_plot" => "fig4_like_current_vs_historical.png",
        "current_x" => "sqrt(s_NN)_GeV; explicit p48 freezeout points",
        "current_point_count" => length(rows(current)),
        "failed_point_policy" => "NaN line gaps; finite failed ratios shown as orange crosses",
        "historical_x" => "T_MeV; median and rho envelope of retained T-rho artifact",
        "current_background" => "FixedMuBConservedCharges quark-only BQS; rhoQ/rhoB=0.4; rhoS=0",
        "historical_background" => "FixedAsymmetricRho; rho_u/rho_d=0.876; rho_s=0",
        "historical_policy" => "x_min_cut; bose_x_min=0.01",
        "comparability" => "unmatched_context_only; no interpolation or quantitative overlay",
        "production_authorized" => false,
        "source_script_sha256" => bytes2hex(SHA.sha256(read(@__FILE__))),
        "output_hashes" => Dict(f => bytes2hex(SHA.sha256(read(joinpath(out,f)))) for f in
            ("current_bqs_fig4_like.png","historical_fixedasymrho_context.png","fig4_like_current_vs_historical.png")),
        "input_hashes" => Dict("current" => bytes2hex(SHA.sha256(read(current))),
                              "historical" => bytes2hex(SHA.sha256(read(historical)))))
    open(joinpath(out, "manifest.json"), "w") do io
        JSON3.write(io, report)
    end
    open(joinpath(out, "README.md"), "w") do io
        write(io, "This directory is diagnostic only. The current panel uses $(length(rows(current))) explicit FixedMuBConservedCharges quark-only BQS points and distinguishes direct finite-q from the q=0 lambda-reference route. Lines guide the eye; they are not extra solved points. Failed points break lines. The historical panel uses the FixedAsymmetricRho T-rho scan as a temperature envelope; it is not a chemical-freezeout sqrt(s_NN) curve and must not be quantitatively overlaid on the current panel.\n")
    end
    println("[gbu-fig4-like] output: $(out)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
end
