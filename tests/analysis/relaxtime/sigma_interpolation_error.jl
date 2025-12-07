#!/usr/bin/env julia

using Printf
using LinearAlgebra

try
    using Plots
catch
    import Pkg
    Pkg.add("Plots")
    using Plots
end

try
    using Interpolations
catch
    import Pkg
    Pkg.add("Interpolations")
    using Interpolations
end

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

include(joinpath(PROJECT_ROOT, "tests", "analysis", "relaxtime", "integrand_vs_s_average_rate.jl"))

# ------------ Config ------------
process_heavy_init = :ssbar_to_uubar   # 初态更重，左端点有 ~1/sqrt 行为
process_heavy_final = :uubar_to_ssbar  # 初态不更重（末态更重），期望线性更合适

n_sigma_points = 64            # 截面积分点（确保参考精度）
reference_nodes = 200          # 参考网格点数（dense 真值）
sample_counts = [4, 6, 8, 12, 16, 24]  # 着重少点高精度

params = build_params()

function s_threshold(process::Symbol)
    pi_sym, pj_sym, _, _ = ASR.parse_particles_from_process(process)
    mi = ASR.get_mass(pi_sym, params.quark_params)
    mj = ASR.get_mass(pj_sym, params.quark_params)
    return (mi + mj)^2
end

function s_threshold_final(process::Symbol)
    _, _, pc_sym, pd_sym = ASR.parse_particles_from_process(process)
    mc = ASR.get_mass(pc_sym, params.quark_params)
    md = ASR.get_mass(pd_sym, params.quark_params)
    return (mc + md)^2
end

function build_s_grid(process::Symbol; δ=1e-4, width=10.0, n=reference_nodes)
    s0 = s_threshold(process)
    sf = s_threshold_final(process)
    s_start = max(s0, sf) + δ
    return collect(range(s_start, s_start + width; length=n))
end

"""Chebyshev-like节点映射到 [a,b]。"""
function cheb_nodes(n::Int, a::Float64, b::Float64)
    k = collect(0:n-1)
    x = cos.((2 .* k .+ 1) .* π ./ (2n))
    return 0.5 .* ((b - a) .* x .+ (b + a))
end

"""指数型左端聚焦节点，lambda>0 越大越偏左。"""
function left_biased_nodes(n::Int, a::Float64, b::Float64; lambda::Float64=3.0)
    t = range(0.0, 1.0; length=n)
    w = (exp.(lambda .* t) .- 1.0) ./ (exp(lambda) - 1.0)
    return a .+ (b - a) .* w
end

function compute_sigma(process::Symbol, s_values::Vector{Float64})
    cache = ASR.CrossSectionCache(process)
    out = Float64[]
    for s in s_values
        push!(out, ASR.get_sigma(cache, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_sigma_points))
    end
    return out
end

function fit_model(model::Symbol, s::Vector{Float64}, y::Vector{Float64}, s0::Float64)
    shift = s .- s0
    if model == :sqrt1
        X = hcat(1.0 ./ sqrt.(shift), shift, ones(length(s)))              # A/√ + B*(s-s0) + C
        return X \ y
    elseif model == :sqrt2
        X = hcat(1.0 ./ sqrt.(shift), shift, shift.^2, ones(length(s)))    # A/√ + B*(s-s0) + C*(s-s0)^2 + D
        return X \ y
    elseif model == :rsqrt
        X = hcat(1.0 ./ sqrt.(shift), shift, ones(length(s)))              # A/√ + rational tail via C on shift in predict
        return X \ y
    elseif model == :poly3
        X = hcat(shift, shift.^2, shift.^3, ones(length(s)))               # cubic around s0
        return X \ y
    elseif model == :poly2
        X = hcat(shift, shift.^2, ones(length(s)))                         # quadratic around s0
        return X \ y
    elseif model == :pade11
        # (a0 + a1*x) / (1 + b1*x); linearize: a0 + a1*x - y - b1*x*y = 0
        X = hcat(ones(length(s)), shift, -(shift .* y))
        return X \ y  # returns (a0, a1, b1)
    elseif model == :spline_lin
        knots = (s,)
        itp = Interpolations.interpolate(knots, y, Gridded(Linear()))
        return Interpolations.extrapolate(itp, Interpolations.Line())
    elseif model == :logpade11
        epsy = max(eps(), minimum(abs.(y[y .!= 0.0]); init=Inf))
        ly = log.(abs.(y) .+ epsy)
        X = hcat(ones(length(s)), shift, -(shift .* ly))
        return (X \ ly, epsy)
    else
        error("unknown model $model")
    end
end

function predict_model(model::Symbol, coeff, s::Vector{Float64}, s0::Float64)
    shift = s .- s0
    if model == :sqrt1
        A, B, C = coeff
        return A .* (1.0 ./ sqrt.(shift)) .+ B .* shift .+ C
    elseif model == :sqrt2
        A, B, C, D = coeff
        return A .* (1.0 ./ sqrt.(shift)) .+ B .* shift .+ C .* shift.^2 .+ D
    elseif model == :rsqrt
        A, B, C = coeff
        return A .* (1.0 ./ sqrt.(shift)) .+ (B .+ C .* shift) ./ (1 .+ shift)
    elseif model == :poly3
        a, b, c, d = coeff
        return a .* shift .+ b .* shift.^2 .+ c .* shift.^3 .+ d
    elseif model == :poly2
        a, b, c = coeff
        return a .* shift .+ b .* shift.^2 .+ c
    elseif model == :pade11
        a0, a1, b1 = coeff
        return (a0 .+ a1 .* shift) ./ (1 .+ b1 .* shift)
    elseif model == :spline_lin
        return coeff.(s)
    elseif model == :logpade11
        (p, epsy) = coeff
        a0, a1, b1 = p
        ly_pred = (a0 .+ a1 .* shift) ./ (1 .+ b1 .* shift)
        return exp.(ly_pred) .- epsy
    else
        error("unknown model $model")
    end
end

function metrics(y_true::Vector{Float64}, y_pred::Vector{Float64})
    n_true = norm(y_true)
    if n_true == 0.0
        return (rel_l2=0.0, rel_max=0.0)
    end
    rel_l2 = norm(y_pred .- y_true) / n_true
    denom = @. max(abs(y_true), eps())
    rel_max = maximum(abs.((y_pred .- y_true) ./ denom))
    return (rel_l2=rel_l2, rel_max=rel_max)
end

function analyze_case(process::Symbol; model::Symbol, biased::Bool, lambda::Float64=3.0)
    s0 = s_threshold(process)
    s_grid = build_s_grid(process)
    σ_ref = compute_sigma(process, s_grid)

    @printf("\nProcess: %s (model=%s, biased=%s)\n", string(process), string(model), string(biased))
    @printf("s0 = %.6f\n", s0)

    for n in sample_counts
        s_sample = biased ? sort(left_biased_nodes(n, s_grid[1], s_grid[end]; lambda=lambda)) : sort(cheb_nodes(n, s_grid[1], s_grid[end]))
        σ_sample = compute_sigma(process, s_sample)

        coeff = fit_model(model, s_sample, σ_sample, s0)
        σ_pred = predict_model(model, coeff, s_grid, s0)

        m = metrics(σ_ref, σ_pred)
        @printf("  n=%-3d -> rel_l2=%.3e rel_max=%.3e\n", n, m.rel_l2, m.rel_max)
    end
end

function best_model_for_process(process::Symbol; biased::Bool, lambda::Float64=3.0)
    models = process == process_heavy_init ? [:sqrt1, :sqrt2, :rsqrt, :poly2] : [:poly2, :poly3, :pade11, :logpade11, :spline_lin]
    s0 = s_threshold(process)
    s_grid = build_s_grid(process)
    σ_ref = compute_sigma(process, s_grid)
    best = nothing
    for model in models
        for n in sample_counts
            s_sample = biased ? sort(left_biased_nodes(n, s_grid[1], s_grid[end]; lambda=lambda)) : sort(cheb_nodes(n, s_grid[1], s_grid[end]))
            σ_sample = compute_sigma(process, s_sample)
            coeff = fit_model(model, s_sample, σ_sample, s0)
            σ_pred = predict_model(model, coeff, s_grid, s0)
            m = metrics(σ_ref, σ_pred)
            if best === nothing || m.rel_l2 < best.rel_l2
                best = (; model, n, rel_l2=m.rel_l2, rel_max=m.rel_max, coeff, s_sample, σ_sample, σ_pred, s_grid, σ_ref, biased)
            end
        end
    end
    return best
end

function plot_fit(result, title_str, outfile)
    (; s_sample, σ_sample, s_grid, σ_ref, σ_pred, model, n) = result
    plt = plot(s_grid, σ_ref; label="reference", lw=2)
    plot!(plt, s_grid, σ_pred; label="fit $(string(model)) (n=$(n))", lw=2, ls=:dash)
    scatter!(plt, s_sample, σ_sample; label="samples", ms=5)
    xlabel!(plt, "s")
    ylabel!(plt, "σ(s)")
    title!(plt, title_str)
    outdir = joinpath(PROJECT_ROOT, "output", "figures")
    mkpath(outdir)
    savefig(plt, joinpath(outdir, outfile))
end

# Run analyses
println("Starting fits...")
lambda_heavy = 4.0   # 更强左聚焦
lambda_light = 4.0
analyze_case(process_heavy_init; model=:sqrt2, biased=true, lambda=lambda_heavy)
analyze_case(process_heavy_final; model=:spline_lin, biased=true, lambda=lambda_light)

best_heavy = best_model_for_process(process_heavy_init; biased=true, lambda=lambda_heavy)
best_light = best_model_for_process(process_heavy_final; biased=true, lambda=lambda_light)

@printf("\nBest heavy-init: model=%s n=%d rel_l2=%.3e rel_max=%.3e (biased=%s)\n", string(best_heavy.model), best_heavy.n, best_heavy.rel_l2, best_heavy.rel_max, string(best_heavy.biased))
@printf("Best heavy-final: model=%s n=%d rel_l2=%.3e rel_max=%.3e (biased=%s)\n", string(best_light.model), best_light.n, best_light.rel_l2, best_light.rel_max, string(best_light.biased))

plot_fit(best_heavy, "σ(s) fit (heavy init)", "sigma_fit_heavy.png")
plot_fit(best_light, "σ(s) fit (heavy final)", "sigma_fit_light.png")

println("Done. Figures in output/figures/.")
