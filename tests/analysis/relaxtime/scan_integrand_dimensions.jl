#!/usr/bin/env julia

try
    using Plots
catch
    import Pkg
    Pkg.add("Plots")
    using Plots
end

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

include(joinpath(PROJECT_ROOT, "tests", "analysis", "relaxtime", "integrand_vs_s_average_rate.jl"))

# 基本设置
process = :ssbar_to_uubar
p_nodes = 6
angle_nodes = 4
phi_nodes = 8
n_sigma_points = 16

params = build_params()
cache = ASR.CrossSectionCache(process)

# 参考节点与权重（固定其他维度）
p_grid, p_w = gauleg(0.0, Λ_inv_fm, p_nodes)
cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)
phi_grid, phi_w = gauleg(0.0, Float64(π), phi_nodes)

ref_pi_idx = cld(length(p_grid), 2)
ref_pj_idx = cld(length(p_grid), 2)
ref_cθi_idx = cld(length(cos_grid), 2)
ref_cθj_idx = cld(length(cos_grid), 2)
ref_φ_idx = cld(length(phi_grid), 2)

ref_pi, ref_w_pi = p_grid[ref_pi_idx], p_w[ref_pi_idx]
ref_pj, ref_w_pj = p_grid[ref_pj_idx], p_w[ref_pj_idx]
ref_cθi, ref_w_cθi = cos_grid[ref_cθi_idx], cos_w[ref_cθi_idx]
ref_cθj, ref_w_cθj = cos_grid[ref_cθj_idx], cos_w[ref_cθj_idx]
ref_φ, ref_w_φ = phi_grid[ref_φ_idx], phi_w[ref_φ_idx]

function weight_product(var::Symbol, w_var::Float64)
    if var == :p_i
        return w_var * ref_w_pj * ref_w_cθi * ref_w_cθj * ref_w_φ
    elseif var == :p_j
        return ref_w_pi * w_var * ref_w_cθi * ref_w_cθj * ref_w_φ
    elseif var == :cθi
        return ref_w_pi * ref_w_pj * w_var * ref_w_cθj * ref_w_φ
    elseif var == :cθj
        return ref_w_pi * ref_w_pj * ref_w_cθi * w_var * ref_w_φ
    elseif var == :φ
        return ref_w_pi * ref_w_pj * ref_w_cθi * ref_w_cθj * w_var
    else
        error("unknown variable $var")
    end
end

function eval_term(p_i::Float64, p_j::Float64, cθi::Float64, cθj::Float64, φ::Float64)
    rec = omega_term_at(process, p_i, p_j, cθi, cθj, φ; params=params, cache=cache, n_sigma_points=n_sigma_points)
    return rec === nothing ? nothing : rec.term
end

function sweep_variable(var::Symbol; n_dense::Int=140)
    if var == :p_i
        dense_x = collect(range(0.0, Λ_inv_fm, length=n_dense))
        node_x = p_grid
        node_w = p_w
        f = x -> eval_term(x, ref_pj, ref_cθi, ref_cθj, ref_φ)
    elseif var == :p_j
        dense_x = collect(range(0.0, Λ_inv_fm, length=n_dense))
        node_x = p_grid
        node_w = p_w
        f = x -> eval_term(ref_pi, x, ref_cθi, ref_cθj, ref_φ)
    elseif var == :cθi
        dense_x = collect(range(0.0, 1.0, length=n_dense))
        node_x = cos_grid
        node_w = cos_w
        f = x -> eval_term(ref_pi, ref_pj, x, ref_cθj, ref_φ)
    elseif var == :cθj
        dense_x = collect(range(0.0, 1.0, length=n_dense))
        node_x = cos_grid
        node_w = cos_w
        f = x -> eval_term(ref_pi, ref_pj, ref_cθi, x, ref_φ)
    elseif var == :φ
        dense_x = collect(range(0.0, Float64(π), length=n_dense))
        node_x = phi_grid
        node_w = phi_w
        f = x -> eval_term(ref_pi, ref_pj, ref_cθi, ref_cθj, x)
    else
        error("unknown variable $var")
    end

    dense_term = Float64[]
    for x in dense_x
        t = f(x)
        push!(dense_term, t === nothing ? 0.0 : t)
    end

    node_weighted = Float64[]
    node_term = Float64[]
    for (x, w) in zip(node_x, node_w)
        t = f(x)
        if t === nothing
            push!(node_term, 0.0)
            push!(node_weighted, 0.0)
        else
            push!(node_term, t)
            push!(node_weighted, t * weight_product(var, w))
        end
    end

    return (; dense_x, dense_term, node_x, node_term, node_weighted)
end

function plot_variable(var::Symbol, label_x::String)
    res = sweep_variable(var)
    plt = plot(res.dense_x, res.dense_term; label="term (raw)", xlabel=label_x, ylabel="term / weighted", lw=2)
    scatter!(plt, res.node_x, res.node_weighted; label="weighted @ GL nodes", ms=5)
    scatter!(plt, res.node_x, res.node_term; label="term @ GL nodes", ms=5, markerstrokecolor=:black, markershape=:diamond)
    return plt
end

outdir = joinpath(PROJECT_ROOT, "output", "figures")
mkpath(outdir)

plots = [
    (:p_i, "p_i [fm⁻¹]", "integrand_vs_p_i.png"),
    (:p_j, "p_j [fm⁻¹]", "integrand_vs_p_j.png"),
    (:cθi, "cosθ_i", "integrand_vs_cθi.png"),
    (:cθj, "cosθ_j", "integrand_vs_cθj.png"),
    (:φ, "φ", "integrand_vs_phi.png"),
]

for (var, label, fname) in plots
    plt = plot_variable(var, label)
    path = joinpath(outdir, fname)
    savefig(plt, path)
    @printf("Saved %s\n", path)
end

println("Done.")
