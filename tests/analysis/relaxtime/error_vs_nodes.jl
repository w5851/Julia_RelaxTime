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

# ---------- Config ----------
process = :ssbar_to_uubar
n_sigma_points = 16
ref_nodes = 64           # reference Gauss-Legendre nodes for 1D high-accuracy integral
test_nodes = [2, 3, 4, 6, 8, 10, 12, 16, 24, 32]

params = build_params()
cache = ASR.CrossSectionCache(process)

p_grid_default, p_w_default = gauleg(0.0, Λ_inv_fm, 6)
cos_grid_default, cos_w_default = gauleg(0.0, 1.0, 4)
phi_grid_default, phi_w_default = gauleg(0.0, Float64(π), 8)

ref_pi = p_grid_default[cld(length(p_grid_default), 2)]
ref_pj = p_grid_default[cld(length(p_grid_default), 2)]
ref_cθi = cos_grid_default[cld(length(cos_grid_default), 2)]
ref_cθj = cos_grid_default[cld(length(cos_grid_default), 2)]
ref_φ = phi_grid_default[cld(length(phi_grid_default), 2)]

ref_w_pi = p_w_default[cld(length(p_w_default), 2)]
ref_w_pj = p_w_default[cld(length(p_w_default), 2)]
ref_w_cθi = cos_w_default[cld(length(cos_w_default), 2)]
ref_w_cθj = cos_w_default[cld(length(cos_w_default), 2)]
ref_w_φ = phi_w_default[cld(length(phi_w_default), 2)]

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
    return rec === nothing ? 0.0 : rec.term
end

function ranges_for(var::Symbol, n::Int)
    if var == :p_i || var == :p_j
        return gauleg(0.0, Λ_inv_fm, n)
    elseif var == :cθi || var == :cθj
        return gauleg(0.0, 1.0, n)
    elseif var == :φ
        return gauleg(0.0, Float64(π), n)
    else
        error("unknown variable $var")
    end
end

function integrand_1d(var::Symbol, x::Float64)
    if var == :p_i
        return eval_term(x, ref_pj, ref_cθi, ref_cθj, ref_φ)
    elseif var == :p_j
        return eval_term(ref_pi, x, ref_cθi, ref_cθj, ref_φ)
    elseif var == :cθi
        return eval_term(ref_pi, ref_pj, x, ref_cθj, ref_φ)
    elseif var == :cθj
        return eval_term(ref_pi, ref_pj, ref_cθi, x, ref_φ)
    elseif var == :φ
        return eval_term(ref_pi, ref_pj, ref_cθi, ref_cθj, x)
    else
        error("unknown variable $var")
    end
end

function integrate_1d(var::Symbol, n::Int)
    xs, ws = ranges_for(var, n)
    total = 0.0
    for (x, w) in zip(xs, ws)
        total += weight_product(var, w) * integrand_1d(var, x)
    end
    return total
end

function relative_error(approx::Float64, reference::Float64)
    if reference == 0.0
        return approx == 0.0 ? 0.0 : Inf
    else
        return abs(approx - reference) / abs(reference)
    end
end

outdir = joinpath(PROJECT_ROOT, "output", "figures")
mkpath(outdir)

vars = [(:p_i, "p_i [fm⁻¹]"), (:p_j, "p_j [fm⁻¹]"), (:cθi, "cosθ_i"), (:cθj, "cosθ_j"), (:φ, "φ")]

for (var, label) in vars
    ref_val = integrate_1d(var, ref_nodes)
    errs = Float64[]
    for n in test_nodes
        val = integrate_1d(var, n)
        push!(errs, relative_error(val, ref_val))
    end
    @printf("%s reference (n=%d): %.6e\n", string(var), ref_nodes, ref_val)
    for (n, e) in zip(test_nodes, errs)
        @printf("  n=%-3d -> rel err=%.3e\n", n, e)
    end

    plt = plot(test_nodes, errs; xlabel="nodes", ylabel="relative error vs ref", yscale=:log10,
        marker=:circle, lw=2, title="Error decay for $(string(var))")
    savefig(plt, joinpath(outdir, "error_vs_nodes_$(string(var)).png"))
end

println("Done. Figures saved to output/figures/error_vs_nodes_*.png")
