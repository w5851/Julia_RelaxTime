module PhaseSpaceSampling

"""Model-agnostic phase-space sampling/integration helpers.

Design intent:
- Keep numeric grid generation and phase-space summation logic in one place.
- Let physics/model code provide integrands (dispersion/distribution/etc).

This module assumes it is included from a parent module that already includes
`integration/GaussLegendre.jl` as a submodule `GaussLegendre`.
"""

using ..GaussLegendre: gauleg

export p_nodes_weights, cos_nodes_weights
export integrate_p, integrate_p_cos

@inline function p_nodes_weights(p_nodes::Int, p_max::Float64, p_grid, p_w)
    if p_grid !== nothing
        p_w !== nothing || error("p_w must be provided when p_grid is provided")
        length(p_grid) == length(p_w) || error("p_grid and p_w length mismatch")
        return p_grid, p_w
    end
    return gauleg(0.0, p_max, p_nodes)
end

@inline function cos_nodes_weights(cos_nodes::Int, cos_grid, cos_w)
    if cos_grid !== nothing
        cos_w !== nothing || error("cos_w must be provided when cos_grid is provided")
        length(cos_grid) == length(cos_w) || error("cos_grid and cos_w length mismatch")
        return cos_grid, cos_w
    end
    return gauleg(-1.0, 1.0, cos_nodes)
end

@inline function integrate_p(
    integrand::Function,
    nodes_p::AbstractVector,
    weights_p::AbstractVector,
)::Float64
    acc = 0.0
    @inbounds for (p, wp) in zip(nodes_p, weights_p)
        acc += wp * integrand(p)
    end
    return acc
end

@inline function integrate_p(
    nodes_p::AbstractVector,
    weights_p::AbstractVector,
    integrand::Function,
)::Float64
    return integrate_p(integrand, nodes_p, weights_p)
end

@inline function integrate_p_cos(
    integrand::Function,
    nodes_p::AbstractVector,
    weights_p::AbstractVector,
    nodes_cos::AbstractVector,
    weights_cos::AbstractVector,
)::Float64
    acc = 0.0
    @inbounds for (p, wp) in zip(nodes_p, weights_p)
        for (c, wc) in zip(nodes_cos, weights_cos)
            acc += wp * wc * integrand(p, c)
        end
    end
    return acc
end

@inline function integrate_p_cos(
    nodes_p::AbstractVector,
    weights_p::AbstractVector,
    nodes_cos::AbstractVector,
    weights_cos::AbstractVector,
    integrand::Function,
)::Float64
    return integrate_p_cos(integrand, nodes_p, weights_p, nodes_cos, weights_cos)
end

end # module
