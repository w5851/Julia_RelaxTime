module AFieldBuilder

using ..GaussLegendre: gauleg, DEFAULT_COSΘ_NODES
using ..OneLoopIntegrals: A
using ..OneLoopIntegralsCorrection: A_aniso

export build_A_triplet, ensure_quark_params_has_A

@inline _xi(thermo_params::NamedTuple)::Float64 = hasproperty(thermo_params, :ξ) ? Float64(thermo_params.ξ) : 0.0

@inline function _validate_inputs(quark_params::NamedTuple, thermo_params::NamedTuple)
    hasproperty(quark_params, :m) || error("quark_params is missing :m")
    hasproperty(quark_params, :μ) || error("quark_params is missing :μ")
    hasproperty(thermo_params, :T) || error("thermo_params is missing :T")
    hasproperty(thermo_params, :Φ) || error("thermo_params is missing :Φ")
    hasproperty(thermo_params, :Φbar) || error("thermo_params is missing :Φbar")
    return nothing
end

function build_A_triplet(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    nodes_p::AbstractVector,
    weights_p::AbstractVector;
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    use_aniso::Bool=true,
)::NamedTuple
    _validate_inputs(quark_params, thermo_params)

    length(nodes_p) == length(weights_p) || error("build_A_triplet: nodes_p and weights_p length mismatch")

    ξ = _xi(thermo_params)
    if use_aniso && abs(ξ) > 0.0
        nodes_cos, weights_cos = gauleg(-1.0, 1.0, cos_nodes)
        A_u = A_aniso(quark_params.m.u, quark_params.μ.u, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
                      ξ, collect(Float64, nodes_p), collect(Float64, weights_p), nodes_cos, weights_cos)
        A_d = A_aniso(quark_params.m.d, quark_params.μ.d, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
                      ξ, collect(Float64, nodes_p), collect(Float64, weights_p), nodes_cos, weights_cos)
        A_s = A_aniso(quark_params.m.s, quark_params.μ.s, thermo_params.T, thermo_params.Φ, thermo_params.Φbar,
                      ξ, collect(Float64, nodes_p), collect(Float64, weights_p), nodes_cos, weights_cos)
        return (u=A_u, d=A_d, s=A_s)
    end

    p_nodes_f64 = collect(Float64, nodes_p)
    p_weights_f64 = collect(Float64, weights_p)
    A_u = A(quark_params.m.u, quark_params.μ.u, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, p_nodes_f64, p_weights_f64)
    A_d = A(quark_params.m.d, quark_params.μ.d, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, p_nodes_f64, p_weights_f64)
    A_s = A(quark_params.m.s, quark_params.μ.s, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, p_nodes_f64, p_weights_f64)
    return (u=A_u, d=A_d, s=A_s)
end

function build_A_triplet(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    p_nodes::Int=16,
    p_max::Float64=20.0,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    use_aniso::Bool=true,
)::NamedTuple
    nodes_p, weights_p = gauleg(0.0, p_max, p_nodes)
    return build_A_triplet(
        quark_params,
        thermo_params,
        nodes_p,
        weights_p;
        cos_nodes=cos_nodes,
        use_aniso=use_aniso,
    )
end

function ensure_quark_params_has_A(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    p_nodes::Int=16,
    p_max::Float64=20.0,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    use_aniso::Bool=true,
    warn_on_auto::Bool=true,
)::NamedTuple
    if hasproperty(quark_params, :A)
        return quark_params
    end

    ξ = _xi(thermo_params)
    if warn_on_auto && abs(ξ) > 0.0
        if use_aniso
            @warn "quark_params missing :A under ξ≠0; auto-building anisotropic A via A_aniso" ξ=ξ p_nodes=p_nodes p_max=p_max cos_nodes=cos_nodes
        else
            @warn "quark_params missing :A under ξ≠0; use_aniso=false so isotropic A will be used" ξ=ξ p_nodes=p_nodes p_max=p_max
        end
    end

    A_vals = build_A_triplet(
        quark_params,
        thermo_params;
        p_nodes=p_nodes,
        p_max=p_max,
        cos_nodes=cos_nodes,
        use_aniso=use_aniso,
    )
    return merge(quark_params, (A=A_vals,))
end

end # module AFieldBuilder
