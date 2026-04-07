using StaticArrays: SVector

@inline function _model_kind_for_shared_core(model::AbstractQCDModel)
    return model isa RPNJLModel ? :RPNJL : :PNJL
end

@inline function _implicit_conditions_with_shared_core(
    model::AbstractQCDModel,
    θ::AbstractVector,
    x::AbstractVector;
    xi,
    p_num::Int,
    t_num::Int,
    kwargs...,
)
    state_n = gap_state_dim(model)
    if state_n == 5 && (length(θ) == 2 || length(θ) == 4)
        T_fm = θ[1]
        mu_vec = if length(θ) == 2
            SVector{3}(θ[2], θ[2], θ[2])
        else
            SVector{3}(θ[2], θ[3], θ[4])
        end
        Tx = typeof(x[1])
        x_state = SVector{5, Tx}(Tuple(x))
        params = GapParams(T_fm, cached_nodes(p_num, t_num), xi;
            p_num=p_num,
            t_num=t_num,
            model_kind=_model_kind_for_shared_core(model),
        )
        Tout = promote_type(Tx, typeof(T_fm), typeof(mu_vec[1]))
        out = Vector{Tout}(undef, 5)
        gap_core_residual!(out, x_state, mu_vec, params)
        return out
    end

    T_fm = θ[1]
    μ_fm = θ[2]
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    r = gap_residual(model, x, T_fm, mu_vec;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        kwargs...)
    return Vector(r)
end

"""
    ImplicitAdapters

将现有 `solve_gap/gap_residual` 与 legacy adapters 组装为最小 `ImplicitProblem`。
"""

function build_pnjl_fixedmu_problem(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    adapters = build_pnjl_fixedmu_adapters(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        kwargs...,
    )
    return ImplicitProblem(
        forward_solve=adapters.forward_solve,
        conditions=adapters.conditions,
        x_dim=gap_state_dim(model),
        theta_dim=2,
    )
end

function build_pnjl_flavor_mu_problem(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    adapters = build_pnjl_flavor_mu_adapters(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        kwargs...,
    )
    return ImplicitProblem(
        forward_solve=adapters.forward_solve,
        conditions=adapters.conditions,
        x_dim=gap_state_dim(model),
        theta_dim=4,
    )
end

function build_njl_problem(
    model::AbstractNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    solver::AbstractGapSolver=NLsolveGapSolver(),
    kwargs...
)
    dim = gap_state_dim(model)
    dim == 2 || dim == 3 || throw(ArgumentError("build_njl_problem supports dim=2/3, got $(dim)"))

    forward_solve = function (θ::AbstractVector)
        T_fm = Float64(θ[1])
        μ_fm = Float64(θ[2])
        mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
        st = solve_gap(model, T_fm, mu_vec; solver=solver, xi=xi, p_num=p_num, t_num=t_num, kwargs...)
        if dim == 2
            return ([st.phi[1], st.phi[2]], nothing)
        end
        return (collect(st.phi), nothing)
    end

    conditions = function (θ::AbstractVector, x::AbstractVector, z)
        _ = z
        return _implicit_conditions_with_shared_core(model, θ, x;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            kwargs...)
    end

    return ImplicitProblem(
        forward_solve=forward_solve,
        conditions=conditions,
        x_dim=dim,
        theta_dim=2,
    )
end
