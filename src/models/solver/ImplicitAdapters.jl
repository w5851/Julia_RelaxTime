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

    return ImplicitProblem(
        forward_solve=forward_solve,
        conditions=conditions,
        x_dim=dim,
        theta_dim=2,
    )
end
