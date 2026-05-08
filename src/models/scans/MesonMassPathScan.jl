"""
    MesonMassPathScan

沿已定义路径逐点运行 meson-mass workflow 的正式入口。

当前最小支持两类路径：

- freeze-out: 直接消费 `(T, μ_B)` 路径点列
- isentropic: 先用 `FixedSigma` 求出路径点，再复用 meson workflow

第一版只固化 meson-mass path scan，不把所有 path family 抽成共享黑箱框架。
"""
module MesonMassPathScan

using Main.Constants_PNJL: ħc_MeV_fm

using ..Models: FixedSigma, create_model, solve_constraint
using ..MesonMassWorkflow
using ..FreezeoutProfiles
using ..FreezeoutPathProfiles
using ..IsentropicPathProfiles
using ..ScanCommon

export run_meson_mass_path_scan
export run_freezeout_meson_mass_scan
export run_isentropic_meson_mass_scan
export DEFAULT_MESON_MASS_PATH_OUTPUT_PATH
export DEFAULT_FREEZEOUT_MESON_MASS_OUTPUT_PATH
export DEFAULT_ISENTROPIC_MESON_MASS_OUTPUT_PATH
export HEADER

const DEFAULT_MESON_MASS_PATH_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "meson_mass", "path_scan",
    "meson_mass_path_scan.csv",
))

const DEFAULT_FREEZEOUT_MESON_MASS_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "meson_mass", "path_scan", "freezeout",
    "freezeout_meson_mass_path_scan.csv",
))

const DEFAULT_ISENTROPIC_MESON_MASS_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "meson_mass", "path_scan", "isentropic",
    "isentropic_meson_mass_path_scan.csv",
))

@inline _hasprop(x, name::Symbol) = x isa NamedTuple ? haskey(x, name) : hasproperty(x, name)
@inline _prop(x, name::Symbol) = x isa NamedTuple ? getfield(x, name) : getproperty(x, name)
@inline _get(x, name::Symbol, default) = _hasprop(x, name) ? _prop(x, name) : default

@inline function _required_real(point, field::Symbol)
    _hasprop(point, field) || throw(ArgumentError("meson mass path point missing required field: $(field)"))
    value = Float64(_prop(point, field))
    isfinite(value) || throw(ArgumentError("meson mass path point field $(field) must be finite, got $(value)"))
    return value
end

@inline function _optional_real(point, field::Symbol, fallback::Float64=NaN)
    _hasprop(point, field) || return fallback
    value = Float64(_prop(point, field))
    isfinite(value) || throw(ArgumentError("meson mass path point field $(field) must be finite, got $(value)"))
    return value
end

@inline function _optional_string(point, field::Symbol, fallback::String="")
    _hasprop(point, field) || return fallback
    value = _prop(point, field)
    value isa AbstractString || throw(ArgumentError("meson mass path point field $(field) must be string, got $(typeof(value))"))
    return String(value)
end

function _normalize_path_point(point, idx::Int)
    T_MeV = _required_real(point, :T_MeV)
    muB_MeV = _optional_real(point, :muB_MeV, NaN)
    sigma_target = _optional_real(point, :sigma_target, NaN)
    if !isfinite(muB_MeV)
        isfinite(sigma_target) || throw(ArgumentError("path point requires either muB_MeV or sigma_target"))
    end

    family = _optional_string(point, :path_family, isfinite(muB_MeV) ? "freezeout" : "isentropic")
    profile = _optional_string(point, :path_profile, "default")
    segment = _optional_string(point, :path_segment, family == "isentropic" ? "fixed_sigma" : "path")
    label = _optional_string(point, :path_label, "$(family):$(profile)")
    sqrt_s_NN_GeV = _optional_real(point, :sqrt_s_NN_GeV, NaN)
    freezeout_profile = _optional_string(point, :freezeout_profile, "")

    return (
        path_family=family,
        path_profile=profile,
        path_segment=segment,
        path_point_index=Int(round(Float64(_get(point, :path_point_index, idx - 1)))),
        path_order_key=Float64(_get(point, :path_order_key, T_MeV)),
        path_label=label,
        freezeout_profile=freezeout_profile,
        sigma_target=sigma_target,
        sqrt_s_NN_GeV=sqrt_s_NN_GeV,
        T_MeV=T_MeV,
        muB_MeV=muB_MeV,
    )
end

function _header_cols(mesons::Tuple{Vararg{Symbol}})
    cols = String[
        "path_family",
        "path_profile",
        "path_segment",
        "path_point_index",
        "path_order_key",
        "path_label",
        "freezeout_profile",
        "sigma_target",
        "sqrt_s_NN_GeV",
        "T_MeV",
        "muB_MeV",
        "muq_MeV",
        "xi",
        "equilibrium_converged",
        "equilibrium_iterations",
        "equilibrium_residual_norm",
        "Phi",
        "Phibar",
        "m_u",
        "m_d",
        "m_s",
    ]

    for m in mesons
        push!(cols, "M_" * String(m))
        push!(cols, "Gamma_" * String(m))
        push!(cols, "converged_" * String(m))
        push!(cols, "residual_" * String(m))

        if m in (:eta, :eta_prime, :sigma, :sigma_prime)
            push!(cols, "thr_uu_" * String(m))
            push!(cols, "thr_ss_" * String(m))
            push!(cols, "thr_min_" * String(m))
            push!(cols, "gap_uu_" * String(m))
            push!(cols, "gap_ss_" * String(m))
            push!(cols, "gap_min_" * String(m))
        else
            push!(cols, "threshold_" * String(m))
            push!(cols, "gap_" * String(m))
        end
    end

    push!(cols, "message")
    return cols
end

const HEADER = join(_header_cols(MesonMassWorkflow.DEFAULT_MESONS), ',')

@inline function _path_group_key(pt)
    return (pt.path_family, pt.path_profile, pt.path_label)
end

@inline function _schedule_points(points)
    sort!(points; by=pt -> (pt.path_family, pt.path_profile, pt.path_label, pt.path_order_key, pt.path_point_index))
    return points
end

@inline function _equilibrium_diag(eq_result)
    if eq_result === nothing
        return (converged=false, iterations=-1, residual_norm=NaN)
    end
    converged = hasproperty(eq_result, :converged) ? Bool(getproperty(eq_result, :converged)) : false
    iterations = hasproperty(eq_result, :iterations) ? Int(getproperty(eq_result, :iterations)) : -1
    residual_norm = hasproperty(eq_result, :residual_norm) ? Float64(getproperty(eq_result, :residual_norm)) : NaN
    return (converged=converged, iterations=iterations, residual_norm=residual_norm)
end

@inline function _constraint_good(result)
    return hasproperty(result, :converged) &&
           Bool(getproperty(result, :converged)) &&
           hasproperty(result, :residual_norm) &&
           isfinite(Float64(getproperty(result, :residual_norm))) &&
           Float64(getproperty(result, :residual_norm)) <= 1e-6
end

function _meson_row_dict(pt, xi::Float64, muB_MeV::Float64, eq_diag, result)
    qp = result.quark_params
    tp = result.thermo_params
    row = Dict{String,Any}(
        "path_family" => pt.path_family,
        "path_profile" => pt.path_profile,
        "path_segment" => pt.path_segment,
        "path_point_index" => pt.path_point_index,
        "path_order_key" => pt.path_order_key,
        "path_label" => pt.path_label,
        "freezeout_profile" => pt.freezeout_profile,
        "sigma_target" => pt.sigma_target,
        "sqrt_s_NN_GeV" => pt.sqrt_s_NN_GeV,
        "T_MeV" => pt.T_MeV,
        "muB_MeV" => muB_MeV,
        "muq_MeV" => muB_MeV / 3.0,
        "xi" => xi,
        "equilibrium_converged" => eq_diag.converged,
        "equilibrium_iterations" => eq_diag.iterations,
        "equilibrium_residual_norm" => eq_diag.residual_norm,
        "Phi" => tp.Φ,
        "Phibar" => tp.Φbar,
        "m_u" => qp.m.u,
        "m_d" => qp.m.d,
        "m_s" => qp.m.s,
    )

    for meson in keys(result.meson_results)
        mr = result.meson_results[meson]
        row["M_" * String(meson)] = mr.mass
        row["Gamma_" * String(meson)] = mr.gamma
        row["converged_" * String(meson)] = mr.converged
        row["residual_" * String(meson)] = mr.residual

        if meson in (:eta, :eta_prime, :sigma, :sigma_prime)
            thr = mr.threshold
            gaps = mr.gaps
            row["thr_uu_" * String(meson)] = thr.uu
            row["thr_ss_" * String(meson)] = thr.ss
            row["thr_min_" * String(meson)] = thr.min
            row["gap_uu_" * String(meson)] = gaps.uu
            row["gap_ss_" * String(meson)] = gaps.ss
            row["gap_min_" * String(meson)] = gaps.min
        else
            row["threshold_" * String(meson)] = mr.threshold
            row["gap_" * String(meson)] = mr.gap
        end
    end

    return row
end

@inline function _format_cell(x)
    x isa Bool && return x ? "true" : "false"
    x isa Real && return isfinite(Float64(x)) ? ScanCommon.fmt(x) : "NaN"
    return string(x)
end

function _row_to_values(cols::Vector{String}, row::Dict{String,Any})
    return map(cols) do col
        value = get(row, col, "")
        col == "message" ? ScanCommon.quote_csv(ScanCommon.clean_message(String(value))) : _format_cell(value)
    end
end

function _write_failure_row(io, cols::Vector{String}, pt, xi::Float64, muB_MeV::Float64, eq_diag, message::AbstractString)
    row = Dict{String,Any}(
        "path_family" => pt.path_family,
        "path_profile" => pt.path_profile,
        "path_segment" => pt.path_segment,
        "path_point_index" => pt.path_point_index,
        "path_order_key" => pt.path_order_key,
        "path_label" => pt.path_label,
        "freezeout_profile" => pt.freezeout_profile,
        "sigma_target" => pt.sigma_target,
        "sqrt_s_NN_GeV" => pt.sqrt_s_NN_GeV,
        "T_MeV" => pt.T_MeV,
        "muB_MeV" => muB_MeV,
        "muq_MeV" => isfinite(muB_MeV) ? muB_MeV / 3.0 : NaN,
        "xi" => xi,
        "equilibrium_converged" => eq_diag.converged,
        "equilibrium_iterations" => eq_diag.iterations,
        "equilibrium_residual_norm" => eq_diag.residual_norm,
        "message" => message,
    )
    println(io, join(_row_to_values(cols, row), ','))
end

function _write_success_row(io, cols::Vector{String}, pt, xi::Float64, muB_MeV::Float64, eq_diag, result)
    row = _meson_row_dict(pt, xi, muB_MeV, eq_diag, result)
    row["message"] = ""
    println(io, join(_row_to_values(cols, row), ','))
end

function _solve_isentropic_equilibrium(model, pt, xi::Float64, seed_guess, mu0;
    p_num::Int,
    t_num::Int,
    iterations::Int,
)
    T_fm = pt.T_MeV / ħc_MeV_fm
    kwargs = (
        seed_guess=seed_guess,
        mu0=mu0,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        iterations=iterations,
        residual_norm_max=1e-6,
    )
    return solve_constraint(
        model,
        FixedSigma(pt.sigma_target),
        T_fm;
        kwargs...,
    )
end

function _decorate_freezeout_points(points, freezeout_profile_name::String)
    label = "freezeout:" * freezeout_profile_name
    return map(points) do pt
        merge(pt, (
            path_family="freezeout",
            path_label=label,
            freezeout_profile=freezeout_profile_name,
        ))
    end
end

function run_meson_mass_path_scan(;
    points,
    xi::Float64=0.0,
    mesons::Tuple{Vararg{Symbol}}=MesonMassWorkflow.DEFAULT_MESONS,
    output_path::AbstractString=DEFAULT_MESON_MASS_PATH_OUTPUT_PATH,
    overwrite::Bool=true,
    p_num::Int=24,
    t_num::Int=8,
    max_iter::Int=40,
    solver_backend::Symbol=:auto,
    solver_kwargs::NamedTuple=(;),
    mass_kwargs::NamedTuple=(;),
)
    raw_points = collect(points)
    isempty(raw_points) && throw(ArgumentError("points must not be empty"))
    normalized = _schedule_points([_normalize_path_point(pt, idx) for (idx, pt) in enumerate(raw_points)])
    cols = _header_cols(mesons)

    mkpath(dirname(output_path))
    !overwrite && isfile(output_path) && throw(ArgumentError("output_path already exists: $(output_path)"))

    model = create_model(:PNJL)
    effective_solver_kwargs = merge((iterations=max_iter,), solver_kwargs)
    effective_mass_kwargs = merge((iterations=max_iter,), mass_kwargs)

    continuation_state = nothing
    previous_group = nothing
    sigma_seed_state = nothing
    sigma_mu0 = 0.0
    open(output_path, "w") do io
        println(io, join(cols, ','))
        for pt in normalized
            group = _path_group_key(pt)
            if previous_group !== nothing && group != previous_group
                continuation_state = nothing
                sigma_seed_state = nothing
                sigma_mu0 = 0.0
            end
            previous_group = group

            eq_diag = (converged=false, iterations=-1, residual_norm=NaN)
            muB_MeV = pt.muB_MeV

            try
                result = if isfinite(pt.muB_MeV)
                    muq_fm = (pt.muB_MeV / 3.0) / ħc_MeV_fm
                    res = MesonMassWorkflow.solve_gap_and_meson_point(
                        pt.T_MeV / ħc_MeV_fm,
                        muq_fm;
                        xi=xi,
                        solver_backend=solver_backend,
                        mesons=mesons,
                        continuation_state=continuation_state,
                        p_num=p_num,
                        t_num=t_num,
                        solver_kwargs=effective_solver_kwargs,
                        mass_kwargs=effective_mass_kwargs,
                    )
                    eq_diag = _equilibrium_diag(res.equilibrium)
                    res
                else
                    sigma_result = _solve_isentropic_equilibrium(
                        model,
                        pt,
                        xi,
                        sigma_seed_state,
                        sigma_mu0;
                        p_num=p_num,
                        t_num=t_num,
                        iterations=max_iter,
                    )
                    eq_diag = _equilibrium_diag(sigma_result)
                    _constraint_good(sigma_result) || throw(ArgumentError("FixedSigma solve did not satisfy convergence contract"))
                    sigma_seed_state = collect(sigma_result.x_state)
                    sigma_mu0 = Float64(sigma_result.mu_vec[1])
                    muB_MeV = 3.0 * sigma_mu0 * ħc_MeV_fm
                    MesonMassWorkflow.solve_gap_and_meson_point(
                        pt.T_MeV / ħc_MeV_fm,
                        sigma_mu0;
                        xi=xi,
                        solver_backend=solver_backend,
                        mesons=mesons,
                        seed_state=sigma_seed_state,
                        continuation_state=continuation_state,
                        p_num=p_num,
                        t_num=t_num,
                        solver_kwargs=effective_solver_kwargs,
                        mass_kwargs=effective_mass_kwargs,
                        flavor_mu_override=collect(sigma_result.mu_vec),
                    )
                end

                continuation_state = result.continuation_state
                _write_success_row(io, cols, pt, xi, muB_MeV, eq_diag, result)
            catch err
                _write_failure_row(io, cols, pt, xi, muB_MeV, eq_diag, sprint(showerror, err))
            end
        end
    end

    return (
        output_path=output_path,
        points=length(normalized),
        xi=xi,
        workflow_entry="Models.run_meson_mass_path_scan",
    )
end

function run_freezeout_meson_mass_scan(;
    sqrt_s_NN_values,
    xi::Float64=0.0,
    freezeout_profile_name::AbstractString="default",
    path_profile_name::AbstractString="baseline_freezeout",
    output_path::AbstractString=DEFAULT_FREEZEOUT_MESON_MASS_OUTPUT_PATH,
    traversal::Symbol=:sqrts_ascending,
    kwargs...,
)
    freezeout_profile = FreezeoutProfiles.load_freezeout_profile(profile=String(freezeout_profile_name))
    path_profile = FreezeoutPathProfiles.load_freezeout_path_profile(profile=String(path_profile_name))
    points = FreezeoutPathProfiles.build_freezeout_path_points(
        sqrt_s_NN_values;
        freezeout_profile=freezeout_profile,
        path_profile=path_profile,
        traversal=traversal,
    )
    decorated = _decorate_freezeout_points(points, freezeout_profile.profile_name)
    result = run_meson_mass_path_scan(;
        points=decorated,
        xi=xi,
        output_path=output_path,
        kwargs...,
    )
    return merge(result, (
        freezeout_profile=freezeout_profile.profile_name,
        path_profile=path_profile.profile_name,
        workflow_entry="Models.run_freezeout_meson_mass_scan",
    ))
end

function run_isentropic_meson_mass_scan(;
    T_MeV_values,
    sigma_target::Real,
    xi::Float64=0.0,
    path_profile_name::AbstractString="default",
    output_path::AbstractString=DEFAULT_ISENTROPIC_MESON_MASS_OUTPUT_PATH,
    traversal::Symbol=:T_ascending,
    kwargs...,
)
    profile = IsentropicPathProfiles.load_isentropic_path_profile(
        sigma_target=sigma_target,
        profile=String(path_profile_name),
    )
    points = IsentropicPathProfiles.build_isentropic_path_points(
        T_MeV_values;
        sigma_target=sigma_target,
        profile=profile,
        traversal=traversal,
    )
    result = run_meson_mass_path_scan(;
        points=points,
        xi=xi,
        output_path=output_path,
        kwargs...,
    )
    return merge(result, (
        sigma_target=profile.sigma_target,
        path_profile=profile.profile_name,
        workflow_entry="Models.run_isentropic_meson_mass_scan",
    ))
end

end # module MesonMassPathScan
