"""Public RMF/RMFT point and scan workflows for GasLiquid."""
module GasLiquidWorkflow

using SHA: sha256
using StaticArrays
import ...Models
using ..GasLiquidEquationSet: GasLiquidCoreParams, GasLiquidState, RMFSolverResult
using ..GasLiquidEquationSet: solve_equilibrium, state_vector
using ..GasLiquidEquationSet: parameter_summary
using ..GasLiquidThermodynamics: pressure_density_entropy_energy

export solve_gas_liquid_point
export solve_gas_liquid_rmf_point
export run_gas_liquid_tmu_scan, run_gas_liquid_trho_scan
export build_gas_liquid_result_row, build_gas_liquid_manifest
export GAS_LIQUID_RMF_SCHEMA_VERSION, GAS_LIQUID_MANIFEST_SCHEMA_VERSION

const GAS_LIQUID_RMF_SCHEMA_VERSION = "gas_liquid_rmf_row_v1"
const GAS_LIQUID_MANIFEST_SCHEMA_VERSION = "gas_liquid_rmf_manifest_v1"
const DEFAULT_SOURCE_MATRIX_IDS = (
    "thesis_2.2_Eq2.40-2.70",
    "DiToro_2006_Eq1-13_Eq21-22",
    "QHDI_nucl-th_9311005_Eq1-4",
)

@inline _create_model(args...; kwargs...) = getfield(Models, :create_model)(args...; kwargs...)

@inline function _as_mev(value_fm::Real, p::GasLiquidCoreParams)
    return Float64(value_fm) * p.hbarc_MeV_fm
end

@inline function _as_fm(value_MeV::Real, p::GasLiquidCoreParams)
    return Float64(value_MeV) / p.hbarc_MeV_fm
end

function _config_hash(p::GasLiquidCoreParams)
    bytes = sha256(string(parameter_summary(p)))
    return join(lowercase(string(b, base=16, pad=2)) for b in bytes)
end

function _default_git_sha()
    configured = strip(get(ENV, "GAS_LIQUID_GIT_SHA", ""))
    !isempty(configured) && return configured
    root = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
    try
        return strip(readchomp(`git -C $root rev-parse HEAD`))
    catch
        return "unknown"
    end
end

function _solver_row(result::RMFSolverResult, thermo, p::GasLiquidCoreParams;
    T_MeV::Real,
    profile::AbstractString,
    run_id::AbstractString,
    point_id::AbstractString,
    formal_status::Symbol=:diagnostic_only,
    source_matrix_ids=DEFAULT_SOURCE_MATRIX_IDS)
    formal_status === :diagnostic_only || throw(ArgumentError("RMF workflow only emits :diagnostic_only; formal promotion requires the documented external gate"))
    mu_B_MeV = _as_mev((result.mu_p + result.mu_n) / 2, p)
    mu_3_MeV = _as_mev((result.mu_p - result.mu_n) / 2, p)
    return (
        schema_version=GAS_LIQUID_RMF_SCHEMA_VERSION,
        run_id=String(run_id),
        point_id=String(point_id),
        mode=result.mode,
        profile=String(profile),
        T_MeV=Float64(T_MeV),
        mu_B_MeV=mu_B_MeV,
        mu_3_MeV=mu_3_MeV,
        mu_p_MeV=_as_mev(result.mu_p, p),
        mu_n_MeV=_as_mev(result.mu_n, p),
        mu_p_star_MeV=_as_mev(thermo.mu_p_star, p),
        mu_n_star_MeV=_as_mev(thermo.mu_n_star, p),
        rho_B_fm3=thermo.rho_B,
        rho_3_fm3=thermo.rho_3,
        rho_p_fm3=thermo.rho_p,
        rho_n_fm3=thermo.rho_n,
        rho_s_fm3=thermo.rho_s_p + thermo.rho_s_n,
        rho_s_p_fm3=thermo.rho_s_p,
        rho_s_n_fm3=thermo.rho_s_n,
        rho_s3_fm3=thermo.rho_s3,
        M_p_MeV=_as_mev(thermo.masses[1], p),
        M_n_MeV=_as_mev(thermo.masses[2], p),
        S_inv_fm=thermo.S,
        D_inv_fm=thermo.D,
        W_inv_fm=thermo.W,
        R_inv_fm=thermo.R,
        omega_fm4=thermo.omega,
        pressure_fm4=thermo.pressure,
        entropy_fm3=thermo.entropy,
        energy_fm4=thermo.energy,
        converged=result.converged,
        solver_status=result.solver_status,
        iterations=result.iterations,
        residual_norm=result.residual_norm,
        failure_reason=result.failure_reason,
        attempts=result.attempts,
        fallback_used=result.fallback_used,
        quadrature_settings=thermo.quadrature,
        source_matrix_ids=Tuple(String.(source_matrix_ids)),
        formal_status=formal_status,
    )
end

function _solve_core_point(T_MeV::Real, p::GasLiquidCoreParams;
    mode::Symbol=:fixed_mu,
    mu_B_MeV::Real=0.0,
    mu_3_MeV::Real=0.0,
    rho_B_fm3=nothing,
    rho_3_fm3=nothing,
    alpha=nothing,
    initial_guess=nothing,
    p_num::Int=96,
    p_max_inv_fm::Real=10.0,
    max_iter::Int=200)
    T_fm = _as_fm(T_MeV, p)
    if mode === :fixed_mu
        result = solve_equilibrium(T_fm, 0.0; params=p, mode=:fixed_mu,
            mu_B=_as_fm(mu_B_MeV, p), mu_3=_as_fm(mu_3_MeV, p),
            initial_guess=initial_guess, p_num=p_num, p_max_inv_fm=p_max_inv_fm, max_iter=max_iter,
            return_result=true)
    elseif mode === :fixed_rho
        rho_B_fm3 === nothing && throw(ArgumentError("rho_B_fm3 is required for :fixed_rho"))
        result = solve_equilibrium(T_fm, 0.0; params=p, mode=:fixed_rho,
            rho_B_target=Float64(rho_B_fm3), rho_3_target=rho_3_fm3 === nothing ? nothing : Float64(rho_3_fm3),
            alpha=alpha, initial_guess=initial_guess, p_num=p_num, p_max_inv_fm=p_max_inv_fm, max_iter=max_iter,
            return_result=true)
    else
        throw(ArgumentError("mode must be :fixed_mu or :fixed_rho"))
    end
    thermo = pressure_density_entropy_energy(result.state, T_fm, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm,
        mu_p_target=result.mu_p, mu_n_target=result.mu_n)
    return result, thermo
end

"""Versioned MeV-facing RMF point workflow.

All rows are diagnostic-only by default. `mode=:fixed_rho` accepts either
`rho_3_fm3` or `alpha=(rho_n-rho_p)/rho_B`.
"""
function solve_gas_liquid_rmf_point(; T_MeV::Real,
    mode::Symbol=:fixed_mu,
    mu_B_MeV::Real=0.0,
    mu_3_MeV::Real=0.0,
    rho_B_fm3=nothing,
    rho_3_fm3=nothing,
    alpha=nothing,
    profile::AbstractString="default",
    model=nothing,
    initial_guess=nothing,
    p_num::Int=96,
    p_max_inv_fm::Real=10.0,
    max_iter::Int=200,
    run_id::AbstractString="diagnostic",
    point_id::AbstractString="point-1",
    formal_status::Symbol=:diagnostic_only,
    source_matrix_ids=DEFAULT_SOURCE_MATRIX_IDS)
    rmf_model = model === nothing ? _create_model(:GasLiquid; profile=String(profile)) : model
    p = rmf_model.params
    result, thermo = _solve_core_point(T_MeV, p; mode=mode, mu_B_MeV=mu_B_MeV, mu_3_MeV=mu_3_MeV,
        rho_B_fm3=rho_B_fm3, rho_3_fm3=rho_3_fm3, alpha=alpha, initial_guess=initial_guess,
        p_num=p_num, p_max_inv_fm=p_max_inv_fm, max_iter=max_iter)
    row = _solver_row(result, thermo, p; T_MeV=T_MeV, profile=profile, run_id=run_id, point_id=point_id,
        formal_status=formal_status, source_matrix_ids=source_matrix_ids)
    return merge(row, (state=result.state, thermo=thermo, model=rmf_model, result=result, config_hash=_config_hash(p)))
end

"""Compatibility workflow retaining the historical `(T_fm,mu_fm)` signature."""
function solve_gas_liquid_point(T_fm::Real, mu_fm::Real; model_kind::Symbol=:GasLiquid, profile::AbstractString="default", kwargs...)
    model = _create_model(model_kind; profile=String(profile))
    p = model.params
    result = solve_equilibrium(T_fm, mu_fm; params=p, mode=:fixed_mu,
        initial_guess=get(kwargs, :initial_guess, nothing), p_num=Int(get(kwargs, :p_num, 96)),
        p_max_inv_fm=Float64(get(kwargs, :p_max_inv_fm, 10.0)), max_iter=Int(get(kwargs, :max_iter, 200)), return_result=true)
    thermo = pressure_density_entropy_energy(result.state, T_fm, p; p_num=Int(get(kwargs, :p_num, 96)),
        p_max_inv_fm=Float64(get(kwargs, :p_max_inv_fm, 10.0)), mu_p_target=result.mu_p, mu_n_target=result.mu_n)
    phi = SVector{3, Float64}(result.state.S, result.state.D, 0.0)
    x_state = Models.MeanFieldState(phi; Phi=1.0, PhiBar=1.0)
    return (
        model_kind=model_kind,
        T_fm=Float64(T_fm),
        mu_fm=Float64(mu_fm),
        x_state=x_state,
        omega=thermo.omega,
        pressure=thermo.pressure,
        rho=thermo.rho_B,
        entropy=thermo.entropy,
        energy=thermo.energy,
        converged=result.converged,
        solver_status=result.solver_status,
        iterations=result.iterations,
        residual_norm=result.residual_norm,
        rho_p=thermo.rho_p,
        rho_n=thermo.rho_n,
        rho_3=thermo.rho_3,
        S_inv_fm=result.state.S,
        D_inv_fm=result.state.D,
        W_inv_fm=result.fields.W,
        R_inv_fm=result.fields.R,
        formal_status=:diagnostic_only,
    )
end

function build_gas_liquid_result_row(result;
    T_MeV=get(result, :T_MeV, nothing),
    profile=get(result, :profile, "default"),
    run_id=get(result, :run_id, "diagnostic"),
    point_id=get(result, :point_id, "point-1"),
    formal_status::Symbol=:diagnostic_only,
    source_matrix_ids=get(result, :source_matrix_ids, DEFAULT_SOURCE_MATRIX_IDS))
    required = (:thermo, :model, :state)
    all(haskey(result, key) for key in required) || throw(ArgumentError("result must be a RMF point result"))
    T_MeV === nothing && throw(ArgumentError("T_MeV is required to build a result row"))
    return _solver_row(result.result, result.thermo, result.model.params;
        T_MeV=T_MeV, profile=profile, run_id=run_id, point_id=point_id,
        formal_status=formal_status, source_matrix_ids=source_matrix_ids)
end

function build_gas_liquid_manifest(rows::AbstractVector;
    profile::AbstractString="default",
    run_id::AbstractString="diagnostic",
    git_sha::AbstractString=_default_git_sha(),
    effective_config_hash::AbstractString="unknown",
    source_matrix_ids=DEFAULT_SOURCE_MATRIX_IDS,
    grid=nothing,
    solver_settings=(;),
    quadrature_settings=(;),
    output_hashes=(),
    formal_status::Symbol=:diagnostic_only)
    formal_status === :diagnostic_only || throw(ArgumentError("RMF manifest only emits :diagnostic_only; formal promotion requires the documented external gate"))
    failed = [row.point_id for row in rows if !row.converged]
    return (
        schema_version=GAS_LIQUID_MANIFEST_SCHEMA_VERSION,
        run_id=String(run_id),
        git_sha=String(git_sha),
        effective_config_hash=String(effective_config_hash),
        profile=String(profile),
        source_matrix_ids=Tuple(String.(source_matrix_ids)),
        grid=grid,
        solver_settings=solver_settings,
        quadrature_settings=quadrature_settings,
        row_count=length(rows),
        failed_points=Tuple(String.(failed)),
        output_hashes=Tuple(String.(output_hashes)),
        formal_status=formal_status,
    )
end

function run_gas_liquid_tmu_scan(T_MeV_grid, mu_B_MeV_grid;
    mu_3_MeV::Real=0.0,
    profile::AbstractString="default",
    model=nothing,
    p_num::Int=96,
    p_max_inv_fm::Real=10.0,
    max_iter::Int=200,
    run_id::AbstractString="diagnostic-tmu",
    formal_status::Symbol=:diagnostic_only)
    rmf_model = model === nothing ? _create_model(:GasLiquid; profile=String(profile)) : model
    rows = Any[]
    seed = nothing
    for (i, T_MeV) in enumerate(T_MeV_grid), (j, mu_B_MeV) in enumerate(mu_B_MeV_grid)
        point = solve_gas_liquid_rmf_point(; T_MeV=T_MeV, mu_B_MeV=mu_B_MeV, mu_3_MeV=mu_3_MeV,
            profile=profile, model=rmf_model, initial_guess=seed, p_num=p_num, p_max_inv_fm=p_max_inv_fm,
            max_iter=max_iter, run_id=run_id, point_id="tmu-$(i)-$(j)", formal_status=formal_status)
        row = build_gas_liquid_result_row(point)
        push!(rows, row)
        point.converged && (seed = collect(state_vector(point.state)))
    end
    manifest = build_gas_liquid_manifest(rows; profile=profile, run_id=run_id, effective_config_hash=_config_hash(rmf_model.params),
        grid=(T_MeV=collect(T_MeV_grid), mu_B_MeV=collect(mu_B_MeV_grid), mu_3_MeV=Float64(mu_3_MeV)),
        solver_settings=(max_iter=max_iter,), quadrature_settings=(p_num=p_num, p_max_inv_fm=Float64(p_max_inv_fm)),
        formal_status=formal_status)
    return (schema_version=GAS_LIQUID_RMF_SCHEMA_VERSION, rows=rows, manifest=manifest, formal_status=formal_status)
end

function run_gas_liquid_trho_scan(T_MeV_grid, rho_B_fm3_grid;
    alpha::Real=0.0,
    profile::AbstractString="default",
    model=nothing,
    p_num::Int=96,
    p_max_inv_fm::Real=10.0,
    max_iter::Int=200,
    run_id::AbstractString="diagnostic-trho",
    formal_status::Symbol=:diagnostic_only)
    rmf_model = model === nothing ? _create_model(:GasLiquid; profile=String(profile)) : model
    rows = Any[]
    seed = nothing
    for (i, T_MeV) in enumerate(T_MeV_grid), (j, rho_B_fm3) in enumerate(rho_B_fm3_grid)
        point = solve_gas_liquid_rmf_point(; T_MeV=T_MeV, mode=:fixed_rho, rho_B_fm3=rho_B_fm3, alpha=alpha,
            profile=profile, model=rmf_model, initial_guess=seed, p_num=p_num, p_max_inv_fm=p_max_inv_fm,
            max_iter=max_iter, run_id=run_id, point_id="trho-$(i)-$(j)", formal_status=formal_status)
        row = build_gas_liquid_result_row(point)
        push!(rows, row)
        point.converged && (seed = collect(state_vector(point.state)))
    end
    manifest = build_gas_liquid_manifest(rows; profile=profile, run_id=run_id, effective_config_hash=_config_hash(rmf_model.params),
        grid=(T_MeV=collect(T_MeV_grid), rho_B_fm3=collect(rho_B_fm3_grid), alpha=Float64(alpha)),
        solver_settings=(max_iter=max_iter,), quadrature_settings=(p_num=p_num, p_max_inv_fm=Float64(p_max_inv_fm)),
        formal_status=formal_status)
    return (schema_version=GAS_LIQUID_RMF_SCHEMA_VERSION, rows=rows, manifest=manifest, formal_status=formal_status)
end

end # module
