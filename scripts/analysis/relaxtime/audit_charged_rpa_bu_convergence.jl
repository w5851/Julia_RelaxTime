"""
Third follow-up diagnostic: node/cutoff, Bose-support and four-algorithm
comparison on the existing baseline chemical freeze-out parameterization.

This script is deliberately analysis-only. It consumes the stable Models
meson workflow once per point, evaluates the four existing density algorithms
on that state, and writes local CSV/README artifacts. No production default is
changed and failed algorithms remain explicit rows.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using Printf: @printf
using .Constants_PNJL: ħc_MeV_fm
using .Models
using .Models: solve_gap_and_meson_point,
               solve_meson_density_from_meson_point,
               solve_strict_bw_meson_density_from_meson_point,
               solve_phase_shift_meson_density_from_meson_point
using .Models.FreezeoutProfiles: load_freezeout_profile, build_freezeout_scan_points
using .Models.FlavorChemicalProfiles: load_flavor_chemical_profile, flavor_mu_profile_fm
using Main.RelaxTime.BUPhaseGates: bose_support_gate,
                                    convergence_gate,
                                    four_density_algorithm_labels

const DEFAULT_OUTPUT_DIR = joinpath(
    PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_bu_convergence",
)
const DEFAULT_SQRTS = [3.0, 7.7, 200.0]

@inline _env_bool(name::AbstractString, default::Bool) = lowercase(get(ENV, name, default ? "true" : "false")) in ("1", "true", "yes")
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))
@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _ratio(numerator::Real, denominator::Real) = denominator > 0.0 ? numerator / denominator : NaN
@inline _fmt(x) = x isa Bool ? (x ? "true" : "false") : string(x)

function _float_list(raw::AbstractString)
    values = Float64[]
    for item in split(raw, ',')
        value = strip(item)
        isempty(value) && continue
        push!(values, parse(Float64, value))
    end
    isempty(values) && throw(ArgumentError("MESON_CONVERGENCE_SQRTS must not be empty"))
    return sort(unique(values))
end

function _sqrts_values()
    raw = strip(get(ENV, "MESON_CONVERGENCE_SQRTS", ""))
    return isempty(raw) ? copy(DEFAULT_SQRTS) : _float_list(raw)
end

Base.@kwdef struct DensityNumerics
    qmax::Float64 = 4.0
    q_nodes::Int = 4
    omega_min::Float64 = 0.05
    omega_max::Float64 = 3.0
    omega_nodes::Int = 8
    eta::Float64 = 1.0e-3
    stable_q_nodes::Int = 32
    pole_iterations::Int = 8
    pole_residual_norm_max::Float64 = 1.0e-5
    phase_measure::Symbol = :single_charge_domega_over_pi
    phase_anchor::Symbol = :high_energy_zero
end

function _numerics(; refined::Bool=false)
    if refined
        return DensityNumerics(
            qmax=_env_float("MESON_CONVERGENCE_REFINED_QMAX", 6.0),
            q_nodes=_env_int("MESON_CONVERGENCE_REFINED_Q_NODES", 8),
            omega_min=_env_float("MESON_CONVERGENCE_REFINED_OMEGA_MIN", 0.05),
            omega_max=_env_float("MESON_CONVERGENCE_REFINED_OMEGA_MAX", 4.0),
            omega_nodes=_env_int("MESON_CONVERGENCE_REFINED_OMEGA_NODES", 16),
            eta=_env_float("MESON_CONVERGENCE_REFINED_ETA", 5.0e-4),
            stable_q_nodes=_env_int("MESON_CONVERGENCE_REFINED_STABLE_Q_NODES", 64),
            pole_iterations=_env_int("MESON_CONVERGENCE_REFINED_POLE_ITERATIONS", 12),
        )
    end
    return DensityNumerics(
        qmax=_env_float("MESON_CONVERGENCE_QMAX", 4.0),
        q_nodes=_env_int("MESON_CONVERGENCE_Q_NODES", 4),
        omega_min=_env_float("MESON_CONVERGENCE_OMEGA_MIN", 0.05),
        omega_max=_env_float("MESON_CONVERGENCE_OMEGA_MAX", 3.0),
        omega_nodes=_env_int("MESON_CONVERGENCE_OMEGA_NODES", 8),
        eta=_env_float("MESON_CONVERGENCE_ETA", 1.0e-3),
        stable_q_nodes=_env_int("MESON_CONVERGENCE_STABLE_Q_NODES", 32),
        pole_iterations=_env_int("MESON_CONVERGENCE_POLE_ITERATIONS", 8),
    )
end

function _validate_numerics(n::DensityNumerics)
    n.qmax > 0.0 && n.q_nodes > 1 && n.omega_max > n.omega_min && n.omega_nodes > 1 ||
        throw(ArgumentError("invalid convergence quadrature settings"))
    n.omega_min >= 0.0 && n.eta > 0.0 && n.stable_q_nodes > 1 ||
        throw(ArgumentError("invalid convergence density settings"))
    return n
end

function _chemical_potentials(flavor_profile, point)
    flavor_mu_profile_fm(flavor_profile, point.muq_MeV)
end

function _algorithm_density(
    algorithm::Symbol,
    charge::Symbol,
    meson_point,
    flavor_mu,
    n::DensityNumerics,
)
    plus = charge === :plus
    pi_channel = plus ? :pi_plus : :pi_minus
    k_channel = plus ? :K_plus : :K_minus
    mu_pi = plus ? flavor_mu.mu_u_fm - flavor_mu.mu_d_fm : flavor_mu.mu_d_fm - flavor_mu.mu_u_fm
    mu_K = plus ? flavor_mu.mu_u_fm - flavor_mu.mu_s_fm : flavor_mu.mu_s_fm - flavor_mu.mu_u_fm
    pi_mass = Float64(meson_point.meson_results[:pi].mass)
    k_mass = Float64(meson_point.meson_results[:K].mass)
    pi_support = bose_support_gate(pi_mass, mu_pi; omega_min=n.omega_min, omega_max=n.omega_max)
    k_support = bose_support_gate(k_mass, mu_K; omega_min=n.omega_min, omega_max=n.omega_max)
    # Stable/BW helpers consume the aggregate mass entries (:pi/:K); only the
    # phase-shift backend has charge-resolved polarization channels.
    density_pi_channel = algorithm === :phase_shift_bu ? pi_channel : :pi
    density_k_channel = algorithm === :phase_shift_bu ? k_channel : :K
    common = (
        pi_channel=density_pi_channel,
        k_channel=density_k_channel,
        μ_pi=mu_pi,
        μ_K=mu_K,
        d_pi=1,
        d_K=1,
    )

    result = if algorithm === :stable_particle_limit
        solve_meson_density_from_meson_point(
            meson_point;
            common...,
            num_q_nodes=n.stable_q_nodes,
        )
    elseif algorithm === :reduced_strict_bw
        solve_strict_bw_meson_density_from_meson_point(
            meson_point;
            common...,
            stage=:stage1_reduced,
            qmax=n.qmax,
            q_nodes=n.q_nodes,
            omega_min=n.omega_min,
            omega_max=n.omega_max,
            omega_nodes=n.omega_nodes,
        )
    elseif algorithm === :q_pole_strict_bw
        solve_strict_bw_meson_density_from_meson_point(
            meson_point;
            common...,
            stage=:stage2_qpole,
            qmax=n.qmax,
            q_nodes=n.q_nodes,
            omega_max=n.omega_max,
            omega_nodes=n.omega_nodes,
            solver_iterations=n.pole_iterations,
            pole_residual_norm_max=n.pole_residual_norm_max,
            pole_require_converged=false,
        )
    elseif algorithm === :phase_shift_bu
        solve_phase_shift_meson_density_from_meson_point(
            meson_point;
            common...,
            scheme=:gbu_reference,
            qmax=n.qmax,
            q_nodes=n.q_nodes,
            omega_min=n.omega_min,
            omega_max=n.omega_max,
            omega_nodes=n.omega_nodes,
            eta=n.eta,
            phase_anchor=n.phase_anchor,
            omega_measure=n.phase_measure,
            density_policy=:strict_normal_domain,
        )
    else
        throw(ArgumentError("unsupported density algorithm $(algorithm)"))
    end

    n_pi_value = Float64(result.n_pi)
    n_K_value = Float64(result.n_K)
    ratio_value = _ratio(n_K_value, n_pi_value)
    valid_density = isfinite(n_pi_value) && isfinite(n_K_value) && n_pi_value > 0.0 && n_K_value >= 0.0
    return (
        status=valid_density ? :ok : :invalid_density,
        message=valid_density ? "" : "density contract failed: require finite n_K>=0 and n_pi>0",
        n_pi=n_pi_value,
        n_K=n_K_value,
        ratio=ratio_value,
        bose_pi_status=pi_support.status,
        bose_K_status=k_support.status,
        bose_pi_min_E_minus_mu=pi_support.min_E_minus_mu,
        bose_K_min_E_minus_mu=k_support.min_E_minus_mu,
        n_pi_result=result,
    )
end

function _quark_conserved_densities(meson_point, n::DensityNumerics)
    eq = meson_point.equilibrium
    model = Models.get_cached_model(:PNJL)
    rho = Models.number_densities(
        model,
        eq.x_state,
        meson_point.thermo_params.T,
        eq.mu_vec;
        p_num=max(n.q_nodes, 4),
        t_num=4,
        xi=0.0,
    )
    rho_u = Float64(rho.quark[1] - rho.antiquark[1])
    rho_d = Float64(rho.quark[2] - rho.antiquark[2])
    rho_s = Float64(rho.quark[3] - rho.antiquark[3])
    return (
        rho_u_fm3=rho_u,
        rho_d_fm3=rho_d,
        rho_s_quark_fm3=rho_s,
        rho_B_fm3=(rho_u + rho_d + rho_s) / 3.0,
        rho_Q_fm3=(2.0 * rho_u - rho_d - rho_s) / 3.0,
        rho_strange_charge_fm3=-rho_s,
    )
end

function _safe_algorithm_density(algorithm, charge, meson_point, flavor_mu, n)
    try
        return _algorithm_density(algorithm, charge, meson_point, flavor_mu, n)
    catch err
        return (
            status=:failed,
            message=sprint(showerror, err),
            n_pi=NaN,
            n_K=NaN,
            ratio=NaN,
            bose_pi_status=:not_evaluated,
            bose_K_status=:not_evaluated,
            bose_pi_min_E_minus_mu=NaN,
            bose_K_min_E_minus_mu=NaN,
            n_pi_result=nothing,
        )
    end
end

function _solve_point(point, continuation_state, flavor_profile)
    flavor_mu = flavor_mu_profile_fm(flavor_profile, point.muq_MeV)
    flavor_override = flavor_profile.apply_to_equilibrium ?
        (flavor_mu.mu_u_fm, flavor_mu.mu_d_fm, flavor_mu.mu_s_fm) : nothing
    kwargs = (
        xi=0.0,
        mesons=(:pi, :K),
        continuation_state=continuation_state,
        mixed_branch_align=:strict_sign_binding,
        p_num=_env_int("MESON_CONVERGENCE_P_NUM", 8),
        t_num=_env_int("MESON_CONVERGENCE_T_NUM", 4),
        solver_kwargs=(; iterations=_env_int("MESON_CONVERGENCE_SOLVER_ITERATIONS", 30)),
        mass_kwargs=(; iterations=20),
        flavor_mu_override=flavor_override,
    )
    try
        return solve_gap_and_meson_point(point.T_fm, point.muB_fm / 3.0; kwargs...), :continuation
    catch first_error
        try
            fallback_kwargs = merge(kwargs, (continuation_state=nothing,))
            return solve_gap_and_meson_point(point.T_fm, point.muB_fm / 3.0; fallback_kwargs...), :fallback_no_continuation
        catch second_error
            return (failure=true, message=string(sprint(showerror, first_error), " | fallback: ", sprint(showerror, second_error))), :failed
        end
    end
end

function _row(point, algorithm, charge, result, n, flavor_profile, source, convergence, meson_point=nothing, conserved=nothing)
    flavor_mu = flavor_mu_profile_fm(flavor_profile, point.muq_MeV)
    charge_plus = charge === :plus
    meson_results = meson_point === nothing ? nothing : meson_point.meson_results
    eq = meson_point === nothing ? nothing : meson_point.equilibrium
    masses = meson_results === nothing ? (pi=NaN, K=NaN) :
        (pi=Float64(meson_results[:pi].mass) * ħc_MeV_fm, K=Float64(meson_results[:K].mass) * ħc_MeV_fm)
    rho = conserved === nothing ? (rho_u_fm3=NaN, rho_d_fm3=NaN, rho_s_quark_fm3=NaN, rho_B_fm3=NaN, rho_Q_fm3=NaN, rho_strange_charge_fm3=NaN) : conserved
    return (
        row_kind=convergence === nothing ? "freezeout" : "convergence",
        sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
        T_MeV=point.T_MeV,
        muB_MeV=point.muB_MeV,
        mu_u_MeV=flavor_mu.mu_u_MeV,
        mu_d_MeV=flavor_mu.mu_d_MeV,
        mu_s_MeV=flavor_mu.mu_s_MeV,
        mu_pi_MeV=charge_plus ? flavor_mu.mu_u_MeV - flavor_mu.mu_d_MeV : flavor_mu.mu_d_MeV - flavor_mu.mu_u_MeV,
        mu_K_MeV=charge_plus ? flavor_mu.mu_u_MeV - flavor_mu.mu_s_MeV : flavor_mu.mu_s_MeV - flavor_mu.mu_u_MeV,
        d_pi=1,
        d_K=1,
        path_profile="baseline_freezeout",
        flavor_profile=flavor_profile.profile_name,
        algorithm=String(algorithm),
        charge=String(charge),
        production_candidate_status="not_authorized",
        status=String(result.status),
        message=String(result.message),
        continuation_source=String(source),
        n_pi=result.n_pi,
        n_K=result.n_K,
        ratio=result.ratio,
        m_pi_MeV=masses.pi,
        m_K_MeV=masses.K,
        rho_u_fm3=rho.rho_u_fm3,
        rho_d_fm3=rho.rho_d_fm3,
        rho_s_quark_fm3=rho.rho_s_quark_fm3,
        rho_B_fm3=rho.rho_B_fm3,
        rho_Q_fm3=rho.rho_Q_fm3,
        rho_strange_charge_fm3=rho.rho_strange_charge_fm3,
        equilibrium_residual_norm=eq === nothing ? NaN : Float64(eq.residual_norm),
        equilibrium_converged=eq === nothing ? false : Bool(eq.converged),
        bose_pi_status=String(result.bose_pi_status),
        bose_K_status=String(result.bose_K_status),
        bose_pi_min_E_minus_mu=result.bose_pi_min_E_minus_mu,
        bose_K_min_E_minus_mu=result.bose_K_min_E_minus_mu,
        qmax_fm_inv=n.qmax,
        q_nodes=n.q_nodes,
        omega_min_fm_inv=n.omega_min,
        omega_max_fm_inv=n.omega_max,
        omega_nodes=n.omega_nodes,
        eta_fm_inv=n.eta,
        stable_q_nodes=n.stable_q_nodes,
        pole_iterations=n.pole_iterations,
        pole_residual_norm_max=n.pole_residual_norm_max,
        p_num=_env_int("MESON_CONVERGENCE_P_NUM", 8),
        t_num=_env_int("MESON_CONVERGENCE_T_NUM", 4),
        phase_measure=String(n.phase_measure),
        phase_anchor=String(n.phase_anchor),
        phase_density_policy="strict_normal_domain",
        convergence_reference_ratio=convergence === nothing ? NaN : convergence.reference,
        convergence_candidate_ratio=convergence === nothing ? NaN : convergence.candidate,
        convergence_relative_difference=convergence === nothing ? NaN : convergence.relative_difference,
        convergence_passed=convergence === nothing ? false : convergence.passed,
    )
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    columns = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(String.(columns), ','))
        for row in rows
            values = map(columns) do column
                value = _fmt(getproperty(row, column))
                occursin(r"[,\"\n]", value) ? string('"', replace(value, '"' => "\"\""), '"') : value
            end
            println(io, join(values, ','))
        end
    end
end

function _write_readme(path, rows, points, low, high)
    successful = filter(row -> row.status == "ok", rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Charged-RPA/BU convergence and freeze-out diagnostic")
        println(io)
        println(io, "- production candidate authorization: `not_authorized`")
        println(io, "- freeze-out path: `default` freezeout profile + `baseline_freezeout` path semantics")
        println(io, "- algorithms: `stable_particle_limit`, `reduced_strict_bw`, `q_pole_strict_bw`, `phase_shift_bu`")
        println(io, "- phase-shift strict candidate: `domega/pi` + `high_energy_zero`; legacy defaults are not changed")
        println(io, "- Bose support is a gate field, not a condensate treatment")
        println(io, "- successful rows: $(length(successful))/$(length(rows)); points: $(length(points))")
        if low !== nothing && high !== nothing
            println(io, "- representative convergence pair: low (q=$(low.q_nodes), omega=$(low.omega_nodes), qmax=$(low.qmax), omega_max=$(low.omega_max)) -> high (q=$(high.q_nodes), omega=$(high.omega_nodes), qmax=$(high.qmax), omega_max=$(high.omega_max))")
        end
        println(io, "- no strict-support, second-sheet pole, full Omega_M feedback, or production claim is made")
        println(io)
        println(io, "| sqrt(s) [GeV] | algorithm | charge | status | K/pi | Bose pi | Bose K | convergence |")
        println(io, "|---:|---|---|---|---:|---|---|---|")
        for row in sort(rows; by=r -> (r.sqrt_s_NN_GeV, r.algorithm, r.charge, r.row_kind))
            convergence = row.row_kind == "convergence" ? string(row.convergence_passed) : "-"
            println(io, "| $(row.sqrt_s_NN_GeV) | $(row.algorithm) | $(row.charge) | $(row.status) | $(row.ratio) | $(row.bose_pi_status) | $(row.bose_K_status) | $(convergence) |")
        end
    end
end

function main()
    output_dir = get(ENV, "MESON_CONVERGENCE_OUTPUT_DIR", DEFAULT_OUTPUT_DIR)
    points = build_freezeout_scan_points(_sqrts_values(); profile=load_freezeout_profile(profile="default"), traversal=:sqrts_descending)
    flavor_profile = load_flavor_chemical_profile(profile=get(ENV, "MESON_CONVERGENCE_FLAVOR_PROFILE", "default"))
    low = _validate_numerics(_numerics())
    high = _validate_numerics(_numerics(refined=true))
    algorithms = four_density_algorithm_labels()
    run_convergence = _env_bool("MESON_CONVERGENCE_RUN_CONVERGENCE", true)
    rows = NamedTuple[]
    continuation_state = nothing
    convergence_done = false

    for point in points
        solved, source = _solve_point(point, continuation_state, flavor_profile)
        if hasproperty(solved, :failure)
            for algorithm in algorithms, charge in (:plus, :minus)
                result = (status=:failed, message=String(solved.message), n_pi=NaN, n_K=NaN, ratio=NaN,
                          bose_pi_status=:not_evaluated, bose_K_status=:not_evaluated,
                          bose_pi_min_E_minus_mu=NaN, bose_K_min_E_minus_mu=NaN)
                push!(rows, _row(point, algorithm, charge, result, low, flavor_profile, source, nothing))
            end
            continue
        end
        continuation_state = solved.continuation_state
        flavor_mu = _chemical_potentials(flavor_profile, point)
        conserved = try
            _quark_conserved_densities(solved, low)
        catch
            nothing
        end

        for algorithm in algorithms, charge in (:plus, :minus)
            result = _safe_algorithm_density(algorithm, charge, solved, flavor_mu, low)
            push!(rows, _row(point, algorithm, charge, result, low, flavor_profile, source, nothing, solved, conserved))
            if run_convergence && !convergence_done
                refined_result = _safe_algorithm_density(algorithm, charge, solved, flavor_mu, high)
                gate = if result.status === :ok && refined_result.status === :ok && isfinite(result.ratio) && isfinite(refined_result.ratio)
                    convergence_gate(result.ratio, refined_result.ratio; rtol=_env_float("MESON_CONVERGENCE_RTOL", 5.0e-2), atol=_env_float("MESON_CONVERGENCE_ATOL", 1.0e-8))
                else
                    (passed=false, reference=result.ratio, candidate=refined_result.ratio, relative_difference=NaN)
                end
                push!(rows, _row(point, algorithm, charge, refined_result, high, flavor_profile, source, gate, solved, conserved))
            end
        end
        convergence_done |= run_convergence
        @printf("[charged-rpa-bu-convergence] sqrt(s)=%.3f GeV source=%s\n", point.sqrt_s_NN_GeV, source)
    end

    isempty(rows) && throw(ArgumentError("no convergence or freeze-out rows were produced"))
    csv_path = joinpath(output_dir, "charged_rpa_bu_convergence_freezeout.csv")
    readme_path = joinpath(output_dir, "README.md")
    _write_csv(csv_path, rows)
    _write_readme(readme_path, rows, points, low, high)
    println("[charged-rpa-bu-convergence] csv=$(csv_path)")
    println("[charged-rpa-bu-convergence] summary=$(readme_path)")
    return (rows=rows, csv=csv_path, summary=readme_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
