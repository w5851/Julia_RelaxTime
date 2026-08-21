"""
    MagneticScan

Production-facing fixed-chemical-potential scan for the magnetic PNJL route.

The magnetic route is intentionally separate from `TmuScan`/`TrhoScan`:
nonzero `eB` requires the dedicated five-field magnetic Omega and its
branch-aware solver.  This module accepts a common flavor chemical potential
`mu_values`, keeps every converged candidate, and writes selected and
candidate-level CSV artifacts.
"""
module MagneticScan

using Printf
using StaticArrays

using Main.Constants_PNJL: ħc_MeV_fm
using ..Models: MagneticGapCandidate, PNJLMagneticModel, create_model
using ..Models: calculate_mass_vec, model_thermo, number_densities, solve_magnetic_gap

export run_magnetic_scan, DEFAULT_T_VALUES, DEFAULT_MU_VALUES, DEFAULT_EB_VALUES
export DEFAULT_OUTPUT_PATH, DEFAULT_CANDIDATES_OUTPUT_PATH

const DEFAULT_T_VALUES = collect(50.0:10.0:200.0)
const DEFAULT_MU_VALUES = collect(0.0:10.0:400.0)
const DEFAULT_EB_VALUES = [0.0, 1.0e4, 2.0e4]
const DEFAULT_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "scan", "magnetic", "magnetic_scan.csv",
))
const DEFAULT_CANDIDATES_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "scan", "magnetic", "magnetic_scan_candidates.csv",
))
const KEY_DIGITS = 8

const SELECTED_HEADER = join((
    "T_MeV",
    "mu_MeV",
    "eB_MeV2",
    "xi",
    "selected_candidate_index",
    "branch_label",
    "method",
    "attempt_count",
    "failed_attempts",
    "omega_fm4",
    "pressure_fm4",
    "rho_norm",
    "entropy_fm3",
    "energy_fm4",
    "rho_u_fm3",
    "rho_d_fm3",
    "rho_s_fm3",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi",
    "PhiBar",
    "M_u_MeV",
    "M_d_MeV",
    "M_s_MeV",
    "residual_norm",
    "n_max",
    "converged",
    "message",
), ",")

const CANDIDATE_HEADER = join((
    "T_MeV",
    "mu_MeV",
    "eB_MeV2",
    "xi",
    "candidate_index",
    "seed_index",
    "branch_label",
    "stability",
    "method",
    "iterations",
    "omega_fm4",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi",
    "PhiBar",
    "M_u_MeV",
    "M_d_MeV",
    "M_s_MeV",
    "residual_u",
    "residual_d",
    "residual_s",
    "residual_Phi",
    "residual_PhiBar",
    "residual_norm",
    "n_max",
    "physical",
    "converged",
    "message",
), ",")

@inline function _validate_real_vector(name::Symbol, values)
    values isa AbstractVector{<:Real} || throw(ArgumentError(
        "$(name) must be AbstractVector{<:Real}, got $(typeof(values))",
    ))
    isempty(values) && throw(ArgumentError("$(name) must not be empty"))
    for (index, value) in pairs(values)
        isfinite(Float64(value)) || throw(ArgumentError(
            "$(name)[$(index)] must be finite, got $(value)",
        ))
    end
    return nothing
end

@inline function _validate_magnetic_scan_inputs(
    T_values,
    mu_values,
    eB_values,
    xi_values,
    model_kind::Symbol,
    solver_mode::Symbol,
    p_num::Int,
    t_num::Int,
    iterations::Int,
)
    model_kind === :PNJLMagnetic || throw(ArgumentError(
        "run_magnetic_scan only supports model_kind=:PNJLMagnetic; got $(model_kind)",
    ))
    solver_mode === :fixed_mu || throw(ArgumentError(
        "run_magnetic_scan only supports solver_mode=:fixed_mu; got $(solver_mode)",
    ))
    _validate_real_vector(:T_values, T_values)
    _validate_real_vector(:mu_values, mu_values)
    _validate_real_vector(:eB_values, eB_values)
    _validate_real_vector(:xi_values, xi_values)
    all(Float64(T) > 0.0 for T in T_values) || throw(ArgumentError(
        "magnetic scan requires T_values > 0 MeV",
    ))
    all(abs(Float64(xi)) <= 1e-14 for xi in xi_values) || throw(ArgumentError(
        "PNJLMagneticModel currently supports only xi=0",
    ))
    p_num >= 4 || throw(ArgumentError("p_num must be >= 4, got $(p_num)"))
    t_num >= 1 || throw(ArgumentError("t_num must be >= 1, got $(t_num)"))
    iterations >= 1 || throw(ArgumentError("iterations must be >= 1, got $(iterations)"))
    return nothing
end

@inline _scan_key(T, mu, eB, xi) = (
    round(Float64(T); digits=KEY_DIGITS),
    round(Float64(mu); digits=KEY_DIGITS),
    round(Float64(eB); digits=KEY_DIGITS),
    round(Float64(xi); digits=KEY_DIGITS),
)

function _load_completed_keys(path::AbstractString)
    completed = Set{NTuple{4, Float64}}()
    isfile(path) || return completed
    for (line_number, line) in enumerate(readlines(path))
        line_number == 1 && continue
        fields = split(line, ','; limit=5)
        length(fields) >= 4 || continue
        try
            push!(completed, _scan_key(
                parse(Float64, fields[1]),
                parse(Float64, fields[2]),
                parse(Float64, fields[3]),
                parse(Float64, fields[4]),
            ))
        catch
            # Leave malformed historical rows for the normal output validator.
        end
    end
    return completed
end

@inline _fmt(value) = value isa AbstractFloat && !isfinite(value) ? "NaN" : @sprintf("%.12g", Float64(value))
@inline _fmt(value::Integer) = string(value)
@inline _fmt(value::Bool) = value ? "true" : "false"
@inline _fmt(value::Symbol) = String(value)

function _quote(value)
    text = String(value)
    return occursin(r"[,\"\r\n]", text) ? "\"" * replace(text, '"' => "\"\"") * "\"" : text
end

@inline _field(value) = _fmt(value)
@inline _field(value::AbstractString) = _quote(value)
@inline _field(value::Symbol) = _quote(value)

function _candidate_mass_mev(model::PNJLMagneticModel, candidate::MagneticGapCandidate)
    masses = calculate_mass_vec(model, candidate.x_state[1:3])
    return masses .* ħc_MeV_fm
end

function _write_failed_selected!(io, T, mu, eB, xi, message)
    values = (
        T, mu, eB, xi, 0, :none, :none, 0, 0,
        (fill(NaN, 17)...),
        0, false, message,
    )
    println(io, join((_field(value) for value in values), ','))
end

function _write_selected!(io, T, mu, eB, xi, result, candidate, thermo, masses, rho_vec, message)
    values = (
        T, mu, eB, xi, result.selected_index, candidate.branch_label, candidate.method,
        result.attempt_count, result.failed_attempts, candidate.omega,
        -candidate.omega, thermo[2], thermo[3], thermo[4], rho_vec[1], rho_vec[2], rho_vec[3],
        candidate.x_state[1], candidate.x_state[2], candidate.x_state[3],
        candidate.x_state[4], candidate.x_state[5], masses[1], masses[2], masses[3],
        candidate.residual_norm, candidate.n_max, result.converged, message,
    )
    println(io, join((_field(value) for value in values), ','))
end

function _write_candidate!(io, T, mu, eB, xi, index, candidate, model, message)
    masses = _candidate_mass_mev(model, candidate)
    values = (
        T, mu, eB, xi, index, candidate.seed_index, candidate.branch_label,
        candidate.stability, candidate.method, candidate.iterations, candidate.omega,
        candidate.x_state[1], candidate.x_state[2], candidate.x_state[3],
        candidate.x_state[4], candidate.x_state[5], masses[1], masses[2], masses[3],
        candidate.residual[1], candidate.residual[2], candidate.residual[3],
        candidate.residual[4], candidate.residual[5], candidate.residual_norm,
        candidate.n_max, candidate.physical, candidate.converged, message,
    )
    println(io, join((_field(value) for value in values), ','))
end

function _derived_candidates_output(output_path::AbstractString)
    base, extension = splitext(String(output_path))
    return base * "_candidates" * extension
end

"""
    run_magnetic_scan(; kwargs...) -> NamedTuple

Run a fixed-`mu` magnetic PNJL scan over `(T, mu, eB)`.  `T` and `mu` are
external MeV values, `eB` is in MeV^2, and `mu` is applied equally to the
three quark flavors.  Only `xi=0` is accepted.  Every distinct converged
candidate is written to `candidates_output_path`; the convenience selected
state and its thermodynamics are written to `output_path`.
"""
function run_magnetic_scan(
    ;
    T_values=DEFAULT_T_VALUES,
    mu_values=DEFAULT_MU_VALUES,
    eB_values=DEFAULT_EB_VALUES,
    xi_values=[0.0],
    model_kind::Symbol=:PNJLMagnetic,
    solver_mode::Symbol=:fixed_mu,
    output_path::AbstractString=DEFAULT_OUTPUT_PATH,
    candidates_output_path::Union{Nothing, AbstractString}=nothing,
    overwrite::Bool=false,
    resume::Bool=true,
    profile::String=get(ENV, "PNJL_PARAM_PROFILE", "default"),
    physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
    p_num::Int=96,
    t_num::Int=8,
    pz_max::Real=0.0,
    n_max::Union{Nothing, Int}=nothing,
    cutoff_N::Int=10,
    method::Symbol=:trust_region,
    fallback_method::Union{Nothing, Symbol}=:newton,
    xtol::Real=1e-9,
    ftol::Real=1e-9,
    iterations::Int=600,
    residual_norm_max::Real=1e-6,
    finite_difference_step::Real=1e-5,
    root_merge_tol::Real=1e-5,
    classify_stability::Bool=false,
    include_default_seeds::Bool=true,
    progress_cb::Union{Nothing, Function}=nothing,
)
    _validate_magnetic_scan_inputs(
        T_values, mu_values, eB_values, xi_values, model_kind, solver_mode,
        p_num, t_num, iterations,
    )
    candidate_path = candidates_output_path === nothing ?
        _derived_candidates_output(output_path) : String(candidates_output_path)
    mkpath(dirname(output_path))
    mkpath(dirname(candidate_path))

    completed = (resume && !overwrite) ? _load_completed_keys(output_path) : Set{NTuple{4, Float64}}()
    selected_mode = (overwrite || !isfile(output_path)) ? "w" : "a"
    candidates_mode = (overwrite || !isfile(candidate_path)) ? "w" : "a"
    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)
    continuation = Dict{Tuple{Float64, Float64}, SVector{5, Float64}}()

    open(output_path, selected_mode) do selected_io
        open(candidate_path, candidates_mode) do candidates_io
            selected_mode == "w" && println(selected_io, SELECTED_HEADER)
            candidates_mode == "w" && println(candidates_io, CANDIDATE_HEADER)

            for xi in xi_values
                for T in T_values
                    for eB in eB_values
                        T_value = Float64(T)
                        mu_key = round(T_value; digits=KEY_DIGITS)
                        eB_key = round(Float64(eB); digits=KEY_DIGITS)
                        seed_key = (mu_key, eB_key)
                        delete!(continuation, seed_key)
                        eB_fm2 = Float64(eB) / ħc_MeV_fm^2
                        model = create_model(
                            :PNJLMagnetic;
                            eB_fm2=eB_fm2,
                            profile=profile,
                            physics_profile=physics_profile,
                            p_num=p_num,
                            pz_max=pz_max,
                            n_max=n_max,
                            cutoff_N=cutoff_N,
                        )

                        for mu in mu_values
                            stats[:total] += 1
                            key = _scan_key(T, mu, eB, xi)
                            if key in completed
                                stats[:skipped] += 1
                                continue
                            end

                            T_fm = T_value / ħc_MeV_fm
                            mu_fm = Float64(mu) / ħc_MeV_fm
                            initial_guess = get(continuation, seed_key, nothing)
                            message = ""
                            result = nothing
                            try
                                result = solve_magnetic_gap(
                                    model,
                                    T_fm,
                                    SVector{3, Float64}(mu_fm, mu_fm, mu_fm);
                                    xi=Float64(xi),
                                    p_num=p_num,
                                    t_num=t_num,
                                    pz_max=pz_max,
                                    n_max=n_max,
                                    cutoff_N=cutoff_N,
                                    initial_guess=initial_guess,
                                    method=method,
                                    fallback_method=fallback_method,
                                    xtol=xtol,
                                    ftol=ftol,
                                    iterations=iterations,
                                    residual_norm_max=residual_norm_max,
                                    finite_difference_step=finite_difference_step,
                                    root_merge_tol=root_merge_tol,
                                    classify_stability=classify_stability,
                                    include_default_seeds=include_default_seeds,
                                )
                                candidate = result.candidates[result.selected_index]
                                thermo = model_thermo(
                                    model,
                                    candidate.x_state,
                                    SVector{3, Float64}(mu_fm, mu_fm, mu_fm),
                                    T_fm;
                                    p_num=p_num,
                                    t_num=t_num,
                                    pz_max=pz_max,
                                    n_max=n_max,
                                    cutoff_N=cutoff_N,
                                    xi=Float64(xi),
                                )
                                densities = number_densities(
                                    model,
                                    candidate.x_state,
                                    T_fm,
                                    SVector{3, Float64}(mu_fm, mu_fm, mu_fm);
                                    p_num=p_num,
                                    t_num=t_num,
                                    pz_max=pz_max,
                                    n_max=n_max,
                                    cutoff_N=cutoff_N,
                                    xi=Float64(xi),
                                )
                                rho_vec = hasproperty(densities, :net) ? densities.net :
                                    densities.quark .- densities.antiquark
                                masses = _candidate_mass_mev(model, candidate)
                                _write_selected!(
                                    selected_io, T, mu, eB, xi, result, candidate,
                                    thermo, masses, rho_vec, message,
                                )
                                for (index, candidate_row) in enumerate(result.candidates)
                                    _write_candidate!(
                                        candidates_io, T, mu, eB, xi, index,
                                        candidate_row, model, message,
                                    )
                                end
                                flush(selected_io)
                                flush(candidates_io)
                                continuation[seed_key] = candidate.x_state
                                stats[:success] += 1
                            catch err
                                err isa InterruptException && rethrow()
                                message = sprint(showerror, err)
                                _write_failed_selected!(selected_io, T, mu, eB, xi, message)
                                flush(selected_io)
                                stats[:failure] += 1
                            end
                            push!(completed, key)
                            if progress_cb !== nothing
                                try
                                    progress_cb((T=T, mu=mu, eB=eB, xi=xi), result)
                                catch
                                    # Progress reporting must not alter scan semantics.
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return (
        total=stats[:total],
        success=stats[:success],
        failure=stats[:failure],
        skipped=stats[:skipped],
        output=String(output_path),
        selected_output=String(output_path),
        candidates_output=String(candidate_path),
    )
end

end # module MagneticScan
