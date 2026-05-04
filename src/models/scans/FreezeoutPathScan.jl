"""
    FreezeoutPathScan

沿 chemical freeze-out baseline 路径组织 `FixedMu` 求解，并复用 TmuScan 的
continuation / multiseed / candidate governance 主链。
"""
module FreezeoutPathScan

using Main.Constants_PNJL: ħc_MeV_fm

using ..Models: FixedMu, build_seed_pool
using ..SeedStrategies: MultiSeed
using ..SeedStrategies: get_all_seeds, auto_phase_hint
using ..SeedStrategies: HADRON_SEED_5, QUARK_SEED_5, HT_GUESS_0p8_SEED_5, HT_GUESS_0p9_SEED_5, HT_GUESS_0p95_SEED_5, WEAK_CHIRAL_CONF_SEED_5
using ..ScanConfig: FreezeoutScanConfig, scan_kwargs
using ..FreezeoutProfiles
using ..FreezeoutPathProfiles
using ..TmuScan
using ..ScanCommon

export run_freezeout_fixedmu_scan, DEFAULT_FREEZEOUT_OUTPUT_PATH

const DEFAULT_FREEZEOUT_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "pnjl", "scan", "freezeout", "freezeout_fixedmu_scan.csv",
))

const HEADER = join((
    "sqrt_s_NN_GeV",
    "muB_MeV",
    "xi",
    "T_MeV",
    "muq_MeV",
    "freezeout_profile",
    "path_profile",
    "path_segment",
    "pressure_fm4",
    "rho",
    "entropy_fm3",
    "energy_fm4",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi1",
    "Phi2",
    "M_u_MeV",
    "M_d_MeV",
    "M_s_MeV",
    "iterations",
    "residual_norm",
    "converged",
    "message",
), ",")

@inline function _validate_real_vector(name::Symbol, values)
    values isa AbstractVector{<:Real} || throw(ArgumentError("$(name) must be AbstractVector{<:Real}, got $(typeof(values))"))
    isempty(values) && throw(ArgumentError("$(name) must not be empty"))
    for (i, v) in pairs(values)
        isfinite(Float64(v)) || throw(ArgumentError("$(name)[$(i)] must be finite Real, got $(v)"))
    end
    return nothing
end

@inline function _validate_inputs(sqrt_s_NN_values, xi_values, solver_backend::Symbol, model_kind::Symbol, traversal::Symbol)
    _validate_real_vector(:sqrt_s_NN_values, sqrt_s_NN_values)
    _validate_real_vector(:xi_values, xi_values)
    TmuScan._validate_tmu_scan_inputs([150.0], [0.0], xi_values, solver_backend, model_kind)
    TmuScan._validate_semantic_mode(:ground_state, nothing)
    TmuScan._validate_auto_pnjl_backend(:models)
    traversal in (:as_given, :sqrts_ascending, :sqrts_descending, :muB_descending, :muB_ascending) ||
        throw(ArgumentError("unsupported freezeout traversal: $(traversal)"))
    return nothing
end

@inline function _schedule_point_key(pt, xi)
    return ScanCommon.key3(pt.sqrt_s_NN_GeV, pt.muB_MeV, xi; digits=6)
end

@inline function _continuation_key(xi)
    return ScanCommon.key2(xi, 0.0; digits=6)
end

function _build_seed_candidates(cache::Dict, xi, T_fm, μ_fm; bootstrap_multiseed::Bool)
    provided_seed_pool = Vector{Vector{Float64}}()
    default_seed_pool = Vector{Vector{Float64}}()
    primary_seed = nothing

    seed_key = _continuation_key(xi)
    if haskey(cache, seed_key)
        primary_seed = Float64.(cache[seed_key])
    end

    if bootstrap_multiseed
        multiseeds = get_all_seeds(MultiSeed(), [T_fm, μ_fm], FixedMu())
        append!(provided_seed_pool, [Float64.(seed) for seed in multiseeds])
    end

    hint = auto_phase_hint(T_fm, μ_fm)
    if hint === :quark
        push!(default_seed_pool, Float64.(QUARK_SEED_5))
        push!(default_seed_pool, Float64.(HADRON_SEED_5))
    else
        push!(default_seed_pool, Float64.(HADRON_SEED_5))
        push!(default_seed_pool, Float64.(QUARK_SEED_5))
    end
    append!(default_seed_pool, (
        Float64.(WEAK_CHIRAL_CONF_SEED_5),
        Float64.(HT_GUESS_0p8_SEED_5),
        Float64.(HT_GUESS_0p9_SEED_5),
        Float64.(HT_GUESS_0p95_SEED_5),
    ))

    if primary_seed === nothing
        if !isempty(provided_seed_pool)
            primary_seed = provided_seed_pool[1]
            provided_seed_pool = provided_seed_pool[2:end]
        elseif !isempty(default_seed_pool)
            primary_seed = default_seed_pool[1]
            default_seed_pool = default_seed_pool[2:end]
        else
            primary_seed = Float64.(HADRON_SEED_5)
        end
    end

    seed_pool = build_seed_pool(FixedMu();
        primary_seed=primary_seed,
        provided_seed_pool=provided_seed_pool,
        default_seed_pool=default_seed_pool,
        seed_extend=(seed, _) -> Float64.(seed),
    )

    candidates = NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}[]
    seen_states = Set{String}()
    for (idx, entry) in enumerate(seed_pool)
        label = TmuScan._seed_pool_source_to_label(entry.source, idx)
        TmuScan._push_unique_candidate!(candidates, seen_states, label, entry.seed)
    end
    return candidates
end

function _write_row(io, pt, xi, profile_name::String, path_profile_name::String, path_segment::String, result, message)
    if result === nothing
        values = (
            ScanCommon.fmt(pt.sqrt_s_NN_GeV),
            ScanCommon.fmt(pt.muB_MeV),
            ScanCommon.fmt(xi),
            ScanCommon.fmt(pt.T_MeV),
            ScanCommon.fmt(pt.muq_MeV),
            profile_name,
            path_profile_name,
            path_segment,
            "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
            "-1", "NaN", "false",
            ScanCommon.quote_csv(message),
        )
        println(io, join(values, ','))
        return
    end

    phi = result.x_state[1:3]
    Phi1 = result.x_state[4]
    Phi2 = result.x_state[5]
    masses_mev = result.masses .* ħc_MeV_fm

    values = (
        ScanCommon.fmt(pt.sqrt_s_NN_GeV),
        ScanCommon.fmt(pt.muB_MeV),
        ScanCommon.fmt(xi),
        ScanCommon.fmt(pt.T_MeV),
        ScanCommon.fmt(pt.muq_MeV),
        profile_name,
        path_profile_name,
        path_segment,
        ScanCommon.fmt(result.pressure),
        ScanCommon.fmt(result.rho_norm),
        ScanCommon.fmt(result.entropy),
        ScanCommon.fmt(result.energy),
        ScanCommon.fmt(phi[1]),
        ScanCommon.fmt(phi[2]),
        ScanCommon.fmt(phi[3]),
        ScanCommon.fmt(Phi1),
        ScanCommon.fmt(Phi2),
        ScanCommon.fmt(masses_mev[1]),
        ScanCommon.fmt(masses_mev[2]),
        ScanCommon.fmt(masses_mev[3]),
        string(result.iterations),
        ScanCommon.fmt(result.residual_norm),
        string(result.converged),
        ScanCommon.quote_csv(message),
    )
    println(io, join(values, ','))
end

function run_freezeout_fixedmu_scan(;
    sqrt_s_NN_values,
    xi_values=[0.0],
    profile_name::AbstractString="default",
    path_profile_name::AbstractString="baseline_freezeout",
    output_path::AbstractString=DEFAULT_FREEZEOUT_OUTPUT_PATH,
    overwrite::Bool=false,
    resume::Bool=true,
    bootstrap_multiseed::Bool=true,
    solver_backend::Symbol=:auto,
    auto_pnjl_backend::Symbol=:models,
    semantic_mode::Symbol=:ground_state,
    selector::Union{Nothing, Function}=nothing,
    model_kind::Symbol=:PNJL,
    traversal::Symbol=:sqrts_ascending,
    p_num::Int=24,
    t_num::Int=8,
    progress_cb::Union{Nothing, Function}=nothing,
    diagnostic_level::Symbol=:none,
    nlsolve_kwargs...
)
    _validate_inputs(sqrt_s_NN_values, xi_values, solver_backend, model_kind, traversal)
    TmuScan._validate_semantic_mode(semantic_mode, selector)
    TmuScan._validate_auto_pnjl_backend(auto_pnjl_backend)

    profile = FreezeoutProfiles.load_freezeout_profile(profile=String(profile_name))
    path_profile = FreezeoutPathProfiles.load_freezeout_path_profile(profile=String(path_profile_name))
    points = FreezeoutPathProfiles.build_freezeout_path_points(
        sqrt_s_NN_values;
        freezeout_profile=profile,
        path_profile=path_profile,
        traversal=traversal,
    )

    mkpath(dirname(output_path))
    completed = (resume && !overwrite && isfile(output_path)) ? ScanCommon.load_completed_keys3(output_path; digits=6) : Set{NTuple{3, Float64}}()
    io_mode = (overwrite || !isfile(output_path)) ? "w" : "a"
    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)
    continuation_seeds = Dict{Tuple{Float64, Float64}, Vector{Float64}}()

    open(output_path, io_mode) do io
        if io_mode == "w"
            println(io, HEADER)
        end

        for xi in xi_values
            delete!(continuation_seeds, _continuation_key(xi))
            for pt in points
                stats[:total] += 1
                key = _schedule_point_key(pt, xi)
                if key in completed
                    stats[:skipped] += 1
                    continue
                end

                candidates = _build_seed_candidates(continuation_seeds, xi, pt.T_fm, pt.muB_fm / 3.0; bootstrap_multiseed=bootstrap_multiseed)
                result, message = ScanCommon.attempt_with_candidates(candidates;
                    solve_point=seed_state -> TmuScan._solve_point(pt.T_fm, pt.muB_fm / 3.0, xi, seed_state;
                        solver_backend=solver_backend,
                        auto_pnjl_backend=auto_pnjl_backend,
                        model_kind=model_kind,
                        diagnostic_level=diagnostic_level,
                        p_num=p_num,
                        t_num=t_num,
                        nlsolve_kwargs...,
                    ),
                    refine=result -> TmuScan._refine_result(pt.T_fm, pt.muB_fm / 3.0, xi, result;
                        solver_backend=solver_backend,
                        auto_pnjl_backend=auto_pnjl_backend,
                        model_kind=model_kind,
                        diagnostic_level=diagnostic_level,
                        p_num=p_num,
                        t_num=t_num,
                        nlsolve_kwargs...,
                    ),
                    promote=TmuScan._promote_success,
                    is_success=TmuScan._is_success,
                    stop_on_first_success=false,
                )

                if result !== nothing && TmuScan._is_success(result)
                    continuation_seeds[_continuation_key(xi)] = copy(result.solution)
                end

                _write_row(io, pt, xi, profile.profile_name, pt.path_profile, pt.path_segment, result, message)
                flush(io)
                push!(completed, key)

                if TmuScan._is_success(result)
                    stats[:success] += 1
                else
                    stats[:failure] += 1
                end

                if progress_cb !== nothing
                    try
                        progress_cb((sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV, T_MeV=pt.T_MeV, muB_MeV=pt.muB_MeV, xi=xi), result)
                    catch
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
        output=output_path,
        freezeout_profile=profile.profile_name,
        path_profile=path_profile.profile_name,
        traversal=traversal,
    )
end

function run_freezeout_fixedmu_scan(config::FreezeoutScanConfig; kwargs...)
    cfg = scan_kwargs(config)
    return run_freezeout_fixedmu_scan(; cfg..., config.nlsolve_kwargs..., kwargs...)
end

end # module FreezeoutPathScan
