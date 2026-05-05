module ServerWarmup

export list_server_warmup_steps, resolve_server_warmup_profile
export run_server_warmup, run_server_warmup_from_env

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(Main, :Models)
    Base.include(Main, joinpath(_PROJECT_ROOT, "src", "models", "Models.jl"))
end
using Main: Models

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, joinpath(_PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end
using Main.Constants_PNJL: ħc_MeV_fm

@inline function _parse_bool_flag(value, default::Bool)
    value === nothing && return default
    if value isa Bool
        return value
    end
    lowered = lowercase(strip(String(value)))
    if lowered in ("1", "true", "yes", "on")
        return true
    elseif lowered in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("invalid boolean flag: $(value)"))
end

function resolve_server_warmup_profile(profile::AbstractString)
    normalized = lowercase(strip(String(profile)))
    normalized in ("none", "point", "service_core") || throw(ArgumentError("unknown server warmup profile: $(profile)"))
    return Symbol(normalized)
end

function list_server_warmup_steps(profile::Symbol)
    if profile === :none
        return Symbol[]
    elseif profile === :point
        return [:pnjl_point, :transport_point]
    elseif profile === :service_core
        return [:pnjl_point, :transport_point, :phase_pipeline]
    end
    throw(ArgumentError("unknown server warmup profile: $(profile)"))
end

function _run_warmup_step(step::Symbol)
    if step === :pnjl_point
        Models.solve_pnjl_point(T_mev=150.0, mu_mev=0.0, xi=0.0, p_num=6, t_num=3)
        return nothing
    elseif step === :transport_point
        Models.run_precompile_capability(:transport_point_api; strict=true)
        return nothing
    elseif step === :phase_pipeline
        outdir = mktempdir()
        Models.run_phase_pipeline(
            :PNJL;
            mode=:research,
            T_grid=[150.0],
            rho_grid=[0.1, 0.2],
            xi=0.0,
            output_dir=outdir,
            profile=:smoke,
            solver_backend=:models,
            reverse_rho=true,
            seed_policy=:hybrid_continuity,
            p_num=6,
            t_num=3,
            iterations=4,
        )
        return nothing
    end
    throw(ArgumentError("unknown server warmup step: $(step)"))
end

function run_server_warmup(profile::Symbol=:none; strict::Bool=false)
    for step in list_server_warmup_steps(profile)
        try
            _run_warmup_step(step)
        catch err
            if strict
                rethrow(err)
            end
            @warn "server warmup step failed" profile=profile step=step exception=(err, catch_backtrace())
        end
    end
    return nothing
end

function run_server_warmup_from_env()
    profile = resolve_server_warmup_profile(get(ENV, "JRT_SERVER_WARMUP_PROFILE", "none"))
    strict = _parse_bool_flag(get(ENV, "JRT_SERVER_WARMUP_STRICT", "0"), false)
    run_server_warmup(profile; strict=strict)
    return nothing
end

end
