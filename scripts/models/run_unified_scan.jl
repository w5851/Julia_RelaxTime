#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PHASE_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

Base.include(@__MODULE__, joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
Base.include(@__MODULE__, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

function _parse_bool(raw::AbstractString, name::AbstractString)
    lowered = lowercase(strip(raw))
    if lowered in ("1", "true", "yes", "on")
        return true
    elseif lowered in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("invalid --$(name)=$(raw), accepted: true|false"))
end

function _parse_real_list(raw::AbstractString, name::AbstractString)
    values = Float64[]
    for token in split(raw, ',')
        cleaned = strip(token)
        isempty(cleaned) && continue
        push!(values, parse(Float64, cleaned))
    end
    isempty(values) && throw(ArgumentError("invalid --$(name)=$(raw), expected comma-separated numbers"))
    return values
end

function _parse_symbol(raw::AbstractString)
    return Symbol(strip(raw))
end

function _parse_model_kind(raw::AbstractString)
    normalized = uppercase(strip(raw))
    if normalized == "PNJL_ANISO"
        return :pnjl_aniso
    end
    return Symbol(normalized)
end

function _parse_tmu_args(args::Vector{String})
    kwargs = Dict{Symbol, Any}()
    for arg in args
        startswith(arg, "--") || throw(ArgumentError("unknown argument: $(arg)"))
        parts = split(arg[3:end], '='; limit=2)
        length(parts) == 2 || throw(ArgumentError("invalid argument format: $(arg), expected --key=value"))
        key = parts[1]
        value = parts[2]
        if key == "model_kind"
            kwargs[:model_kind] = _parse_model_kind(value)
        elseif key == "T_values"
            kwargs[:T_values] = _parse_real_list(value, key)
        elseif key == "mu_values"
            kwargs[:mu_values] = _parse_real_list(value, key)
        elseif key == "xi_values"
            kwargs[:xi_values] = _parse_real_list(value, key)
        elseif key == "output_path"
            kwargs[:output_path] = value
        elseif key == "overwrite"
            kwargs[:overwrite] = _parse_bool(value, key)
        elseif key == "resume"
            kwargs[:resume] = _parse_bool(value, key)
        elseif key == "use_phase_aware"
            kwargs[:use_phase_aware] = _parse_bool(value, key)
        elseif key == "solver_backend"
            kwargs[:solver_backend] = Symbol(lowercase(value))
        elseif key == "p_num"
            kwargs[:p_num] = parse(Int, value)
        elseif key == "t_num"
            kwargs[:t_num] = parse(Int, value)
        elseif key == "iterations"
            kwargs[:iterations] = parse(Int, value)
        else
            throw(ArgumentError("unknown scan tmu option: --$(key)"))
        end
    end
    return kwargs
end

function _parse_trho_args(args::Vector{String})
    kwargs = Dict{Symbol, Any}()
    for arg in args
        startswith(arg, "--") || throw(ArgumentError("unknown argument: $(arg)"))
        parts = split(arg[3:end], '='; limit=2)
        length(parts) == 2 || throw(ArgumentError("invalid argument format: $(arg), expected --key=value"))
        key = parts[1]
        value = parts[2]
        if key == "model_kind"
            kwargs[:model_kind] = _parse_model_kind(value)
        elseif key == "T_values"
            kwargs[:T_values] = _parse_real_list(value, key)
        elseif key == "rho_values"
            kwargs[:rho_values] = _parse_real_list(value, key)
        elseif key == "xi_values"
            kwargs[:xi_values] = _parse_real_list(value, key)
        elseif key == "output_path"
            kwargs[:output_path] = value
        elseif key == "overwrite"
            kwargs[:overwrite] = _parse_bool(value, key)
        elseif key == "resume"
            kwargs[:resume] = _parse_bool(value, key)
        elseif key == "reverse_rho"
            kwargs[:reverse_rho] = _parse_bool(value, key)
        elseif key == "seed_policy"
            kwargs[:seed_policy] = _parse_symbol(lowercase(value))
        elseif key == "constraint_mode"
            kwargs[:constraint_mode] = _parse_symbol(lowercase(value))
        elseif key == "solver_backend"
            kwargs[:solver_backend] = _parse_symbol(lowercase(value))
        elseif key == "p_num"
            kwargs[:p_num] = parse(Int, value)
        elseif key == "t_num"
            kwargs[:t_num] = parse(Int, value)
        elseif key == "iterations"
            kwargs[:iterations] = parse(Int, value)
        else
            throw(ArgumentError("unknown scan trho option: --$(key)"))
        end
    end
    return kwargs
end

function _run_scan_tmu(args::Vector{String})
    kwargs = _parse_tmu_args(args)
    stats = Main.Models.run_tmu_scan(; kwargs...)
    println("[scan tmu] total=$(stats.total) success=$(stats.success) failure=$(stats.failure) skipped=$(stats.skipped)")
    println("[scan tmu] output=$(stats.output)")
    return nothing
end

function _run_scan_trho(args::Vector{String})
    kwargs = _parse_trho_args(args)
    stats = Main.Models.run_trho_scan(; kwargs...)
    println("[scan trho] total=$(stats.total) success=$(stats.success) failure=$(stats.failure) skipped=$(stats.skipped)")
    println("[scan trho] output=$(stats.output)")
    return nothing
end

function _run_workflow_phase(args::Vector{String})
    isfile(PHASE_SCRIPT) || throw(ArgumentError("phase workflow script not found: $(PHASE_SCRIPT)"))
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(PHASE_SCRIPT) $(args)`
    run(cmd)
    return nothing
end

function _usage()
    println("Usage:")
    println("  julia scripts/models/run_unified_scan.jl scan tmu [--key=value ...]")
    println("  julia scripts/models/run_unified_scan.jl scan trho [--key=value ...]")
    println("  julia scripts/models/run_unified_scan.jl workflow phase [--key=value ...]")
end

function main(args=ARGS)
    isempty(args) && throw(ArgumentError("missing subcommand; expected: scan|workflow"))
    if args[1] in ("-h", "--help")
        _usage()
        return nothing
    end
    length(args) >= 2 || throw(ArgumentError("missing second-level subcommand"))

    if args[1] == "scan"
        if args[2] == "tmu"
            return _run_scan_tmu(args[3:end])
        elseif args[2] == "trho"
            return _run_scan_trho(args[3:end])
        end
        throw(ArgumentError("unknown scan subcommand: $(args[2]); expected: tmu|trho"))
    elseif args[1] == "workflow"
        if args[2] == "phase"
            return _run_workflow_phase(args[3:end])
        end
        throw(ArgumentError("unknown workflow subcommand: $(args[2]); expected: phase"))
    end

    throw(ArgumentError("unknown subcommand: $(args[1]); expected: scan|workflow"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
