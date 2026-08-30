#!/usr/bin/env julia

"""
Status: core-candidate
Purpose: unified scan/workflow CLI entry for current `Models`-based pipelines.
Replacement: authoritative current entry; prefer this over removed legacy scan wrappers such as `scripts/pnjl/run_tmu_scan.jl`.
"""

using Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PHASE_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

Pkg.activate(PROJECT_ROOT)

Base.include(@__MODULE__, joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
Base.include(@__MODULE__, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

function _parse_bool(raw::AbstractString, name::AbstractString)
    lowered = lowercase(strip(raw))
    if lowered in ("1", "true", "yes", "on")
        return true
    elseif lowered in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("invalid --$(name)=$(raw), accepted (case-insensitive, trimmed): 1|0|true|false|yes|no|on|off"))
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
    raw_trimmed = strip(raw)
    lowered = lowercase(raw_trimmed)
    if lowered == "pnjl_aniso"
        return :pnjl_aniso
    end

    for kind in Models.registered_model_kinds()
        if lowercase(String(kind)) == lowered
            return kind
        end
    end

    return Symbol(raw_trimmed)
end

function _assert_required_keys(kwargs::Dict{Symbol, Any}, keys::Tuple{Vararg{Symbol}}, command::AbstractString)
    missing = Symbol[]
    for key in keys
        haskey(kwargs, key) || push!(missing, key)
    end
    if !isempty(missing)
        throw(ArgumentError("missing required options for $(command): $(join(string.(missing), ", "))"))
    end
    return nothing
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

function _parse_magnetic_args(args::Vector{String})
    kwargs = Dict{Symbol, Any}()
    for arg in args
        startswith(arg, "--") || throw(ArgumentError("unknown argument: $(arg)"))
        parts = split(arg[3:end], '='; limit=2)
        length(parts) == 2 || throw(ArgumentError("invalid argument format: $(arg), expected --key=value"))
        key = parts[1]
        value = parts[2]
        if key == "model_kind"
            kwargs[:model_kind] = _parse_model_kind(value)
        elseif key == "solver_mode"
            kwargs[:solver_mode] = Symbol(lowercase(value))
        elseif key == "T_values"
            kwargs[:T_values] = _parse_real_list(value, key)
        elseif key == "mu_values"
            kwargs[:mu_values] = _parse_real_list(value, key)
        elseif key == "eB_values"
            kwargs[:eB_values] = _parse_real_list(value, key)
        elseif key == "xi_values"
            kwargs[:xi_values] = _parse_real_list(value, key)
        elseif key == "output_path"
            kwargs[:output_path] = value
        elseif key == "candidates_output_path"
            kwargs[:candidates_output_path] = value
        elseif key == "overwrite"
            kwargs[:overwrite] = _parse_bool(value, key)
        elseif key == "resume"
            kwargs[:resume] = _parse_bool(value, key)
        elseif key == "profile"
            kwargs[:profile] = value
        elseif key == "physics_profile"
            kwargs[:physics_profile] = value
        elseif key == "magnetic_profile"
            kwargs[:magnetic_profile] = value
        elseif key == "p_num"
            kwargs[:p_num] = parse(Int, value)
        elseif key == "t_num"
            kwargs[:t_num] = parse(Int, value)
        elseif key == "pz_max"
            kwargs[:pz_max] = parse(Float64, value)
        elseif key == "n_max"
            kwargs[:n_max] = parse(Int, value)
        elseif key == "n_max_policy"
            kwargs[:n_max_policy] = Symbol(lowercase(value))
        elseif key == "thermal_tail_factor"
            kwargs[:thermal_tail_factor] = parse(Float64, value)
        elseif key == "n_max_floor"
            kwargs[:n_max_floor] = parse(Int, value)
        elseif key == "n_max_cap"
            kwargs[:n_max_cap] = parse(Int, value)
        elseif key == "cutoff_N"
            kwargs[:cutoff_N] = parse(Int, value)
        elseif key == "method"
            kwargs[:method] = Symbol(lowercase(value))
        elseif key == "fallback_method"
            kwargs[:fallback_method] = lowercase(value) == "none" ? nothing : Symbol(lowercase(value))
        elseif key == "iterations"
            kwargs[:iterations] = parse(Int, value)
        elseif key == "classify_stability"
            kwargs[:classify_stability] = _parse_bool(value, key)
        else
            throw(ArgumentError("unknown scan magnetic option: --$(key)"))
        end
    end
    return kwargs
end

function _run_scan_tmu(args::Vector{String})
    kwargs = _parse_tmu_args(args)
    _assert_required_keys(kwargs, (:model_kind, :T_values, :mu_values, :xi_values, :output_path), "scan tmu")
    stats = Models.run_scan_pipeline(:tmu; kwargs...)
    println("[scan tmu] total=$(stats.total) success=$(stats.success) failure=$(stats.failure) skipped=$(stats.skipped)")
    println("[scan tmu] output=$(stats.output)")
    return nothing
end

function _run_scan_trho(args::Vector{String})
    kwargs = _parse_trho_args(args)
    _assert_required_keys(kwargs, (:model_kind, :T_values, :rho_values, :xi_values, :output_path), "scan trho")
    stats = Models.run_scan_pipeline(:trho; kwargs...)
    println("[scan trho] total=$(stats.total) success=$(stats.success) failure=$(stats.failure) skipped=$(stats.skipped)")
    println("[scan trho] output=$(stats.output)")
    return nothing
end

function _run_scan_magnetic(args::Vector{String})
    kwargs = _parse_magnetic_args(args)
    _assert_required_keys(kwargs, (:model_kind, :T_values, :mu_values, :eB_values, :xi_values, :output_path), "scan magnetic")
    stats = Models.run_scan_pipeline(:magnetic; kwargs...)
    println("[scan magnetic] total=$(stats.total) success=$(stats.success) failure=$(stats.failure) skipped=$(stats.skipped)")
    println("[scan magnetic] selected_output=$(stats.selected_output)")
    println("[scan magnetic] candidates_output=$(stats.candidates_output)")
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
    println("  julia --project=. scripts/models/run_unified_scan.jl scan tmu [--key=value ...]")
    println("  julia --project=. scripts/models/run_unified_scan.jl scan trho [--key=value ...]")
    println("  julia --project=. scripts/models/run_unified_scan.jl scan magnetic [--key=value ...]")
    println("  julia --project=. scripts/models/run_unified_scan.jl workflow phase [--key=value ...]")
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
        elseif args[2] == "magnetic"
            return _run_scan_magnetic(args[3:end])
        end
        throw(ArgumentError("unknown scan subcommand: $(args[2]); expected: tmu|trho|magnetic"))
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
