#!/usr/bin/env julia

"""Run a resumable diagnostic T-μB charged-GBU screening shard.

This entrypoint evaluates fixed quark-only BQS backgrounds with the current
infinite-thermal kernel. It intentionally uses a screening resolution and does
not claim full production acceptance. T rows can be distributed across GitHub
Actions jobs; each row continues in μB order so the previous equilibrium seed
is retained. Point JSON files are atomic and carry the immutable scan identity.
"""

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
pushfirst!(LOAD_PATH, ROOT)
isdefined(Main, :Models) || Base.include(Main, joinpath(ROOT, "src", "models", "Models.jl"))
isdefined(Main, :ChargedGBUResearchWorkflow) || Base.include(Main,
    joinpath(ROOT, "src", "models", "workflow_apps", "ChargedGBUResearchWorkflow.jl"))

module ChargedGBUContourScan

using CSV, Dates, JSON3, SHA, TOML, Printf
using Main.Models

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const Workflow = Main.ChargedGBUResearchWorkflow
include("charged_gbu_q0_reference.jl")
const Reference = ChargedGBUQ0Reference
const CHANNELS = Workflow.R.CHANNELS
const DEFAULT_CONFIG = Workflow.DEFAULT_CONFIG
const DEFAULT_T_GRID = "40:220:10"
const DEFAULT_MUB_GRID = "0:800:50"
const DENSITY_ROUTES = ("direct_finite_q", "q0_lambda_reference")

function density_route(value)
    value in DENSITY_ROUTES || throw(ArgumentError("unsupported density_route: $(value)"))
    return String(value)
end

coordinate_contract(route) = (internal_frequency="lambda=omega+(mu1-mu2)",
    bose_frequency="external omega", extrapolation=route == "q0_lambda_reference" ?
        "lambda0=sqrt(lambda^2-q^2) for lambda>=q; phase=0 otherwise" : "none; direct finite-q loop")

_float(value, label) = begin
    result = try parse(Float64, String(value)) catch
        throw(ArgumentError("$(label) must be a finite number, got $(value)"))
    end
    isfinite(result) || throw(ArgumentError("$(label) must be finite, got $(value)"))
    result
end
_elapsed_seconds(start_ns::UInt64) = (Base.time_ns() - start_ns) / 1.0e9

function parse_grid(spec::AbstractString; label::AbstractString="grid")
    fields = split(strip(spec), ':')
    length(fields) == 3 || throw(ArgumentError("$(label) must use start:stop:step, got $(spec)"))
    first_value, last_value, step = (_float(fields[i], label) for i in 1:3)
    step > 0 && last_value >= first_value || throw(ArgumentError("invalid $(label) bounds/step"))
    count = floor(Int, (last_value - first_value) / step + 1e-10) + 1
    endpoint = first_value + (count - 1) * step
    isapprox(endpoint, last_value; atol=1e-8, rtol=1e-10) ||
        throw(ArgumentError("$(label) does not land on its stop"))
    values = [first_value + (i - 1) * step for i in 1:count]
    values[end] = last_value
    return values
end

function parse_channels(spec::AbstractString)
    names = [strip(name) for name in split(spec, ',') if !isempty(strip(name))]
    isempty(names) && throw(ArgumentError("channels must not be empty"))
    symbols = Symbol.(names)
    all(in(CHANNELS), symbols) || throw(ArgumentError("unsupported charged channel"))
    unique(symbols) == symbols || throw(ArgumentError("channels must be unique"))
    return symbols
end

function shard_indices(total::Integer, shard_index::Integer, shard_count::Integer)
    total >= 0 && shard_count >= 1 && 0 <= shard_index < shard_count ||
        throw(ArgumentError("invalid shard specification"))
    return [i for i in 1:total if mod(i - 1, shard_count) == shard_index]
end

function point(T_MeV::Real, muB_MeV::Real)
    T = Float64(T_MeV); muB = Float64(muB_MeV)
    isfinite(T) && T > 0 && isfinite(muB) && muB >= 0 || throw(ArgumentError("invalid point"))
    h = Main.Constants_PNJL.ħc_MeV_fm
    return (sqrt_s_NN_GeV=NaN, muB_GeV=muB / 1000, T_GeV=T / 1000,
        muB_MeV=muB, T_MeV=T, muq_MeV=muB / 3, muB_fm=muB / h, T_fm=T / h)
end

function screening_settings()
    # Deliberately diagnostic resolution; accepted production settings remain in
    # charged_gbu_infinite_v1.toml and are never lowered by this scanner.
    return (mesh=64, cut_nodes=32, tail_nodes=32, omega_nodes=64, q_nodes=8, qmax=8.0)
end

_jsonsafe(x::AbstractFloat) = isfinite(x) ? x : nothing
_jsonsafe(x::NamedTuple) = Dict(String(k) => _jsonsafe(v) for (k, v) in pairs(x))
_jsonsafe(x::AbstractDict) = Dict(String(k) => _jsonsafe(v) for (k, v) in pairs(x))
_jsonsafe(x::AbstractArray) = [_jsonsafe(v) for v in x]
_jsonsafe(x) = x
_as_float(x) = x === nothing ? NaN : Float64(x)

function _write_atomic(path, payload)
    mkpath(dirname(path)); partial = path * ".partial"
    open(partial, "w") do io
        JSON3.write(io, _jsonsafe(payload)); write(io, '\n')
    end
    mv(partial, path; force=true)
end

_read_json(path) = JSON3.read(read(path, String))

function _number_tag(value::Real)
    replace(@sprintf("%.6f", Float64(value)), '.' => 'p', '-' => 'm')
end

function point_path(output, row_index, col_index, T_MeV, muB_MeV)
    joinpath(output, "points", "point_T$(row_index)_mu$(col_index)_$(_number_tag(T_MeV))_$(_number_tag(muB_MeV)).json")
end

function source_hashes(config_path)
    paths = [joinpath(ROOT, "src", "models", "Models.jl"),
        joinpath(ROOT, "src", "models", "solver", "runtime", "ConstraintSolverFixedMuBConservedCharges.jl"),
        joinpath(ROOT, "src", "models", "workflow_apps", "ChargedGBUResearchWorkflow.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_qgate.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_thermal.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_profile.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "causal_gbu_infinite_yield.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "charged_gbu_q0_reference.jl"),
        joinpath(ROOT, "scripts", "analysis", "relaxtime", "run_charged_gbu_contour_scan.jl"), config_path]
    return Dict(replace(relpath(p, ROOT), '\\' => '/') => bytes2hex(sha256(read(p))) for p in paths)
end

function _scan_identity(config, T_values, muB_values, channels, settings, shard_index, shard_count, hashes;
        route="direct_finite_q")
    payload = JSON3.write((config=config, T_values=T_values, muB_values=muB_values,
        channels=String.(channels), settings=settings, shard_index=shard_index,
        density_route=density_route(route), coordinate_contract=coordinate_contract(route),
        shard_count=shard_count, source_hashes=hashes, schema="charged_gbu_contour_scan_v2"))
    bytes2hex(sha256(payload))
end

function _channel_record(bg, channel, settings; route="direct_finite_q")
    density_route(route)
    started = Base.time_ns()
    try
        q_values, q_weights = Workflow.R.gauleg(0.0, settings.qmax, settings.q_nodes)
        ref = if route == "q0_lambda_reference"
            kernel0 = Workflow.I.kernel(bg, channel, 0.; cut_nodes=settings.cut_nodes,
                split_inv_fm=max(36.0, 28.0 * bg.T + 2.0))
            Reference.reference(Workflow.P.profile(kernel0; mesh=settings.mesh, tail_nodes=settings.tail_nodes))
        else
            nothing
        end
        rows = NamedTuple[]
        for (q, weight) in zip(q_values, q_weights)
            shell = if route == "q0_lambda_reference"
                Reference.shell(ref, q; nodes=settings.omega_nodes)
            else
                kernel = Workflow.I.kernel(bg, channel, q; cut_nodes=settings.cut_nodes,
                    split_inv_fm=max(36.0, q + 28.0 * bg.T + 2.0))
                profile = Workflow.P.profile(kernel; mesh=settings.mesh, tail_nodes=settings.tail_nodes)
                Workflow.Y.shell(profile; nodes=settings.omega_nodes)
            end
            push!(rows, (q_inv_fm=q, weight=weight, density=shell.density,
                root_count=shell.root_count, passed=isfinite(shell.density)))
        end
        density = sum(row.weight * row.density for row in rows)
        return (status=isfinite(density) && all(row.passed for row in rows) ? "screened" : "nonfinite",
            passed=isfinite(density) && all(row.passed for row in rows), density_inv_fm3=density,
            elapsed_s=_elapsed_seconds(started), rows=rows,
            density_route=route,
            gate_scope="finite_shell_output_only;full_production_gates_omitted", reason="")
    catch err
        err isa InterruptException && rethrow()
        return (status="evaluation_failed", passed=false, density_inv_fm3=NaN,
            density_route=route,
            elapsed_s=_elapsed_seconds(started), rows=NamedTuple[], gate_scope="screening", reason=sprint(showerror, err))
    end
end

function _ratios(channels, records)
    get_record(ch) = haskey(records, ch) ? records[ch] : nothing
    function one(num, den)
        a, b = get_record(num), get_record(den)
        a !== nothing && b !== nothing && a.passed && b.passed && isfinite(a.density_inv_fm3) &&
            isfinite(b.density_inv_fm3) && b.density_inv_fm3 > 0 || return (value=NaN, passed=false)
        return (value=a.density_inv_fm3 / b.density_inv_fm3, passed=true)
    end
    plus, minus = one(:K_plus, :pi_plus), one(:K_minus, :pi_minus)
    return (Kplus_over_pi_plus=plus.value, Kminus_over_pi_minus=minus.value,
        plus_passed=plus.passed, minus_passed=minus.passed)
end

function _point_record(identity, row_index, col_index, T_MeV, muB_MeV, bg, seed, records, channels, elapsed_s, failure_reason)
    status = bg === nothing ? "background_failed" : all(records[ch].passed for ch in channels) ? "screened" : "gate_failed"
    return (schema="charged_gbu_contour_point_v2", scan_identity=identity, generated_at_utc=string(Dates.now(Dates.UTC)),
        row_index=row_index, col_index=col_index, T_MeV=Float64(T_MeV), muB_MeV=Float64(muB_MeV), status=status,
        failure_reason=failure_reason,
        background=bg === nothing ? nothing : (residual=bg.residual, elapsed_s=elapsed_s, seed=seed), channels=records,
        ratios=bg === nothing ? (Kplus_over_pi_plus=NaN, Kminus_over_pi_minus=NaN, plus_passed=false, minus_passed=false) :
            _ratios(channels, records))
end

function _summary_row(data)
    row = Dict{Symbol,Any}(:row_index=>Int(data.row_index), :col_index=>Int(data.col_index),
        :T_MeV=>Float64(data.T_MeV), :muB_MeV=>Float64(data.muB_MeV), :status=>String(data.status),
        :background_residual=>data.background === nothing ? NaN : _as_float(data.background.residual),
        :background_elapsed_s=>data.background === nothing ? NaN : _as_float(data.background.elapsed_s),
        :Kplus_over_pi_plus=>_as_float(data.ratios.Kplus_over_pi_plus), :Kminus_over_pi_minus=>_as_float(data.ratios.Kminus_over_pi_minus),
        :plus_passed=>Bool(data.ratios.plus_passed), :minus_passed=>Bool(data.ratios.minus_passed))
    for ch in CHANNELS
        record = hasproperty(data.channels, ch) ? getproperty(data.channels, ch) : nothing
        prefix = String(ch); row[Symbol(prefix * "_density_inv_fm3")] = record === nothing ? NaN : _as_float(record.density_inv_fm3)
        row[Symbol(prefix * "_passed")] = record !== nothing && Bool(record.passed)
        row[Symbol(prefix * "_elapsed_s")] = record === nothing ? NaN : _as_float(record.elapsed_s)
    end
    return (; row...)
end

function parse_args(args=ARGS)
    options = Dict{Symbol,Any}(:t_grid=>DEFAULT_T_GRID, :muB_grid=>DEFAULT_MUB_GRID,
        :channels=>join(String.(CHANNELS), ","), :output=>nothing, :shard_index=>0, :shard_count=>1,
        :resume=>false, :config=>DEFAULT_CONFIG, :density_route=>"direct_finite_q")
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--help", "-h")
            return nothing
        elseif arg == "--resume"
            options[:resume] = true
        elseif arg in ("--t-grid", "--muB-grid", "--channels", "--output", "--shard-index", "--shard-count", "--config", "--density-route")
            i < length(args) || throw(ArgumentError("missing value for $(arg)")); i += 1
            key = arg == "--t-grid" ? :t_grid : arg == "--muB-grid" ? :muB_grid : arg == "--channels" ? :channels :
                arg == "--output" ? :output : arg == "--shard-index" ? :shard_index : arg == "--shard-count" ? :shard_count :
                arg == "--density-route" ? :density_route : :config
            options[key] = key in (:shard_index, :shard_count) ? parse(Int, args[i]) : args[i]
        else
            throw(ArgumentError("unknown option $(arg)"))
        end
        i += 1
    end
    options[:output] === nothing && throw(ArgumentError("--output is required"))
    density_route(options[:density_route])
    return options
end

print_help() = println("Usage: julia --project=. scripts/analysis/relaxtime/run_charged_gbu_contour_scan.jl --output DIR [--t-grid 40:220:10] [--muB-grid 0:800:50] [--channels pi_plus,K_plus] [--shard-index 0] [--shard-count 4] [--density-route direct_finite_q|q0_lambda_reference] [--resume]")

function run_scan(options)
    config_path = abspath(String(options[:config]))
    c = Workflow.validate_config(TOML.parsefile(config_path))
    T_values = parse_grid(options[:t_grid]; label="T grid"); muB_values = parse_grid(options[:muB_grid]; label="muB grid")
    channels = parse_channels(options[:channels]); settings = screening_settings(); output = abspath(String(options[:output]))
    route = density_route(options[:density_route])
    options[:shard_count] <= length(T_values) || throw(ArgumentError("shard_count cannot exceed T-row count"))
    ispath(output) && !options[:resume] && throw(ArgumentError("output exists; use --resume")); mkpath(output)
    hashes = source_hashes(config_path); identity = _scan_identity(c, T_values, muB_values, channels, settings, options[:shard_index], options[:shard_count], hashes; route=route)
    manifest_path = joinpath(output, "manifest.json")
    if options[:resume] && isfile(manifest_path)
        previous = _read_json(manifest_path)
        hasproperty(previous, :scan_identity) && String(previous.scan_identity) == identity ||
            throw(ArgumentError("resume manifest identity mismatch: $(manifest_path)"))
    end
    selected = shard_indices(length(T_values), options[:shard_index], options[:shard_count]); model = Models.create_model(:PNJL); point_paths = String[]
    for row_index in selected
        seed = nothing
        for (col_index, muB_MeV) in enumerate(muB_values)
            path = point_path(output, row_index, col_index, T_values[row_index], muB_MeV)
            if options[:resume] && isfile(path)
                old = try _read_json(path) catch; nothing end
                if old !== nothing && hasproperty(old, :scan_identity) && String(old.scan_identity) == identity
                    if hasproperty(old, :background) && old.background !== nothing && hasproperty(old.background, :seed)
                        seed = Float64.(old.background.seed)
                    else
                        seed = nothing
                    end
                    push!(point_paths, path); continue
                end
                old === nothing && throw(ArgumentError("resume point JSON is unreadable: $(path)"))
                throw(ArgumentError("resume point identity mismatch: $(path)"))
            end
            bg, records, elapsed_s, background_seed, failure_reason = nothing, Dict{Symbol,Any}(), NaN, nothing, ""
            try
                started = Base.time_ns(); solved = Workflow.background(model, point(T_values[row_index], muB_MeV), c; seed=seed)
                elapsed_s = _elapsed_seconds(started); bg = solved.bg; seed = solved.seed; background_seed = solved.seed
                for channel in channels; records[channel] = _channel_record(bg, channel, settings; route=route); end
            catch err
                err isa InterruptException && rethrow(); seed = nothing; failure_reason = sprint(showerror, err)
                @warn "charged GBU contour point failed" row_index col_index T_MeV=T_values[row_index] muB_MeV failure_reason
            end
            record = _point_record(identity, row_index, col_index, T_values[row_index], muB_MeV, bg, background_seed, records, channels, elapsed_s, failure_reason)
            _write_atomic(path, record); push!(point_paths, path)
        end
    end
    rows = [_summary_row(_read_json(path)) for path in point_paths if isfile(path)]; sort!(rows; by=x->(x.row_index, x.col_index))
    CSV.write(joinpath(output, "contour_points.csv"), rows)
    manifest = (schema="charged_gbu_contour_scan_v2", scan_identity=identity, generated_at_utc=string(Dates.now(Dates.UTC)),
        git_head=readchomp(`git -C $ROOT rev-parse HEAD`), config=replace(relpath(config_path, ROOT), '\\'=>'/'),
        source_hashes=hashes, T_grid=T_values, muB_grid=muB_values, channels=String.(channels), settings=settings,
        density_route=route, coordinate_contract=coordinate_contract(route),
        shard_index=options[:shard_index], shard_count=options[:shard_count], point_count=length(rows),
        successful_points=count(row->row.status == "screened", rows), production_default=false, diagnostic_only=true)
    _write_atomic(joinpath(output, "manifest.json"), manifest); return manifest
end

function main(args=ARGS)
    options = parse_args(args); options === nothing && return print_help(); run_scan(options)
end

end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && ChargedGBUContourScan.main()
