#!/usr/bin/env julia

isdefined(Main, :ChargedGBUContourScan) || Base.include(Main,
    joinpath(@__DIR__, "run_charged_gbu_contour_scan.jl"))

"""Compare published chemical freezeout profiles using the retained q=0 route.

This analysis runner reuses the contour scanner's q0_lambda_reference kernel and
the unified quark-only BQS equilibrium. It does not replace the finite-q
production default. Each requested energy is solved directly, never interpolated
from an old contour. Numerical output is separate from plotting and scoring.
"""
module ChargedGBUFreezeoutComparison

using CSV, Dates, JSON3, SHA, TOML, Printf
using Main.Models

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const Scan = Main.ChargedGBUContourScan
const Workflow = Scan.Workflow
const FP = Models.FreezeoutProfiles
const ROUTE = "q0_lambda_reference"
const CHANNELS = (:pi_plus, :K_plus)
const DEFAULT_PROFILES = ["cleymans2006", "andronic2018"]
const PARAMETER_POLICIES = ("published_central", "literature_envelope_exploration")

function validate_parameter_policy(policy::AbstractString)
    policy in PARAMETER_POLICIES || throw(ArgumentError(
        "parameter_policy must be published_central or literature_envelope_exploration"))
    return String(policy)
end
const DEFAULT_ENERGIES = sort!(unique!(vcat(collect(2.5:0.25:12.0),
    [2.35, 2.4, 2.7, 3.32, 3.83, 4.29, 4.87, 6.27, 7.62, 7.7, 8.77, 11.5, 12.32, 14.5,
     17.27, 19.6, 22., 27., 32., 39., 50., 62.4, 80., 100., 130., 160.,
     200., 300., 500., 1000., 2000., 2760., 5020.])))

function numerical_settings(resolution::AbstractString)
    base = Scan.screening_settings(; route=ROUTE)
    resolution == "screening" && return base
    resolution == "standard" && return merge(base,
        (mesh=128, cut_nodes=64, tail_nodes=64, omega_nodes=96, q_nodes=16))
    resolution == "fine" && return merge(base,
        (mesh=256, cut_nodes=96, tail_nodes=96, omega_nodes=128, q_nodes=32))
    throw(ArgumentError("resolution must be screening, standard or fine"))
end

function comparison_points(profile_names, energies; profile_dir=FP.freezeout_profile_dir())
    names = String.(profile_names)
    !isempty(names) && length(unique(names)) == length(names) ||
        throw(ArgumentError("profiles must be nonempty and unique"))
    es = Workflow.energy_grid(energies) # descending energy for continuation
    profiles = [FP.load_freezeout_profile(; profile=name, dir=profile_dir) for name in names]
    rows = NamedTuple[]
    for (name, profile) in zip(names, profiles)
        profile.profile_name == name || throw(ArgumentError("profile_name differs from requested file: $(name)"))
        points = FP.build_freezeout_scan_points(es; profile=profile, traversal=:sqrts_descending)
        all(p -> isfinite(p.T_MeV) && p.T_MeV > 0, points) ||
            throw(ArgumentError("nonpositive temperature on requested $(name) chemical freezeout path"))
        for (index, pt) in enumerate(points)
            push!(rows, (profile_name=name, family=profile.family, source_tag=profile.source_tag,
                energy_index=index, point=pt))
        end
    end
    return rows
end

function source_hashes(config_path, names, profile_dir)
    hashes = Scan.source_hashes(config_path)
    paths = [@__FILE__, joinpath(ROOT, "src", "models", "scans", "FreezeoutProfiles.jl"),
        joinpath(ROOT, "src", "config", "ConfigLoader.jl"), joinpath(ROOT, "Project.toml")]
    isfile(joinpath(ROOT, "Manifest.toml")) && push!(paths, joinpath(ROOT, "Manifest.toml"))
    append!(paths, [joinpath(profile_dir, name * ".toml") for name in names])
    for path in paths
        hashes[replace(relpath(abspath(path), ROOT), '\\'=>'/')] = Workflow.hashfile(path)
    end
    return hashes
end

function summarize(record)
    ch = record.channels
    getchannel(name) = hasproperty(ch, name) ? getproperty(ch, name) : nothing
    pi_record = getchannel(:pi_plus); k_record = getchannel(:K_plus)
    number(r, key) = r === nothing ? NaN : Scan._as_float(getproperty(r, key))
    countwarnings(r) = r === nothing ? 0 : Int(r.endpoint_warning_shells)
    return (profile_name=String(record.profile_name), family=String(record.family),
        source_tag=String(record.source_tag), sqrt_s_NN_GeV=Float64(record.sqrt_s_NN_GeV),
        T_MeV=Float64(record.T_MeV), muB_MeV=Float64(record.muB_MeV),
        Kplus_over_pi_plus=Scan._as_float(record.ratios.Kplus_over_pi_plus),
        plus_passed=Bool(record.ratios.plus_passed), status=String(record.status),
        pi_plus_density_inv_fm3=number(pi_record, :density_inv_fm3),
        K_plus_density_inv_fm3=number(k_record, :density_inv_fm3),
        pi_plus_endpoint_warning_shells=countwarnings(pi_record),
        K_plus_endpoint_warning_shells=countwarnings(k_record),
        background_residual=record.background === nothing ? NaN : Float64(record.background.bg.residual),
        failure_reason=String(record.failure_reason),
        pi_plus_failure_reason=pi_record === nothing ? "" : String(pi_record.reason),
        K_plus_failure_reason=k_record === nothing ? "" : String(k_record.reason),
        density_route=ROUTE, endpoint_policy="finite_window", production_authorized=false)
end

function run_scan(; output, profiles=DEFAULT_PROFILES, energies=DEFAULT_ENERGIES,
        resolution="standard", resume=false, config=Workflow.DEFAULT_CONFIG,
        profile_dir=FP.freezeout_profile_dir(), parameter_policy="published_central")
    policy = validate_parameter_policy(parameter_policy)
    profile_dir = abspath(profile_dir)
    points = comparison_points(profiles, energies; profile_dir=profile_dir)
    settings = numerical_settings(resolution)
    config_path = abspath(config); c = Workflow.validate_config(TOML.parsefile(config_path))
    output = abspath(output)
    hashes = source_hashes(config_path, profiles, profile_dir)
    inputs = (schema="charged_gbu_freezeout_comparison_v1", points=points, settings=settings,
        config=c, density_route=ROUTE, source_hashes=hashes, parameter_policy=policy)
    identity = bytes2hex(sha256(JSON3.write(inputs)))
    statepath = joinpath(output, "run.json")
    if ispath(output)
        resume || throw(ArgumentError("output exists; use --resume or a new output directory"))
        previous = Workflow.readjson(statepath)
        String(previous.scan_identity) == identity || throw(ArgumentError("resume identity mismatch"))
    else
        resume && throw(ArgumentError("cannot resume an absent output directory"))
        mkpath(output)
        # The full base commit supplies unchanged transitive dependencies; retain
        # exact bytes of the numerical inputs, including this uncommitted runner.
        for path in keys(hashes)
            isabspath(path) && throw(ArgumentError("snapshot key must be relative: $(path)"))
            target = startswith(path, "../") ? joinpath(output, "source_snapshot", "external_profiles", basename(path)) :
                joinpath(output, "source_snapshot", path)
            mkpath(dirname(target)); cp(joinpath(ROOT, path), target)
        end
        Workflow.writejson(statepath, merge(inputs, (scan_identity=identity,
            git_head=readchomp(`git -C $ROOT rev-parse HEAD`),
            git_status=readchomp(`git -C $ROOT status --porcelain`),
            generated_at_utc=string(Dates.now(Dates.UTC)), julia_version=string(VERSION),
            coordinate_contract=Scan.coordinate_contract(ROUTE),
            background_contract=(type="FixedMuBConservedCharges", nS_inv_fm3=0., nQ_over_nB=.4,
                xi=0., meson_feedback=false),
            diagnostic_only=true, production_authorized=false)))
    end
    model = Models.create_model(:PNJL)
    seed = nothing; previous_profile = ""; paths = String[]
    for item in points
        pt = item.point
        if previous_profile != item.profile_name
            seed = nothing; previous_profile = item.profile_name
        end
        path = joinpath(output, "points", item.profile_name,
            @sprintf("energy_%03d.json", item.energy_index))
        mkpath(dirname(path))
        if resume && isfile(path)
            old = Workflow.readjson(path)
            String(old.scan_identity) == identity || throw(ArgumentError("point identity mismatch: $(path)"))
            seed = old.background === nothing ? nothing : Float64.(old.background.seed)
            push!(paths, path); continue
        end
        started = Base.time_ns()
        solved = nothing; records = Dict{Symbol,Any}(); reason = ""
        try
            solved = Workflow.background(model, pt, c; seed=seed)
            seed = solved.seed
            for channel in CHANNELS
                records[channel] = Scan._channel_record(solved.bg, channel, settings; route=ROUTE)
            end
        catch err
            err isa InterruptException && rethrow()
            reason = sprint(showerror, err); seed = nothing
        end
        ratios = Scan._ratios(CHANNELS, records)
        status = solved === nothing ? "background_failed" :
            ratios.plus_passed ? "screened" : "gate_failed"
        record = (schema="charged_gbu_freezeout_point_v1", scan_identity=identity,
            profile_name=item.profile_name, family=item.family, source_tag=item.source_tag,
            sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV, T_MeV=pt.T_MeV, muB_MeV=pt.muB_MeV,
            energy_index=item.energy_index, background=solved, channels=records, ratios=ratios,
            status=status, failure_reason=reason, elapsed_s=Scan._elapsed_seconds(started))
        Workflow.writejson(path, record); push!(paths, path)
        @printf("%s sqrt(s)=%.6g T=%.4f muB=%.4f K/pi=%.7g %s (%.2f s)\n",
            item.profile_name, pt.sqrt_s_NN_GeV, pt.T_MeV, pt.muB_MeV,
            ratios.Kplus_over_pi_plus, status, record.elapsed_s)
        flush(stdout)
    end
    rows = [summarize(Workflow.readjson(path)) for path in paths]
    sort!(rows; by=x -> (x.profile_name, x.sqrt_s_NN_GeV))
    csvpath = joinpath(output, "freezeout_ratios.csv")
    CSV.write(csvpath, rows)
    manifest = (schema="charged_gbu_freezeout_comparison_manifest_v1", scan_identity=identity,
        run_file="run.json", run_sha256=Workflow.hashfile(statepath),
        generated_at_utc=string(Dates.now(Dates.UTC)), resolution=resolution,
        profiles=profiles, point_count=length(rows), successful_points=count(r -> r.plus_passed, rows),
        warning_points=count(r -> r.pi_plus_endpoint_warning_shells+r.K_plus_endpoint_warning_shells > 0, rows),
        density_route=ROUTE, settings=settings, diagnostic_only=true, production_authorized=false,
        parameter_policy=policy,
        gate_scope="q0 reference geometry and shell finiteness; full finite-q production gates not applicable",
        warning_policy="retain q0 static-imaginary diagnostics; failed points remain missing",
        outputs=Dict(replace(relpath(p, output), '\\'=>'/')=>Workflow.hashfile(p) for p in vcat(paths, [csvpath])))
    Workflow.writejson(joinpath(output, "manifest.json"), manifest)
    return manifest
end

function parse_args(args)
    opts = Dict{Symbol,Any}(:profiles=>DEFAULT_PROFILES, :energies=>DEFAULT_ENERGIES,
        :resolution=>"standard", :resume=>false, :parameter_policy=>"published_central")
    i = 1
    while i <= length(args)
        arg = args[i]
        arg in ("--help", "-h") && return nothing
        if arg == "--resume"
            opts[:resume] = true
        elseif arg in ("--output", "--profiles", "--sqrts-list", "--resolution", "--config", "--profile-dir", "--parameter-policy")
            i < length(args) || throw(ArgumentError("missing value for $(arg)")); i += 1
            key = arg == "--sqrts-list" ? :energies : Symbol(replace(arg[3:end], '-'=>'_'))
            opts[key] = key == :energies ? parse.(Float64, split(args[i], ',')) :
                key == :profiles ? String.(split(args[i], ',')) : args[i]
        else
            throw(ArgumentError("unknown option: $(arg)"))
        end
        i += 1
    end
    haskey(opts, :output) || throw(ArgumentError("--output is required"))
    validate_parameter_policy(opts[:parameter_policy])
    return opts
end

function main(args=ARGS)
    opts = parse_args(args)
    if opts === nothing
        println("Usage: julia --project=. scripts/analysis/relaxtime/run_charged_gbu_freezeout_comparison.jl --output DIR [--profiles cleymans2006,andronic2018] [--sqrts-list 7.7,19.6,200] [--resolution screening|standard|fine] [--profile-dir DIR] [--parameter-policy published_central|literature_envelope_exploration] [--resume]")
        return
    end
    run_scan(; opts...)
end

end # module

abspath(PROGRAM_FILE) == abspath(@__FILE__) && ChargedGBUFreezeoutComparison.main()
