module TransportSamplingScanConfig

using JSON3, SHA

const SOURCE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function reject_retired_sampling_flags(args)
    retired=("--sigma-grid-n","--tau-threshold-subtraction","--tau-no-threshold-subtraction",
        "--no-tau-threshold-subtraction","--tau-asym-window","--tau-asym-fit-min-points",
        "--tau-asym-extra-points","--tau-interpolation-mode","--tr-p-max")
    for arg in args
        flag=first(split(arg,'=';limit=2))
        flag in retired && throw(ArgumentError("$flag was retired; configure the thermal-weighted sampling profile instead"))
    end
    return nothing
end

sampling_diagnostics_path(opts) = splitext(opts.output)[1] * "_sampling.csv"
sampling_contract_path(opts) = splitext(opts.output)[1] * "_sampling_contract.json"

function calculation_identity(opts)
    profile = profile_for_options(opts)
    profile === nothing && return nothing
    files = String[]
    for base in ("src", "scripts/relaxtime", "config/models/pnjl", "config/physics")
        for (dir, _, names) in walkdir(joinpath(SOURCE_ROOT, base))
            for name in names
                endswith(name,".jl") || endswith(name,".toml") || continue
                push!(files,joinpath(dir,name))
            end
        end
    end
    append!(files,[joinpath(SOURCE_ROOT,"Project.toml"),joinpath(SOURCE_ROOT,"Manifest.toml")])
    sources = [(path=replace(relpath(path,SOURCE_ROOT),'\\'=>'/'),
        sha256=bytes2hex(sha256(read(path)))) for path in sort!(unique!(files))]
    params = model_for_options(opts).params
    model_parameters = NamedTuple{fieldnames(typeof(params))}(Tuple(getfield(params,k) for k in fieldnames(typeof(params))))
    return (schema_version=1,profile_id=profile.id,profile_sha256=profile.source_sha256,
        propagator_xi_policy=String(opts.propagator_xi_policy),compute_bulk=opts.compute_bulk,
        model_parameters=model_parameters,sources=sources)
end

"""Bind a result CSV to its numerical configuration before writing any points."""
function prepare_sampling_output!(opts)
    path = sampling_contract_path(opts)
    profile = profile_for_options(opts)
    if profile === nothing
        isfile(path) && !opts.overwrite &&
            throw(ArgumentError("result belongs to a sampling profile; use that profile or a new output path"))
        return nothing
    end
    identity = calculation_identity(opts)
    digest = bytes2hex(sha256(JSON3.write(identity)))
    occupied = isfile(opts.output) && filesize(opts.output)>0
    if !opts.overwrite
        occupied && !opts.resume && throw(ArgumentError("result already exists; use --resume or a new output path"))
        if isfile(path)
            saved = JSON3.read(read(path,String))
            saved["identity_sha256"] == digest ||
                throw(ArgumentError("sampling configuration or numerical source changed; use a new output path"))
        elseif occupied || isfile(sampling_diagnostics_path(opts))
            throw(ArgumentError("existing result has no sampling identity; cannot resume it with this profile"))
        end
    elseif isfile(sampling_diagnostics_path(opts))
        rm(sampling_diagnostics_path(opts))
    end
    mkpath(dirname(path))
    open(path,"w") do io
        JSON3.pretty(io,(identity_sha256=digest,identity=identity))
    end
    return path
end

function profile_for_options(opts)
    hasproperty(opts,:sampling_profile) && opts.sampling_profile !== nothing || return nothing
    return Main.Models.TransportSamplingProfiles.load_profile(opts.sampling_profile)
end

function model_for_options(opts)
    hasproperty(opts,:model) && return opts.model
    profile = profile_for_options(opts)
    return profile === nothing ? Main.PNJL_MODEL : Main.Models.TransportSamplingProfiles.profile_model(profile)
end

function residual_norm_max(opts)
    hasproperty(opts,:residual_norm_max) && return opts.residual_norm_max
    profile = profile_for_options(opts)
    return profile === nothing ? 1e-4 : profile.background.residual_norm_max
end

function configure_profile!(opts,args;phase_guided=false)
    opts[:sampling_profile] === nothing && return opts
    p = Main.Models.TransportSamplingProfiles.load_profile(opts[:sampling_profile])
    settings = [(:p_num,p.background.p_nodes,"--p-num"),(:t_num,p.background.cos_nodes,"--t-num"),
        (:tau_p_nodes,p.tau_integral.p_nodes,"--tau-p-nodes"),
        (:tau_angle_nodes,p.tau_integral.cos_nodes,"--tau-angle-nodes"),
        (:tau_phi_nodes,p.tau_integral.phi_nodes,"--tau-phi-nodes"),
        (phase_guided ? :tau_n_sigma_points : :tau_n_sigma,p.tau_integral.t_nodes,"--tau-n-sigma"),
        (:sigma_cache_policy,:thermal_weighted,"--sigma-cache-policy")]
    phase_guided || append!(settings,[(:tr_p_nodes,p.terminal_integral.p_nodes,"--tr-p-nodes"),
        (:integration_mode,:semi_infinite,"--mode")])
    phase_guided && push!(settings,(:phase_anchor_policy,:direct_coexistence,"--phase-anchor-policy"))
    for (key,value,flag) in settings
        flag in args && opts[key] != value && throw(ArgumentError("$flag conflicts with sampling profile $(p.id)"))
        opts[key] = value
    end
    opts[:sampling_profile] = p.source_path
    return opts
end

function record_profile!(cfg,opts)
    profile = profile_for_options(opts)
    profile === nothing && return cfg
    cfg["sampling_profile"] = Main.Models.TransportSamplingProfiles.effective_config(profile)
    cfg["terminal_integration_mode"] = "semi_infinite"
    cfg["phase_reference_usage"] = "seed hint; stable multiseed selection uses the sampling-profile background"
    cfg["tr_p_max_fm"] = nothing
    return cfg
end

end # module TransportSamplingScanConfig
