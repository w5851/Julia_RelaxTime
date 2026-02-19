#!/usr/bin/env julia

using TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS_CFG_ROOT = joinpath(ROOT, "config", "models")

const DOMAINS = (
    :njl,
    :pnjl,
    :rpnjl,
)

const REQUIRED_KEYS = Dict(
    :njl => [
        "model.label",
        "model.N_color",
        "model.N_flavor",
        "model.Lambda_MeV",
        "model.G_over_Lambda2",
        "model.K_over_Lambda5",
        "model.m_ud0_MeV",
        "model.m_s0_MeV",
    ],
    :pnjl => [
        "model.N_color",
        "model.N_flavor",
        "model.Lambda_MeV",
        "model.G_over_Lambda2",
        "model.K_over_Lambda5",
        "model.m_ud0_MeV",
        "model.m_s0_MeV",
        "polyakov.T0_MeV",
        "polyakov.a0",
        "polyakov.a1",
        "polyakov.a2",
        "polyakov.b3",
        "polyakov.b4",
    ],
    :rpnjl => [
        "rpnjl.Lambda_MeV",
        "rpnjl.G_over_Lambda2",
        "rpnjl.K_over_Lambda5",
        "rpnjl.m_ud0_MeV",
        "rpnjl.m_s0_MeV",
        "rpnjl.g1_MeV_inv8",
        "rpnjl.g2_MeV_inv8",
        "rpnjl.kappa",
        "rpnjl.T0_MeV",
        "rpnjl.a0",
        "rpnjl.a1",
        "rpnjl.a2",
        "rpnjl.b3",
        "rpnjl.b4",
    ],
)

function flatten_keys(data::Dict{String, Any}; prefix::String="")
    out = Set{String}()
    for (k, v) in data
        path = isempty(prefix) ? k : string(prefix, ".", k)
        if v isa Dict{String, Any}
            union!(out, flatten_keys(v; prefix=path))
        else
            push!(out, path)
        end
    end
    return out
end

function deep_merge(base::Dict{String, Any}, override::Dict{String, Any})
    result = deepcopy(base)
    for (k, v) in override
        if haskey(result, k) && result[k] isa Dict{String, Any} && v isa Dict{String, Any}
            result[k] = deep_merge(result[k], v)
        else
            result[k] = v
        end
    end
    return result
end

function model_dir(domain::Symbol)
    return joinpath(MODELS_CFG_ROOT, String(domain))
end

function parse_profile(path::String)
    return isfile(path) ? TOML.parsefile(path) : Dict{String, Any}()
end

function profile_path(domain::Symbol, profile::String)
    return joinpath(model_dir(domain), string(profile, ".toml"))
end

function summarize_overrides(default_cfg::Dict{String, Any}, override_cfg::Dict{String, Any})
    default_keys = flatten_keys(default_cfg)
    override_keys = flatten_keys(override_cfg)
    unknown_override = sort(collect(setdiff(override_keys, default_keys)))
    return (override_keys=sort(collect(override_keys)), unknown_override=unknown_override)
end

function assert_required_keys!(violations::Vector{String}, domain::Symbol, merged_cfg::Dict{String, Any})
    merged_keys = flatten_keys(merged_cfg)
    for key in REQUIRED_KEYS[domain]
        if !(key in merged_keys)
            push!(violations, string(domain, ": merged config missing required key ", key))
        end
    end
end

function main()
    violations = String[]
    infos = String[]

    println("[config-profile-matrix] checking njl/pnjl/rpnjl")

    for domain in DOMAINS
        default_path = profile_path(domain, "default")
        unittest_path = profile_path(domain, "unittest")

        isfile(default_path) || begin
            push!(violations, string(domain, ": missing default profile at ", default_path))
            continue
        end

        default_cfg = parse_profile(default_path)
        unittest_cfg = parse_profile(unittest_path)
        merged_cfg = deep_merge(default_cfg, unittest_cfg)

        stats = summarize_overrides(default_cfg, unittest_cfg)

        assert_required_keys!(violations, domain, merged_cfg)

        if !isempty(stats.unknown_override)
            push!(violations,
                string(domain, ": unittest contains unknown override keys: ", join(stats.unknown_override, ", "))
            )
        end

        push!(infos,
            string(
                "  - ",
                domain,
                ": default=",
                default_path,
                ", unittest=",
                isfile(unittest_path) ? unittest_path : "<missing: inherit default>",
                ", overrides=",
                length(stats.override_keys)
            )
        )

        if !isempty(stats.override_keys)
            push!(infos, string("    overrides keys: ", join(stats.override_keys, ", ")))
        end
    end

    for line in infos
        println(line)
    end

    if !isempty(violations)
        println("[config-profile-matrix] FAILED")
        for item in violations
            println(" - " * item)
        end
        exit(1)
    end

    println("[config-profile-matrix] OK")
    println("  rule: unittest profile is treated as override-only layer on top of default")
end

main()
