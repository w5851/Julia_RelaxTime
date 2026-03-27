#!/usr/bin/env julia

using Statistics

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "src", "config", "ConfigLoader.jl"))
if !isdefined(Main, :ConfigLoader)
    Base.include(Main, _CONFIG_LOADER_PATH)
end

using Main.ConfigLoader: load_config, reset_config_loader_cache!

const SAMPLE_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "config", "physics"))
const SAMPLE_PROFILE = "default"
const WARMUP = 20
const SAMPLES = 1000

function _measure_once(config_dir, default_config; profile)
    timed = @timed load_config(config_dir, default_config; profile=profile)
    allocs = timed.gcstats.poolalloc + timed.gcstats.bigalloc + timed.gcstats.malloc + timed.gcstats.realloc
    return (time_ns=Int(round(timed.time * 1.0e9)), allocs=allocs, bytes=timed.bytes)
end

function _run_scenario(scenario::String, config_dir, default_config; profile, reset_each_call::Bool)
    times = Vector{Int}(undef, SAMPLES)
    allocs = Vector{Int}(undef, SAMPLES)
    bytes = Vector{Int}(undef, SAMPLES)

    for _ in 1:WARMUP
        reset_each_call && reset_config_loader_cache!()
        _ = _measure_once(config_dir, default_config; profile=profile)
    end

    GC.gc()
    for i in 1:SAMPLES
        reset_each_call && reset_config_loader_cache!()
        measured = _measure_once(config_dir, default_config; profile=profile)
        times[i] = measured.time_ns
        allocs[i] = measured.allocs
        bytes[i] = measured.bytes
    end

    return (
        scenario=scenario,
        median_time_ns=Int(round(median(times))),
        allocs=Int(round(median(allocs))),
        bytes=Int(round(median(bytes))),
    )
end

function main()
    default_config = Dict{String, Any}()
    println("env.julia_version=$(VERSION)")
    println("env.threads=$(Threads.nthreads())")
    println("env.profile=$(SAMPLE_PROFILE)")
    println("env.sample_path=$(joinpath(SAMPLE_CONFIG_DIR, string(SAMPLE_PROFILE, ".toml")))")
    println("env.warmup=$(WARMUP)")
    println("env.samples=$(SAMPLES)")
    println("scenario,median_time_ns,allocs,bytes")

    reset_config_loader_cache!()
    before = _run_scenario("before", SAMPLE_CONFIG_DIR, default_config; profile=SAMPLE_PROFILE, reset_each_call=true)

    reset_config_loader_cache!()
    _ = load_config(SAMPLE_CONFIG_DIR, default_config; profile=SAMPLE_PROFILE)
    after = _run_scenario("after", SAMPLE_CONFIG_DIR, default_config; profile=SAMPLE_PROFILE, reset_each_call=false)

    println("$(before.scenario),$(before.median_time_ns),$(before.allocs),$(before.bytes)")
    println("$(after.scenario),$(after.median_time_ns),$(after.allocs),$(after.bytes)")
end

main()
