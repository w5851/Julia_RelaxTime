using BenchmarkTools

const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end

const PC = Main.PolarizationCache

@inline function _round_cache_sigdigits(x::Float64)
    (x == 0.0 || isnan(x) || isinf(x)) && return x
    return round(x; sigdigits=12)
end

struct SigdigitsPolarizationKey
    channel_code::UInt8
    k0::Float64
    k_norm::Float64
    m1::Float64
    m2::Float64
    μ1::Float64
    μ2::Float64
    T::Float64
    Φ::Float64
    Φbar::Float64
    ξ::Float64
    A1::Float64
    A2::Float64
    num_s_quark::UInt8

    function SigdigitsPolarizationKey(channel::Symbol, k0::Float64, k_norm::Float64,
                                      m1::Float64, m2::Float64, μ1::Float64, μ2::Float64,
                                      T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64,
                                      A1::Float64, A2::Float64, num_s_quark::Int)
        return new(
            PC._channel_code(channel),
            _round_cache_sigdigits(k0),
            _round_cache_sigdigits(k_norm),
            _round_cache_sigdigits(m1),
            _round_cache_sigdigits(m2),
            _round_cache_sigdigits(μ1),
            _round_cache_sigdigits(μ2),
            _round_cache_sigdigits(T),
            _round_cache_sigdigits(Φ),
            _round_cache_sigdigits(Φbar),
            _round_cache_sigdigits(ξ),
            _round_cache_sigdigits(A1),
            _round_cache_sigdigits(A2),
            UInt8(num_s_quark),
        )
    end
end

@inline function Base.isequal(a::SigdigitsPolarizationKey, b::SigdigitsPolarizationKey)
    return a.channel_code == b.channel_code &&
           isequal(a.k0, b.k0) &&
           isequal(a.k_norm, b.k_norm) &&
           isequal(a.m1, b.m1) &&
           isequal(a.m2, b.m2) &&
           isequal(a.μ1, b.μ1) &&
           isequal(a.μ2, b.μ2) &&
           isequal(a.T, b.T) &&
           isequal(a.Φ, b.Φ) &&
           isequal(a.Φbar, b.Φbar) &&
           isequal(a.ξ, b.ξ) &&
           isequal(a.A1, b.A1) &&
           isequal(a.A2, b.A2) &&
           a.num_s_quark == b.num_s_quark
end

@inline function Base.hash(k::SigdigitsPolarizationKey, h::UInt)
    h = hash(k.channel_code, h)
    h = hash(k.k0, h)
    h = hash(k.k_norm, h)
    h = hash(k.m1, h)
    h = hash(k.m2, h)
    h = hash(k.μ1, h)
    h = hash(k.μ2, h)
    h = hash(k.T, h)
    h = hash(k.Φ, h)
    h = hash(k.Φbar, h)
    h = hash(k.ξ, h)
    h = hash(k.A1, h)
    h = hash(k.A2, h)
    h = hash(k.num_s_quark, h)
    return h
end

const _BASE_SAMPLE = (
    channel=:P,
    k0=100.0,
    k_norm=50.0,
    m1=5.0,
    m2=5.0,
    μ1=300.0,
    μ2=300.0,
    T=150.0,
    Φ=0.5,
    Φbar=0.5,
    ξ=0.2,
    A1=-50.0,
    A2=-50.0,
    num_s=0,
)

function _sample_grid()
    scales = (0.0, 1e-14, 5e-13, 1e-12, 5e-12, 1e-10)
    channels = (:P, :S)
    samples = NamedTuple[]
    for channel in channels
        for δk0 in scales, δk in scales, δμ in (0.0, 1e-13, 1e-11)
            push!(samples, (
                channel=channel,
                k0=_BASE_SAMPLE.k0 * (1 + δk0),
                k_norm=_BASE_SAMPLE.k_norm * (1 + δk),
                m1=_BASE_SAMPLE.m1,
                m2=_BASE_SAMPLE.m2,
                μ1=_BASE_SAMPLE.μ1 * (1 + δμ),
                μ2=_BASE_SAMPLE.μ2,
                T=_BASE_SAMPLE.T,
                Φ=_BASE_SAMPLE.Φ,
                Φbar=_BASE_SAMPLE.Φbar,
                ξ=_BASE_SAMPLE.ξ,
                A1=_BASE_SAMPLE.A1,
                A2=_BASE_SAMPLE.A2,
                num_s=_BASE_SAMPLE.num_s,
            ))
        end
    end
    return samples
end

const SAMPLE_GRID = _sample_grid()
const SCALAR_VALUES = vcat(
    [s.k0 for s in SAMPLE_GRID],
    [s.k_norm for s in SAMPLE_GRID],
    [s.μ1 for s in SAMPLE_GRID],
    [s.T for s in SAMPLE_GRID],
)

@inline function _current_key(sample)
    return PC.PolarizationKey(
        sample.channel,
        sample.k0,
        sample.k_norm,
        sample.m1,
        sample.m2,
        sample.μ1,
        sample.μ2,
        sample.T,
        sample.Φ,
        sample.Φbar,
        sample.ξ,
        sample.A1,
        sample.A2,
        sample.num_s,
    )
end

@inline function _sigdigits_key(sample)
    return SigdigitsPolarizationKey(
        sample.channel,
        sample.k0,
        sample.k_norm,
        sample.m1,
        sample.m2,
        sample.μ1,
        sample.μ2,
        sample.T,
        sample.Φ,
        sample.Φbar,
        sample.ξ,
        sample.A1,
        sample.A2,
        sample.num_s,
    )
end

function _make_current_dict(samples)
    d = Dict{PC.PolarizationKey, Tuple{Float64, Float64}}()
    @inbounds for (i, sample) in pairs(samples)
        d[_current_key(sample)] = (Float64(i), -Float64(i))
    end
    return d
end

function _make_sigdigits_dict(samples)
    d = Dict{SigdigitsPolarizationKey, Tuple{Float64, Float64}}()
    @inbounds for (i, sample) in pairs(samples)
        d[_sigdigits_key(sample)] = (Float64(i), -Float64(i))
    end
    return d
end

const CURRENT_DICT = _make_current_dict(SAMPLE_GRID)
const SIGDIGITS_DICT = _make_sigdigits_dict(SAMPLE_GRID)

function _lookup_current(samples)
    acc = 0.0
    @inbounds for sample in samples
        acc += get(CURRENT_DICT, _current_key(sample), (0.0, 0.0))[1]
    end
    return acc
end

function _lookup_sigdigits(samples)
    acc = 0.0
    @inbounds for sample in samples
        acc += get(SIGDIGITS_DICT, _sigdigits_key(sample), (0.0, 0.0))[1]
    end
    return acc
end

function _bucket_stats(samples)
    current_keys = [_current_key(s) for s in samples]
    sigdigits_keys = [_sigdigits_key(s) for s in samples]
    current_unique = length(Set(current_keys))
    sigdigits_unique = length(Set(sigdigits_keys))
    return (
        sample_count=length(samples),
        current_unique=current_unique,
        sigdigits_unique=sigdigits_unique,
        unique_delta=sigdigits_unique - current_unique,
    )
end

function main()
    stats = _bucket_stats(SAMPLE_GRID)
    println("PolarizationCache quantization evaluation")
    println("  sample_count      = $(stats.sample_count)")
    println("  current_unique    = $(stats.current_unique)")
    println("  sigdigits_unique  = $(stats.sigdigits_unique)")
    println("  unique_delta      = $(stats.unique_delta)")

    current_scalar = @benchmark PC._round_cache.(SCALAR_VALUES)
    sigdigits_scalar = @benchmark _round_cache_sigdigits.(SCALAR_VALUES)
    current_key = @benchmark [_current_key(s) for s in SAMPLE_GRID]
    sigdigits_key = @benchmark [_sigdigits_key(s) for s in SAMPLE_GRID]
    current_lookup = @benchmark _lookup_current(SAMPLE_GRID)
    sigdigits_lookup = @benchmark _lookup_sigdigits(SAMPLE_GRID)

    println("\n[scalar quantization]")
    show(stdout, MIME("text/plain"), current_scalar)
    println()
    show(stdout, MIME("text/plain"), sigdigits_scalar)
    println()

    println("\n[key construction]")
    show(stdout, MIME("text/plain"), current_key)
    println()
    show(stdout, MIME("text/plain"), sigdigits_key)
    println()

    println("\n[cache-hit lookup]")
    show(stdout, MIME("text/plain"), current_lookup)
    println()
    show(stdout, MIME("text/plain"), sigdigits_lookup)
    println()
end

main()