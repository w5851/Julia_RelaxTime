#!/usr/bin/env julia

using Dates

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SAMPLING_LIB_PATH = joinpath(@__DIR__, "xi_smoothness_sampling_lib.jl")

Base.include(Main, SAMPLING_LIB_PATH)
using Main.XiSmoothnessSampling: SampleRow, sample_params

@inline function _fmt_anchor_float(x::Float64)
    return isfinite(x) ? string(x) : ""
end

function _print_usage()
    println("Usage: julia --project=. scripts/relaxtime/generate_xi_smoothness_params.jl [options]")
    println("Options:")
    println("  --seed <int>           RNG seed (default: 42)")
    println("  --output <path>        output CSV path")
    println("  --random-count <int>   random sample count (default: 12)")
    println("  --near-count <int>     near-phase-line sample count (default: 12)")
    println("  --tmin <float>         T lower bound in MeV (default: 50)")
    println("  --tmax <float>         T upper bound in MeV (default: 270)")
    println("  --muqmin <float>       mu_q lower bound in MeV (default: 0)")
    println("  --muqmax <float>       mu_q upper bound in MeV (default: 360)")
    println("  --boundary-csv <path>  boundary csv path (default: v2 accepted table)")
    println("  --crossover-csv <path> crossover csv path (default: v2 accepted table)")
end

function _parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :seed => 42,
        :output => nothing,
        :random_count => 12,
        :near_count => 12,
        :tmin => 50.0,
        :tmax => 270.0,
        :muqmin => 0.0,
        :muqmax => 360.0,
        :boundary_csv => joinpath(
            PROJECT_ROOT, "data", "reference", "pnjl", "issue130_phase_reference_v2",
            "accepted", "tables", "maxwell_surface_accepted_phase_map_v1.csv",
        ),
        :crossover_csv => joinpath(
            PROJECT_ROOT, "data", "reference", "pnjl", "issue130_phase_reference_v2",
            "accepted", "tables", "crossover_surface_accepted_phase_map_v1.csv",
        ),
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            _print_usage()
            exit(0)
        end

        if !startswith(arg, "--")
            throw(ArgumentError("unknown argument: $(arg)"))
        end

        key = arg[3:end]
        i == length(args) && throw(ArgumentError("missing value for $(arg)"))
        i += 1
        val = args[i]

        if key == "seed"
            cfg[:seed] = parse(Int, val)
        elseif key == "output"
            cfg[:output] = val
        elseif key == "random-count"
            cfg[:random_count] = parse(Int, val)
        elseif key == "near-count"
            cfg[:near_count] = parse(Int, val)
        elseif key == "tmin"
            cfg[:tmin] = parse(Float64, val)
        elseif key == "tmax"
            cfg[:tmax] = parse(Float64, val)
        elseif key == "muqmin"
            cfg[:muqmin] = parse(Float64, val)
        elseif key == "muqmax"
            cfg[:muqmax] = parse(Float64, val)
        elseif key == "boundary-csv"
            cfg[:boundary_csv] = val
        elseif key == "crossover-csv"
            cfg[:crossover_csv] = val
        else
            throw(ArgumentError("unknown option: --$(key)"))
        end
        i += 1
    end

    cfg[:random_count] >= 0 || throw(ArgumentError("--random-count must be non-negative"))
    cfg[:near_count] >= 0 || throw(ArgumentError("--near-count must be non-negative"))
    total = cfg[:random_count] + cfg[:near_count]
    total > 0 || throw(ArgumentError("total sample count must be positive"))

    if cfg[:output] === nothing
        total = cfg[:random_count] + cfg[:near_count]
        cfg[:output] = joinpath(
            PROJECT_ROOT,
            "data",
            "outputs",
            "results",
            "relaxtime",
            "xi_smoothness_sampling",
            "sampling",
            "params_$(total)_seed$(cfg[:seed]).csv",
        )
    end

    return cfg
end

function _write_output(path::String, rows::Vector{SampleRow}; cfg::Dict{Symbol, Any})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# schema: xi_smoothness_sampling_v1")
        println(io, "# generated_utc: " * Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"))
        println(io, "# seed: " * string(cfg[:seed]))
        println(io, "# random_count: " * string(cfg[:random_count]))
        println(io, "# near_count: " * string(cfg[:near_count]))
        println(io, "# T_range_MeV: [" * string(cfg[:tmin]) * ", " * string(cfg[:tmax]) * "]")
        println(io, "# muq_range_MeV: [" * string(cfg[:muqmin]) * ", " * string(cfg[:muqmax]) * "]")
        println(io, "sample_id,source,T_MeV,muq_MeV,muB_MeV,anchor_type,anchor_T_MeV,anchor_muq_MeV,delta_T,delta_muq,rng_seed")
        for r in rows
            anchor_type = r.source == "random_uniform" ? "" : r.anchor_type
            println(io, join((
                r.sample_id,
                r.source,
                string(r.T_MeV),
                string(r.muq_MeV),
                string(r.muB_MeV),
                anchor_type,
                _fmt_anchor_float(r.anchor_T_MeV),
                _fmt_anchor_float(r.anchor_muq_MeV),
                _fmt_anchor_float(r.delta_T),
                _fmt_anchor_float(r.delta_muq),
                string(r.rng_seed),
            ), ","))
        end
    end
end

function main()
    cfg = _parse_args(ARGS)
    total = cfg[:random_count] + cfg[:near_count]

    rows = sample_params(
        total;
        seed=cfg[:seed],
        random_count=cfg[:random_count],
        near_count=cfg[:near_count],
        T_range=(cfg[:tmin], cfg[:tmax]),
        muq_range=(cfg[:muqmin], cfg[:muqmax]),
        boundary_csv=String(cfg[:boundary_csv]),
        crossover_csv=String(cfg[:crossover_csv]),
    )

    output = String(cfg[:output])
    _write_output(output, rows; cfg=cfg)
    println("Wrote sampling params CSV: " * output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
