#!/usr/bin/env julia
# Fixed-T aliases for the shared thermal-weighted scan; common CSV/provenance.
if !isdefined(Main,:ScanOptions)
    Base.include(Main,joinpath(@__DIR__,"run_gap_transport_scan.jl"))
end

module FixedTemperatureTransportScan

function parse_args(args::Vector{String})
    T_MeV=150.0;muB_MeV=800.0
    xi_min=0.0;xi_max=1.0;xi_step=0.1
    explicit_xi=nothing;xi0_only=false
    forwarded=String["--output",joinpath("data","outputs","results","relaxtime","scan","transport_vs_xi_T150_muB800.csv"),
        "--no-compute-bulk"]
    i=1
    while i<=length(args)
        arg=args[i]
        function value()
            i<length(args) || throw(ArgumentError("missing value for $arg"))
            i+=1
            return args[i]
        end
        if arg=="--T-MeV"
            T_MeV=parse(Float64,value())
        elseif arg=="--muB-MeV"
            muB_MeV=parse(Float64,value())
        elseif arg=="--xi-min"
            xi_min=parse(Float64,value())
        elseif arg=="--xi-max"
            xi_max=parse(Float64,value())
        elseif arg=="--xi-step"
            xi_step=parse(Float64,value())
        elseif arg=="--xi-list"
            explicit_xi=parse.(Float64,split(value(),','))
        elseif arg=="--xi"
            explicit_xi===nothing && (explicit_xi=Float64[])
            push!(explicit_xi,parse(Float64,value()))
        elseif arg=="--xi0-only"
            xi0_only=true
        elseif arg in ("--sigma-cache-dir","--ignore-lock","--dump-debug","--debug-out")
            throw(ArgumentError("$arg was retired; use the sampling profile and --channel-diagnostics-output"))
        elseif arg in ("-h","--help")
            println("Fixed-T scan: --T-MeV 150 --muB-MeV 800; --xi-min 0 --xi-max 1 --xi-step 0.1.")
            println("--xi-list / --xi / --xi0-only select explicit xi; --out aliases --output.")
            Main.GapTransportScanCLI.print_usage()
            return nothing
        else
            push!(forwarded,arg=="--out" ? "--output" : arg)
        end
        i+=1
    end
    all(isfinite,(T_MeV,muB_MeV,xi_min,xi_max,xi_step)) && T_MeV>0 && xi_step>0 && xi_min<=xi_max ||
        throw(ArgumentError("finite T>0 and a valid xi interval/step are required"))
    xis=xi0_only ? [0.0] : explicit_xi===nothing ? collect(xi_min:xi_step:xi_max) : explicit_xi
    !isempty(xis) && all(x->isfinite(x) && x> -1,xis) || throw(ArgumentError("xi values must be finite and > -1"))
    append!(forwarded,["--tmin",string(T_MeV),"--tmax",string(T_MeV),
        "--mubmin",string(muB_MeV),"--mubmax",string(muB_MeV),"--xi-list",join(xis,',')])
    return Main.GapTransportScanCLI.parse_args(forwarded)
end

function main(args=copy(ARGS))
    opts=parse_args(args)
    opts===nothing && return
    reference=Main._load_runtime_phase_reference(opts)
    ctx=Main.ProvenanceMetadata.new_run_context("scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl",args)
    return Main.run_scan(opts,ctx;phase_reference=reference)
end

end # module FixedTemperatureTransportScan

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    FixedTemperatureTransportScan.main()
end
