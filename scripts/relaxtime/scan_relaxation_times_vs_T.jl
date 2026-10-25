#!/usr/bin/env julia
# Temperature-scan aliases for the shared qualified transport entrypoint.
# Results now use the common gap_transport_scan CSV schema; all six tau columns,
# channel diagnostics and the sampling identity are retained.
if !isdefined(Main,:ScanOptions)
    Base.include(Main,joinpath(@__DIR__,"run_gap_transport_scan.jl"))
end

module RelaxationTemperatureScan

function parse_args(args::Vector{String})
    forwarded=String["--output",joinpath("data","outputs","results","relaxtime","scan","relaxation_times_vs_T.csv"),
        "--tmin","80","--tmax","220","--tstep","10","--no-compute-bulk"]
    muB_values=[0.0,800.0]
    i=1
    while i<=length(args)
        arg=args[i]
        if arg=="--mub-list"
            i<length(args) || throw(ArgumentError("missing value for --mub-list"))
            i+=1
            muB_values=parse.(Float64,split(args[i],','))
        elseif arg in ("-h","--help")
            println("Temperature scan: --out aliases --output; --diagnostics-out aliases --channel-diagnostics-output.")
            println("--mub-list 0,800 selects an explicit muB list. Remaining options use run_gap_transport_scan.jl.")
            Main.GapTransportScanCLI.print_usage()
            return nothing
        else
            push!(forwarded,arg=="--out" ? "--output" : arg=="--diagnostics-out" ? "--channel-diagnostics-output" : arg)
        end
        i+=1
    end
    !isempty(muB_values) && all(isfinite,muB_values) || throw(ArgumentError("invalid --mub-list"))
    muB_values=sort!(unique!(muB_values))
    append!(forwarded,["--mubmin",string(first(muB_values)),"--mubmax",string(last(muB_values))])
    opts=Main.GapTransportScanCLI.parse_args(forwarded)
    return (scan_options=opts,muB_values=muB_values)
end

function main(args=copy(ARGS))
    parsed=parse_args(args)
    parsed===nothing && return
    opts=parsed.scan_options
    reference=Main._load_runtime_phase_reference(opts)
    ctx=Main.ProvenanceMetadata.new_run_context("scripts/relaxtime/scan_relaxation_times_vs_T.jl",args)
    return Main.run_scan(opts,ctx;phase_reference=reference,muB_values=parsed.muB_values)
end

end # module RelaxationTemperatureScan

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    RelaxationTemperatureScan.main()
end
