"""Explicit infinite-thermal charged GBU research-production entrypoint."""
module ChargedGBUFreezeoutCLI
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
function parse_args(args)
    opts=Dict{Symbol,Any}(:output=>joinpath(ROOT,"data","outputs","results","relaxtime","meson_density",
        "charged_gbu_infinite","freezeout_20260907"),:resume=>false,:make_plot=>true)
    i=1
    while i<=length(args)
        a=args[i]
        if a=="--resume";opts[:resume]=true
        elseif a=="--no-plot";opts[:make_plot]=false
        elseif a in ("--help","-h");return nothing
        elseif a=="--workers"
            i<length(args) || error("missing worker count");i+=1;opts[:workers]=parse(Int,args[i])
        elseif a in ("--output","--config","--sqrts-list")
            i<length(args) || error("missing value for $(a)");i+=1
            key=a=="--output" ? :output : a=="--config" ? :config : :energies
            opts[key]=a=="--sqrts-list" ? parse.(Float64,split(args[i],',')) : args[i]
        else;error("unknown option: $(a)")
        end
        i+=1
    end
    return opts
end
function main(args=ARGS)
    opts=parse_args(args)
    if opts===nothing
        println("Usage: julia --project=. scripts/relaxtime/run_charged_gbu_freezeout_scan.jl [--output DIR] [--config TOML] [--sqrts-list 3,5,...] [--workers 1..4] [--resume] [--no-plot]")
        return
    end
    isdefined(Main,:Models) || Base.include(Main,joinpath(ROOT,"src","models","Models.jl"))
    result=Base.invokelatest(() -> Main.Models.run_charged_gbu_freezeout_scan(;opts...))
    println(result)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
