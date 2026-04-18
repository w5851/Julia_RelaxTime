"""
Unified relaxtime scan entrypoint.

Subcommands:
- gap-transport   -> scripts/relaxtime/run_gap_transport_scan.jl
- tau-vs-t        -> scripts/relaxtime/scan_relaxation_times_vs_T.jl
- manual-workflow -> scripts/relaxtime/run_manual_relaxation_scan_workflow.jl
"""

module SingleEntryScanCLI

export parse_entry_args, resolve_target_script, run_subcommand, main

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const _SUBCOMMAND_MAP = Dict{String,Symbol}(
    "gap-transport" => :gap_transport,
    "tau-vs-t" => :tau_vs_t,
    "manual-workflow" => :manual_workflow,
)

function print_usage(io::IO=stdout)
    println(io, "Usage: julia --project=. scripts/relaxtime/run_scan.jl <subcommand> [options]")
    println(io)
    println(io, "Subcommands:")
    println(io, "  gap-transport    Scan (T, mu_B, xi) and compute gap + transport")
    println(io, "  tau-vs-t         Scan relaxation times versus temperature")
    println(io, "  manual-workflow  Run manual relaxation scan workflow")
    println(io)
    println(io, "Use '<subcommand> --help' to view subcommand options.")
end

function parse_entry_args(args::Vector{String})::Tuple{Symbol,Vector{String}}
    isempty(args) && throw(ArgumentError("missing subcommand"))
    first = args[1]
    if first in ("-h", "--help")
        throw(ArgumentError("help requested"))
    end
    sub = get(_SUBCOMMAND_MAP, first, nothing)
    sub === nothing && throw(ArgumentError("unknown subcommand: $(first)"))
    return sub, args[2:end]
end

function resolve_target_script(subcommand::Symbol)::String
    if subcommand === :gap_transport
        return joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_gap_transport_scan.jl")
    elseif subcommand === :tau_vs_t
        return joinpath(PROJECT_ROOT, "scripts", "relaxtime", "scan_relaxation_times_vs_T.jl")
    elseif subcommand === :manual_workflow
        return joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_manual_relaxation_scan_workflow.jl")
    end
    throw(ArgumentError("unsupported subcommand: $(subcommand)"))
end

function _julia_cmd(argv::Vector{String})::Cmd
    jc = Base.julia_cmd()
    cmd = Cmd(vcat(jc.exec, argv))
    return jc.env === nothing ? cmd : setenv(cmd, jc.env)
end

function run_subcommand(subcommand::Symbol, pass_args::Vector{String})::Int
    target = resolve_target_script(subcommand)
    isfile(target) || throw(ArgumentError("target script not found: $(target)"))

    rel_target = relpath(target, PROJECT_ROOT)
    argv = vcat(String["--project=.", rel_target], pass_args)

    proc = run(ignorestatus(_julia_cmd(argv)))
    return proc.exitcode
end

function main(args::Vector{String}=copy(ARGS))::Int
    if isempty(args) || args[1] in ("-h", "--help")
        print_usage(stdout)
        return 0
    end

    sub, rest = parse_entry_args(args)
    return run_subcommand(sub, rest)
end

end # module SingleEntryScanCLI

if abspath(PROGRAM_FILE) == @__FILE__
    exit(Main.SingleEntryScanCLI.main(copy(ARGS)))
end
