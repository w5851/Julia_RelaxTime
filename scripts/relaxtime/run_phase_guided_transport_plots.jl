#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "plot_scan_csv.py")

using JSON3

struct Options
    case_dir::String
    csv_path::String
    fig_dir::String
    python_cmd::String
    overwrite::Bool
end

function print_usage(io::IO=stdout)
    println(io, "Usage: julia --project=. scripts/relaxtime/run_phase_guided_transport_plots.jl --case-dir <dir> [options]\n")
    println(io, "Options:")
    println(io, "  --case-dir <dir>   canonical case directory containing phase_guided_transport_scan.csv")
    println(io, "  --csv <path>       optional explicit CSV path (default <case-dir>/phase_guided_transport_scan.csv)")
    println(io, "  --fig-dir <dir>    optional figure output dir (default <case-dir>/figures)")
    println(io, "  --python <cmd>     python executable (default python)")
    println(io, "  --overwrite        clear existing figures directory before plotting")
    println(io, "  -h, --help         show help")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :case_dir => nothing,
        :csv_path => nothing,
        :fig_dir => nothing,
        :python_cmd => "python",
        :overwrite => false,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $(arg)")
            i += 1
            return args[i]
        end

        if arg == "--case-dir"
            opts[:case_dir] = require_value()
        elseif arg == "--csv"
            opts[:csv_path] = require_value()
        elseif arg == "--fig-dir"
            opts[:fig_dir] = require_value()
        elseif arg == "--python"
            opts[:python_cmd] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $(arg)")
        end
        i += 1
    end

    isnothing(opts[:case_dir]) && error("--case-dir is required")
    case_dir = normpath(abspath(String(opts[:case_dir])))
    csv_path = isnothing(opts[:csv_path]) ? joinpath(case_dir, "phase_guided_transport_scan.csv") : normpath(abspath(String(opts[:csv_path])))
    fig_dir = isnothing(opts[:fig_dir]) ? joinpath(case_dir, "figures") : normpath(abspath(String(opts[:fig_dir])))

    isdir(case_dir) || error("case directory not found: $(case_dir)")
    isfile(csv_path) || error("case csv not found: $(csv_path)")
    isfile(PLOT_SCRIPT) || error("plot script not found: $(PLOT_SCRIPT)")

    return Options(case_dir, csv_path, fig_dir, String(opts[:python_cmd]), Bool(opts[:overwrite]))
end

function _python_bin(cmd::String)
    if Sys.which(cmd) !== nothing
        return String(Sys.which(cmd))
    end
    cmd == "python" && Sys.which("python3") !== nothing && return String(Sys.which("python3"))
    error("python interpreter not found: $(cmd)")
end

function _run_plot(py::String, csv_path::String, fig_dir::String, ys::String; multi_y::Bool=false)
    args = String[
        PLOT_SCRIPT,
        "--mode", "lines",
        "--csv", csv_path,
        "--x", "xi",
        "--ys", ys,
        "--split", "scan_group",
        "--out-dir", fig_dir,
        "--check",
    ]
    multi_y && append!(args, ["--multi-y"])
    cmd = Cmd(vcat([py], args))
    run(Cmd(cmd; dir=PROJECT_ROOT))
end

function _write_plot_manifest(fig_dir::String)
    files = String[]
    for path in sort(readdir(fig_dir; join=true))
        isdir(path) || continue
        for sub in sort(readdir(path; join=true))
            isfile(sub) || continue
            endswith(lowercase(sub), ".png") || continue
            push!(files, replace(relpath(sub, fig_dir), '\\' => '/'))
        end
    end
    payload = Dict(
        "schema_version" => "v1",
        "figures" => files,
        "count" => length(files),
    )
    open(joinpath(fig_dir, "plot_manifest.json"), "w") do io
        write(io, JSON3.write(payload))
    end
end

function _append_readme(case_dir::String, fig_dir::String)
    readme_path = joinpath(case_dir, "README.md")
    isfile(readme_path) || return

    files = String[]
    for path in sort(readdir(fig_dir; join=true))
        isdir(path) || continue
        for sub in sort(readdir(path; join=true))
            isfile(sub) || continue
            endswith(lowercase(sub), ".png") || continue
            push!(files, replace(relpath(sub, case_dir), '\\' => '/'))
        end
    end

    text = read(readme_path, String)
    marker = "## Generated Figures"
    if occursin(marker, text)
        text = first(split(text, marker; limit=2))
        text = rstrip(text) * "\n\n"
    else
        text = rstrip(text) * "\n\n"
    end

    open(readme_path, "w") do io
        write(io, text)
        println(io, marker)
        println(io, "- tau flavor-resolved same-panel figures by `scan_group`")
        println(io, "- individual `eta`, `sigma`, `zeta`, `eta_over_s`, `sigma_over_T` figures by `scan_group`")
        println(io, "- plot manifest: `figures/plot_manifest.json`")
        if !isempty(files)
            println(io)
            println(io, "Generated PNG files:")
            for file in files
                println(io, "- `", file, "`")
            end
        end
    end
end

function run_phase_guided_plots(opts::Options)
    py = _python_bin(opts.python_cmd)
    if opts.overwrite && isdir(opts.fig_dir)
        rm(opts.fig_dir; recursive=true, force=true)
    end
    mkpath(opts.fig_dir)

    _run_plot(py, opts.csv_path, opts.fig_dir, "tau_u,tau_d,tau_s"; multi_y=true)
    for y in ("eta", "sigma", "zeta", "eta_over_s", "sigma_over_T")
        _run_plot(py, opts.csv_path, opts.fig_dir, y; multi_y=false)
    end

    _write_plot_manifest(opts.fig_dir)
    _append_readme(opts.case_dir, opts.fig_dir)
    println("Phase-guided transport plotting finished. Output: $(opts.fig_dir)")
end

function main(args::Vector{String}=copy(ARGS))
    opts = parse_args(args)
    run_phase_guided_plots(opts)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
