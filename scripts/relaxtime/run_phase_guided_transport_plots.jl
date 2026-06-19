#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "plot_scan_csv.py")

using JSON3
using SHA: sha256

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
    println(io, "  --fig-dir <dir>    optional figure output dir (default data/outputs/figures/relaxtime/transport/phase_guided/<mode>/<case_name>)")
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
    fig_dir = isnothing(opts[:fig_dir]) ? _default_fig_dir(case_dir) : normpath(abspath(String(opts[:fig_dir])))

    isdir(case_dir) || error("case directory not found: $(case_dir)")
    isfile(csv_path) || error("case csv not found: $(csv_path)")
    isfile(PLOT_SCRIPT) || error("plot script not found: $(PLOT_SCRIPT)")

    return Options(case_dir, csv_path, fig_dir, String(opts[:python_cmd]), Bool(opts[:overwrite]))
end

function _default_fig_dir(case_dir::String)
    mode_dir = basename(dirname(case_dir))
    case_name = basename(case_dir)
    config_path = joinpath(case_dir, "effective_config.json")
    if isfile(config_path)
        cfg = JSON3.read(read(config_path, String))
        config_view = haskey(cfg, "options") ? cfg["options"] : cfg
        if haskey(config_view, "mode")
            mode_value = String(config_view["mode"])
            if mode_value == "mode_a_fixed_muB_phase_scaled"
                mode_dir = "mode_a_fixed_muB_phase_scaled"
            elseif mode_value == "mode_b_fixed_T_sparse_muB"
                mode_dir = "mode_b_fixed_T_sparse_muB"
            end
        end
        haskey(config_view, "case_name") && (case_name = String(config_view["case_name"]))
    end
    return normpath(joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "transport", "phase_guided", mode_dir, case_name))
end

function _python_bin(cmd::String)
    if Sys.which(cmd) !== nothing
        return String(Sys.which(cmd))
    end
    cmd == "python" && Sys.which("python3") !== nothing && return String(Sys.which("python3"))
    error("python interpreter not found: $(cmd)")
end

function _run_plot(py::String, csv_path::String, fig_dir::String, ys::String; split::String, group::String)
    args = String[
        PLOT_SCRIPT,
        "--mode", "lines",
        "--csv", csv_path,
        "--x", "xi",
        "--ys", ys,
        "--split", split,
        "--group", group,
        "--out-dir", fig_dir,
        "--check",
    ]
    cmd = Cmd(vcat([py], args))
    run(Cmd(cmd; dir=PROJECT_ROOT))
end

function _detect_mode(case_dir::String, csv_path::String)
    config_path = joinpath(case_dir, "effective_config.json")
    if isfile(config_path)
        cfg = JSON3.read(read(config_path, String))
        config_view = haskey(cfg, "options") ? cfg["options"] : cfg
        if haskey(config_view, "mode")
            mode_value = String(config_view["mode"])
            mode_value == "mode_a_fixed_muB_phase_scaled" && return :mode_a_fixed_muB_phase_scaled
            mode_value == "mode_b_fixed_T_sparse_muB" && return :mode_b_fixed_T_sparse_muB
        end
    end
    open(csv_path, "r") do io
        header_seen = false
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            if !header_seen
                header_seen = true
                continue
            end
            parts = split(s, ',')
            length(parts) >= 5 || continue
            mode_value = strip(parts[5])
            mode_value == "mode_a_fixed_muB_phase_scaled" && return :mode_a_fixed_muB_phase_scaled
            mode_value == "mode_b_fixed_T_sparse_muB" && return :mode_b_fixed_T_sparse_muB
        end
    end
    error("unable to detect phase-guided mode from csv: $(csv_path)")
end

function _plot_layout(mode::Symbol)
    if mode == :mode_a_fixed_muB_phase_scaled
        return (; split="plot_panel", group="plot_series_label", panel_desc="fixed mu_B panel", series_desc="alpha_T lines with explicit T")
    elseif mode == :mode_b_fixed_T_sparse_muB
        return (; split="plot_panel", group="plot_series_label", panel_desc="fixed T panel", series_desc="mu_B lines")
    end
    error("unsupported phase-guided mode: $(mode)")
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

function _sha256_file(path::String)
    open(path, "r") do io
        return bytes2hex(sha256(io))
    end
end

function _plot_output_summary(fig_dir::String)
    png_files = String[]
    sidecar_files = String[]
    min_png_bytes = nothing
    for path in sort(readdir(fig_dir; join=true))
        isdir(path) || continue
        for sub in sort(readdir(path; join=true))
            isfile(sub) || continue
            rel = replace(relpath(sub, fig_dir), '\\' => '/')
            if endswith(lowercase(sub), ".png")
                push!(png_files, rel)
                size = filesize(sub)
                min_png_bytes = min_png_bytes === nothing ? size : min(min_png_bytes, size)
            elseif endswith(lowercase(sub), ".png.provenance.json")
                push!(sidecar_files, rel)
            end
        end
    end
    return Dict(
        "png_count" => length(png_files),
        "provenance_sidecar_count" => length(sidecar_files),
        "min_png_bytes" => min_png_bytes,
    )
end

function _remove_plot_sidecars(fig_dir::String)
    for path in sort(readdir(fig_dir; join=true))
        isdir(path) || continue
        for sub in sort(readdir(path; join=true))
            isfile(sub) || continue
            endswith(lowercase(sub), ".png.provenance.json") && rm(sub; force=true)
        end
    end
end

function _refresh_case_manifest(case_dir::String, fig_dir::String)
    manifest_path = joinpath(case_dir, "manifest.json")
    isfile(manifest_path) || return

    manifest = JSON3.read(read(manifest_path, String), Dict{String,Any})
    files = get(manifest, "result_files", Any[])
    hashes = Dict{String,Any}(get(manifest, "hashes", Dict{String,Any}()))
    for name in files
        name_string = String(name)
        name_string == "manifest.json" && continue
        path = joinpath(case_dir, name_string)
        isfile(path) && (hashes[name_string] = _sha256_file(path))
    end

    plot_manifest_path = joinpath(fig_dir, "plot_manifest.json")
    manifest["hashes"] = hashes
    manifest["figure_summary"] = _plot_output_summary(fig_dir)
    if isfile(plot_manifest_path)
        manifest["figure_hashes"] = Dict(
            "plot_manifest.json" => _sha256_file(plot_manifest_path),
        )
    end

    open(manifest_path, "w") do io
        JSON3.pretty(io, manifest)
        println(io)
    end
end

function _append_readme(case_dir::String, fig_dir::String, layout)
    readme_path = joinpath(case_dir, "README.md")
    isfile(readme_path) || return

    files = String[]
    for path in sort(readdir(fig_dir; join=true))
        isdir(path) || continue
        for sub in sort(readdir(path; join=true))
            isfile(sub) || continue
            endswith(lowercase(sub), ".png") || continue
            push!(files, replace(relpath(sub, PROJECT_ROOT), '\\' => '/'))
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
        println(io, "- `tau_u`, `tau_d`, `tau_s`, `eta`, `sigma`, `zeta`, `eta_over_s`, `sigma_over_T` all plot against `xi`")
        println(io, "- panel rule: `$(layout.panel_desc)`")
        println(io, "- line rule: `$(layout.series_desc)`")
        println(io, "- plot manifest: `$(replace(relpath(joinpath(fig_dir, "plot_manifest.json"), PROJECT_ROOT), '\\' => '/'))`")
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
    mode = _detect_mode(opts.case_dir, opts.csv_path)
    layout = _plot_layout(mode)
    if opts.overwrite && isdir(opts.fig_dir)
        rm(opts.fig_dir; recursive=true, force=true)
    end
    mkpath(opts.fig_dir)

    for y in ("tau_u", "tau_d", "tau_s", "eta", "sigma", "zeta", "eta_over_s", "sigma_over_T")
        _run_plot(py, opts.csv_path, opts.fig_dir, y; split=layout.split, group=layout.group)
    end

    _remove_plot_sidecars(opts.fig_dir)
    _write_plot_manifest(opts.fig_dir)
    _append_readme(opts.case_dir, opts.fig_dir, layout)
    _refresh_case_manifest(opts.case_dir, opts.fig_dir)
    println("Phase-guided transport plotting finished. Output: $(opts.fig_dir)")
end

function main(args::Vector{String}=copy(ARGS))
    opts = parse_args(args)
    run_phase_guided_plots(opts)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
