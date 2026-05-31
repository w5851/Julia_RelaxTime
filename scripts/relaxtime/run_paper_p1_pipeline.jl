#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

Base.@kwdef mutable struct PaperP1Options
    stages::Vector{Symbol} = Symbol[:assets]
    result_dir::String = joinpath("data", "outputs", "results", "relaxtime", "paper_p1", "orchestrated")
    figure_dir::String = joinpath("data", "outputs", "figures", "relaxtime", "paper_p1", "orchestrated")
    xi_values::Vector{Float64} = [-0.3, 0.0, 0.3]
    mott_mub_values::Vector{Float64} = [0.0, 200.0, 400.0, 600.0, 800.0]
    mott_tmin::Float64 = 100.0
    mott_tmax::Float64 = 300.0
    mott_tstep::Float64 = 5.0
    mott_slice_plan::Union{Nothing,String} = nothing
    mott_p_num::Int = 12
    mott_t_num::Int = 6
    mott_max_iter::Int = 40
    mott_equilibrium_branch_mode::Symbol = :stable
    mott_grid_csv::Union{Nothing,String} = nothing
    sigma_values::Vector{Float64} = [20.0, 30.0, 45.0]
    isentropic_xi_values::Vector{Float64} = [0.0]
    isentropic_t_values::Vector{Float64} = collect(100.0:10.0:300.0)
    isentropic_csvs::Vector{String} = String[]
    phase_dirs::Vector{String} = String[]
    phase_output_root::String = joinpath("data", "processed", "pnjl", "phase_diagram", "paper_p1_orchestrated")
    phase_reference_root::String = joinpath("data", "reference", "pnjl", "paper_p1_orchestrated")
    use_phase_reference::Bool = false
    phase_tag::String = "paper_p1"
    phase_tmin::Float64 = 30.0
    phase_tmax::Float64 = 240.0
    phase_tstep::Float64 = 5.0
    phase_rhomin::Float64 = 0.0
    phase_rhomax::Float64 = 4.0
    phase_rhostep::Float64 = 0.05
    phase_p_num::Int = 24
    phase_t_num::Int = 8
    phase_iterations::Int = 80
    phase_mu_scale::Float64 = 3.0
    formats::String = "png,pdf"
    dry_run::Bool = false
    overwrite::Bool = false
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_paper_p1_pipeline.jl [options]")
    println()
    println("Stages are independent and can be combined:")
    println("  --stage=mott       run Mott mu_B slices and combine CSV")
    println("  --stage=phase      run dense phase reference for selected xi")
    println("  --stage=isentropic run fixed-sigma meson-mass paths")
    println("  --stage=assets     build CSV/figure assets from existing or produced inputs")
    println("  --stage=all        expands to mott,phase,isentropic,assets")
    println()
    println("Common options:")
    println("  --result-dir <path>         base result directory")
    println("  --figure-dir <path>         figure asset directory")
    println("  --xi-list <csv>             xi values")
    println("  --dry-run                   print commands without running")
    println("  --overwrite                 pass overwrite where supported")
    println()
    println("Input override options for assets stage:")
    println("  --mott-grid-csv <path>")
    println("  --phase-dir <path>          repeatable")
    println("  --phase-reference-root <path>")
    println("  --phase-tag <name>")
    println("  --isentropic-csv <path>     repeatable")
    println()
    println("Mott options:")
    println("  --mott-mub-list <csv>")
    println("  --mott-tmin/--mott-tmax/--mott-tstep <MeV>")
    println("  --mott-slice-plan <csv>       CSV with muB_MeV,T_min_MeV,T_max_MeV,T_step_MeV")
    println("  --mott-p-num/--mott-t-num/--mott-max-iter <int>")
    println("  --mott-equilibrium-branch-mode <stable|continuation>")
    println()
    println("Phase options:")
    println("  --phase-output-root <path>")
    println("  --phase-reference-root <path>")
    println("  --phase-tag <name>")
    println("  --phase-tmin/--phase-tmax/--phase-tstep <MeV>")
    println("  --phase-rho-min/--phase-rho-max/--phase-rho-step <value>")
    println("  --phase-p-num/--phase-t-num/--phase-iterations <int>")
    println()
    println("Isentropic options:")
    println("  --sigma-list <csv>")
    println("  --isentropic-xi-list <csv>   independent xi list for fixed-sigma paths")
    println("  --isentropic-T-list <csv>")
end

function parse_float_list(raw::AbstractString)
    vals = Float64[]
    for token in split(raw, ',')
        s = strip(token)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("empty numeric list: $(raw)"))
    return vals
end

function parse_stage_list(raw::AbstractString)
    stages = Symbol[]
    for token in split(raw, ',')
        s = Symbol(lowercase(strip(token)))
        s === :all && return Symbol[:mott, :phase, :isentropic, :assets]
        s in (:mott, :phase, :isentropic, :assets) || throw(ArgumentError("unsupported stage: $(token)"))
        push!(stages, s)
    end
    return stages
end

function parse_args(args::Vector{String})
    opts = PaperP1Options()
    opts.stages = Symbol[]
    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $(arg)"))
            i += 1
            return args[i]
        end

        if arg == "--stage"
            append!(opts.stages, parse_stage_list(require_value()))
        elseif startswith(arg, "--stage=")
            append!(opts.stages, parse_stage_list(split(arg, "="; limit=2)[2]))
        elseif arg == "--result-dir"
            opts.result_dir = require_value()
        elseif arg == "--figure-dir" || arg == "--out-dir"
            opts.figure_dir = require_value()
        elseif arg == "--xi-list"
            opts.xi_values = parse_float_list(require_value())
        elseif arg == "--mott-grid-csv"
            opts.mott_grid_csv = require_value()
        elseif arg == "--phase-dir"
            push!(opts.phase_dirs, require_value())
        elseif arg == "--isentropic-csv"
            push!(opts.isentropic_csvs, require_value())
        elseif arg == "--mott-mub-list"
            opts.mott_mub_values = parse_float_list(require_value())
        elseif arg == "--mott-tmin"
            opts.mott_tmin = parse(Float64, require_value())
        elseif arg == "--mott-tmax"
            opts.mott_tmax = parse(Float64, require_value())
        elseif arg == "--mott-tstep"
            opts.mott_tstep = parse(Float64, require_value())
        elseif arg == "--mott-slice-plan"
            opts.mott_slice_plan = require_value()
        elseif arg == "--mott-p-num"
            opts.mott_p_num = parse(Int, require_value())
        elseif arg == "--mott-t-num"
            opts.mott_t_num = parse(Int, require_value())
        elseif arg == "--mott-max-iter"
            opts.mott_max_iter = parse(Int, require_value())
        elseif arg == "--mott-equilibrium-branch-mode"
            mode = Symbol(lowercase(require_value()))
            mode in (:stable, :continuation) || throw(ArgumentError("mott equilibrium branch mode must be stable or continuation"))
            opts.mott_equilibrium_branch_mode = mode
        elseif arg == "--phase-output-root"
            opts.phase_output_root = require_value()
        elseif arg == "--phase-reference-root"
            opts.phase_reference_root = require_value()
            opts.use_phase_reference = true
        elseif arg == "--phase-tag"
            opts.phase_tag = require_value()
            opts.use_phase_reference = true
        elseif arg == "--phase-tmin"
            opts.phase_tmin = parse(Float64, require_value())
        elseif arg == "--phase-tmax"
            opts.phase_tmax = parse(Float64, require_value())
        elseif arg == "--phase-tstep"
            opts.phase_tstep = parse(Float64, require_value())
        elseif arg == "--phase-rho-min" || arg == "--phase-rhomin"
            opts.phase_rhomin = parse(Float64, require_value())
        elseif arg == "--phase-rho-max" || arg == "--phase-rhomax"
            opts.phase_rhomax = parse(Float64, require_value())
        elseif arg == "--phase-rho-step" || arg == "--phase-rhostep"
            opts.phase_rhostep = parse(Float64, require_value())
        elseif arg == "--phase-p-num"
            opts.phase_p_num = parse(Int, require_value())
        elseif arg == "--phase-t-num"
            opts.phase_t_num = parse(Int, require_value())
        elseif arg == "--phase-iterations"
            opts.phase_iterations = parse(Int, require_value())
        elseif arg == "--sigma-list"
            opts.sigma_values = parse_float_list(require_value())
        elseif arg == "--isentropic-xi-list"
            opts.isentropic_xi_values = parse_float_list(require_value())
        elseif arg == "--isentropic-T-list"
            opts.isentropic_t_values = parse_float_list(require_value())
        elseif arg == "--formats"
            opts.formats = require_value()
        elseif arg == "--dry-run"
            opts.dry_run = true
        elseif arg == "--overwrite"
            opts.overwrite = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $(arg)"))
        end
        i += 1
    end
    isempty(opts.stages) && (opts.stages = Symbol[:assets])
    opts.stages = unique(opts.stages)
    return opts
end

function _run_cmd(cmd::Cmd, opts::PaperP1Options)
    println(join(cmd.exec, " "))
    opts.dry_run || run(cmd)
end

function _julia_script_cmd(script::String, args::Vector{String})
    return Cmd(vcat([Base.julia_cmd().exec..., "--project=$(PROJECT_ROOT)", script], args))
end

function _python_cmd(args::Vector{String})
    python = get(ENV, "PYTHON", "")
    if isempty(python)
        py = Sys.which("python")
        py === nothing && (py = Sys.which("python3"))
        py === nothing && throw(ArgumentError("python executable not found"))
        python = py
    end
    return Cmd(vcat([python], args))
end

function _float_tag(x::Float64)
    return replace(replace(string(round(Int, x)), "-" => "m"), "." => "p")
end

function _parse_csv_line(line::AbstractString)
    return [strip(String(part)) for part in split(String(line), ',')]
end

function load_mott_slice_plan(path::String)
    isfile(path) || throw(ArgumentError("mott slice plan not found: $(path)"))
    rows = NamedTuple[]
    header = String[]
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "#") && continue
        if isempty(header)
            header = _parse_csv_line(line)
            required = ("muB_MeV", "T_min_MeV", "T_max_MeV", "T_step_MeV")
            required_label = join(required, ",")
            all(name -> name in header, required) ||
                throw(ArgumentError("mott slice plan must include columns: $(required_label)"))
            continue
        end
        vals = _parse_csv_line(line)
        row = Dict{String,String}()
        for (idx, key) in enumerate(header)
            row[key] = idx <= length(vals) ? vals[idx] : ""
        end
        push!(rows, (
            muB=Float64(parse(Float64, row["muB_MeV"])),
            T_min=Float64(parse(Float64, row["T_min_MeV"])),
            T_max=Float64(parse(Float64, row["T_max_MeV"])),
            T_step=Float64(parse(Float64, row["T_step_MeV"])),
        ))
    end
    isempty(rows) && throw(ArgumentError("mott slice plan has no data rows: $(path)"))
    return rows
end

function mott_slices(opts::PaperP1Options)
    if opts.mott_slice_plan !== nothing
        return load_mott_slice_plan(opts.mott_slice_plan)
    end
    return [(muB=Float64(muB), T_min=opts.mott_tmin, T_max=opts.mott_tmax, T_step=opts.mott_tstep) for muB in opts.mott_mub_values]
end

function write_mott_config(path::String, opts::PaperP1Options, slice)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "schema_version = \"v1\"")
        println(io, "profile_name = \"paper_p1_mott_muB$(_float_tag(slice.muB))\"")
        println(io)
        println(io, "[scan.mott_phase]")
        println(io, "muB_MeV = $(slice.muB)")
        println(io, "xi_list = [", join(opts.xi_values, ", "), "]")
        println(io, "T_min_MeV = $(slice.T_min)")
        println(io, "T_max_MeV = $(slice.T_max)")
        println(io, "T_step_MeV = $(slice.T_step)")
        println(io, "p_num = $(opts.mott_p_num)")
        println(io, "t_num = $(opts.mott_t_num)")
        println(io, "max_iter = $(opts.mott_max_iter)")
        println(io, "equilibrium_branch_mode = \"$(opts.mott_equilibrium_branch_mode)\"")
        println(io, "resume = true")
        println(io, "overwrite = false")
    end
end

function combine_mott_csv(paths::Vector{String}, output_path::String)
    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)
    header_written = false
    metadata_written = false
    open(output_path, "w") do out
        for path in paths
            for line in eachline(path)
                s = strip(line)
                isempty(s) && continue
                if startswith(s, "#")
                    if !metadata_written
                        println(out, line)
                    end
                    continue
                end
                if startswith(s, "run_id,") || startswith(s, "T_MeV,")
                    if !header_written
                        println(out, line)
                        header_written = true
                        metadata_written = true
                    end
                    continue
                end
                println(out, line)
            end
        end
    end
end

function _selector_metadata(mode::Symbol)
    if mode === :stable
        return (
            equilibrium_branch_mode="stable",
            equilibrium_selector_policy="pressure_max_under_constraints",
            equilibrium_selector_tiebreak="residual_norm_then_seed_index",
        )
    end
    return (
        equilibrium_branch_mode=String(mode),
        equilibrium_selector_policy="continuation_seed",
        equilibrium_selector_tiebreak="not_applicable",
    )
end

function _json_string(value::AbstractString)
    escaped = replace(replace(String(value), "\\" => "\\\\"), "\"" => "\\\"")
    return "\"" * escaped * "\""
end

function write_mott_combined_manifest(path::String, opts::PaperP1Options, slices, csvs::Vector{String}, output_csv::String)
    mkpath(dirname(path))
    selector = _selector_metadata(opts.mott_equilibrium_branch_mode)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"schema_version\": \"paper_p1_mott_combined_v1\",")
        println(io, "  \"produced_by\": \"scripts/relaxtime/run_paper_p1_pipeline.jl\",")
        println(io, "  \"output_csv\": ", _json_string(output_csv), ",")
        println(io, "  \"scan\": {")
        println(io, "    \"mott_phase\": {")
        println(io, "      \"xi_list\": [$(join(opts.xi_values, ", "))],")
        println(io, "      \"equilibrium_branch_mode\": \"$(selector.equilibrium_branch_mode)\",")
        println(io, "      \"equilibrium_selector_policy\": \"$(selector.equilibrium_selector_policy)\",")
        println(io, "      \"equilibrium_selector_tiebreak\": \"$(selector.equilibrium_selector_tiebreak)\",")
        println(io, "      \"p_num\": $(opts.mott_p_num),")
        println(io, "      \"t_num\": $(opts.mott_t_num),")
        println(io, "      \"max_iter\": $(opts.mott_max_iter)")
        println(io, "    }")
        println(io, "  },")
        println(io, "  \"slices\": [")
        for (idx, slice) in enumerate(slices)
            comma = idx == length(slices) ? "" : ","
            println(io, "    {\"muB_MeV\": $(slice.muB), \"T_min_MeV\": $(slice.T_min), \"T_max_MeV\": $(slice.T_max), \"T_step_MeV\": $(slice.T_step), \"source_csv\": $(_json_string(csvs[idx]))}$(comma)")
        end
        println(io, "  ]")
        println(io, "}")
    end
end

function run_mott_stage!(opts::PaperP1Options)
    config_dir = joinpath(opts.result_dir, "configs")
    csvs = String[]
    slices = mott_slices(opts)
    if opts.mott_slice_plan !== nothing
        println("mott slice plan: $(opts.mott_slice_plan)")
    end
    for slice in slices
        tag = lpad(string(round(Int, slice.muB)), 4, '0')
        cfg_path = joinpath(config_dir, "mott_muB$(tag).toml")
        outdir = joinpath(opts.result_dir, "mott_muB$(tag)")
        write_mott_config(cfg_path, opts, slice)
        cmd = _julia_script_cmd(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_mott_phase_scan.jl"),
            ["--config", cfg_path, "--outdir", outdir, opts.overwrite ? "--overwrite" : "--resume"])
        _run_cmd(cmd, opts)
        push!(csvs, joinpath(outdir, "mott_phase_scan.csv"))
    end
    output_csv = joinpath(opts.result_dir, "mott_grid_combined.csv")
    if opts.dry_run
        println("combine_mott_csv -> $(output_csv)")
        println("write_mott_combined_manifest -> $(joinpath(opts.result_dir, "mott_grid_combined_manifest.json"))")
    else
        combine_mott_csv(csvs, output_csv)
        write_mott_combined_manifest(joinpath(opts.result_dir, "mott_grid_combined_manifest.json"), opts, slices, csvs, output_csv)
    end
    opts.mott_grid_csv = output_csv
end

function run_phase_stage!(opts::PaperP1Options)
    cmd = _julia_script_cmd(joinpath(PROJECT_ROOT, "scripts", "pnjl", "build_dense_phase_reference.jl"),
        [
            "--xi-list", join(opts.xi_values, ","),
            "--T-min", string(opts.phase_tmin),
            "--T-max", string(opts.phase_tmax),
            "--T-step", string(opts.phase_tstep),
            "--rho-min", string(opts.phase_rhomin),
            "--rho-max", string(opts.phase_rhomax),
            "--rho-step", string(opts.phase_rhostep),
            "--p-num", string(opts.phase_p_num),
            "--t-num", string(opts.phase_t_num),
            "--iterations", string(opts.phase_iterations),
            "--tag", opts.phase_tag,
            "--output-root", opts.phase_output_root,
            "--reference-root", opts.phase_reference_root,
            opts.overwrite ? "--overwrite" : "",
        ])
    _run_cmd(Cmd(filter(!isempty, cmd.exec)), opts)
    opts.phase_dirs = [joinpath(opts.phase_output_root, "xi_" * replace(replace(@sprintf("%.3f", xi), "." => "p"), "-" => "m")) for xi in opts.xi_values]
    opts.use_phase_reference = true
end

function run_isentropic_stage!(opts::PaperP1Options)
    opts.isentropic_csvs = String[]
    t_list = join(opts.isentropic_t_values, ",")
    for xi in opts.isentropic_xi_values, sigma in opts.sigma_values
        xi_tag = replace(replace(@sprintf("%.3f", xi), "." => "p"), "-" => "m")
        sigma_tag = replace(replace(string(sigma), "." => "p"), "-" => "m")
        output = joinpath(opts.result_dir, "isentropic_xi$(xi_tag)_sigma$(sigma_tag).csv")
        cmd = _julia_script_cmd(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_isentropic_meson_mass_scan.jl"),
            ["--sigma-target", string(sigma), "--T-list", t_list, "--xi", string(xi), "--output", output, opts.overwrite ? "--overwrite" : ""])
        _run_cmd(Cmd(filter(!isempty, cmd.exec)), opts)
        push!(opts.isentropic_csvs, output)
    end
end

function run_assets_stage!(opts::PaperP1Options)
    opts.mott_grid_csv === nothing && throw(ArgumentError("assets stage requires --mott-grid-csv or stage=mott"))
    args = String[
        joinpath(PROJECT_ROOT, "scripts", "relaxtime", "build_paper_p1_figure_assets.py"),
        "--mott-grid-csv", opts.mott_grid_csv,
        "--out-dir", opts.figure_dir,
        "--phase-mu-scale", string(opts.phase_mu_scale),
        "--formats", opts.formats,
    ]
    for path in opts.isentropic_csvs
        append!(args, ["--isentropic-csv", path])
    end
    if opts.use_phase_reference
        append!(args, ["--phase-reference-root", opts.phase_reference_root, "--phase-reference-tag", opts.phase_tag])
    else
        for path in opts.phase_dirs
            append!(args, ["--phase-dir", path])
        end
    end
    _run_cmd(_python_cmd(args), opts)
end

function main(args::Vector{String})
    opts = parse_args(args)
    :mott in opts.stages && run_mott_stage!(opts)
    :phase in opts.stages && run_phase_stage!(opts)
    :isentropic in opts.stages && run_isentropic_stage!(opts)
    :assets in opts.stages && run_assets_stage!(opts)
    return 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main(ARGS))
end
