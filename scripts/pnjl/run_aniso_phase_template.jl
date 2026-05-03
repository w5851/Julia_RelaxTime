#!/usr/bin/env julia

"""
各向异性 PNJL + 相图产出可复现实验模板。

用途：
- 统一执行 T-μ 扫描（各 xi）
- 统一执行相结构计算（boundary/CEP/spinodal/crossover）
- 可选调用 Python 绘图脚本生成相图

示例：
  julia --project=. scripts/pnjl/run_aniso_phase_template.jl --profile=smoke --xi-values=0.0,0.2
  julia --project=. scripts/pnjl/run_aniso_phase_template.jl --profile=full --xi-values=0.0,0.2,0.4 --skip-plot
"""

using Dates
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SCRIPT_TMU = joinpath(PROJECT_ROOT, "scripts", "models", "run_unified_scan.jl")
const SCRIPT_PHASE = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")
const SCRIPT_PLOT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "plot_phase_diagram.py")

function _parse_xi_values(raw::String)
    vals = Float64[]
    for token in split(raw, ',')
        t = strip(token)
        isempty(t) && continue
        push!(vals, parse(Float64, t))
    end
    isempty(vals) && error("xi-values 不能为空")
    return vals
end

function _profile_config(profile::String)
    p = lowercase(profile)
    if p == "smoke"
        return (
            tmu=(T_min=100, T_max=140, T_step=20, mu_min=0, mu_max=240, mu_step=60, p_num=12, t_num=4),
            phase=(T_min=80, T_max=160, T_step=20, rho_min=0.0, rho_max=2.0, rho_step=0.2),
        )
    elseif p == "full"
        return (
            tmu=(T_min=50, T_max=200, T_step=10, mu_min=0, mu_max=400, mu_step=10, p_num=24, t_num=8),
            phase=(T_min=30, T_max=350, T_step=10, rho_min=0.0, rho_max=4.0, rho_step=0.05),
        )
    else
        error("未知 profile: $profile（支持 smoke/full）")
    end
end

function _help()
    println("""
各向异性 PNJL + 相图产出可复现实验模板

用法：
  julia --project=. scripts/pnjl/run_aniso_phase_template.jl [options]

选项：
  --profile=smoke|full    预置网格规模（默认 smoke）
  --xi-values=0.0,0.2     xi 列表（默认 0.0,0.2）
  --tag=YYYYmmdd_HHMMSS   运行标签（默认自动时间戳）
  --skip-plot             跳过 python 相图绘制
  --python=python         Python 命令名（默认 python）
  --help                  显示帮助
""")
end

function _range_values(min_v::Real, max_v::Real, step_v::Real)
    step_v > 0 || error("step must be positive")
    max_v >= min_v || error("max must be >= min")
    values = collect(Float64(min_v):Float64(step_v):Float64(max_v))
    isempty(values) && error("range cannot be empty")
    return values
end

function _csv_values(values::AbstractVector{<:Real})
    return join((@sprintf("%.6g", v) for v in values), ",")
end

function _run_cmd(args::Vector{String}; title::String)
    println("\n>>> ", title)
    println(join(args, " "))
    run(Cmd(args; dir=PROJECT_ROOT))
end

function main(args=ARGS)
    profile = "smoke"
    xi_values = [0.0, 0.2]
    tag = Dates.format(now(), "yyyymmdd_HHMMSS")
    skip_plot = false
    python_cmd = "python"

    for arg in args
        if arg == "--help" || arg == "-h"
            _help()
            return
        elseif startswith(arg, "--profile=")
            profile = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--xi-values=")
            xi_values = _parse_xi_values(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--tag=")
            tag = split(arg, "=", limit=2)[2]
        elseif arg == "--skip-plot"
            skip_plot = true
        elseif startswith(arg, "--python=")
            python_cmd = split(arg, "=", limit=2)[2]
        else
            @warn "未知参数" arg
        end
    end

    cfg = _profile_config(profile)

    out_root = joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl", "experiment_templates", "aniso_phase", tag)
    out_scan_dir = joinpath(out_root, "scan")
    out_phase_dir = joinpath(out_root, "phase")
    out_fig_dir = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "pnjl", "experiment_templates", "aniso_phase", tag)
    mkpath(out_scan_dir)
    mkpath(out_phase_dir)
    mkpath(out_fig_dir)

    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())

    println("="^68)
    println("各向异性 PNJL 相图模板执行")
    println("profile = ", profile)
    println("xi-values = ", xi_values)
    println("tag = ", tag)
    println("out scan = ", out_scan_dir)
    println("out phase = ", out_phase_dir)
    println("out fig = ", out_fig_dir)
    println("="^68)

    for xi in xi_values
        xi_str = @sprintf("%.3f", xi)

        tmu_out = joinpath(out_scan_dir, "tmu_scan_xi$(xi_str).csv")
        T_values = _csv_values(_range_values(cfg.tmu.T_min, cfg.tmu.T_max, cfg.tmu.T_step))
        mu_values = _csv_values(_range_values(cfg.tmu.mu_min, cfg.tmu.mu_max, cfg.tmu.mu_step))
        _run_cmd([
            julia_bin,
            "--project=.",
            SCRIPT_TMU,
            "scan",
            "tmu",
            "--model_kind=pnjl_aniso",
            "--T_values=$(T_values)",
            "--mu_values=$(mu_values)",
            "--xi_values=$(xi)",
            "--p_num=$(cfg.tmu.p_num)",
            "--t_num=$(cfg.tmu.t_num)",
            "--output_path=$(tmu_out)",
            "--overwrite",
        ]; title="T-μ 扫描 xi=$(xi)")

        _run_cmd([
            julia_bin,
            "--project=.",
            SCRIPT_PHASE,
            "--xi=$(xi)",
            "--T_min=$(cfg.phase.T_min)",
            "--T_max=$(cfg.phase.T_max)",
            "--T_step=$(cfg.phase.T_step)",
            "--rho_min=$(cfg.phase.rho_min)",
            "--rho_max=$(cfg.phase.rho_max)",
            "--rho_step=$(cfg.phase.rho_step)",
            "--output_dir=$(out_phase_dir)",
        ]; title="相结构计算 xi=$(xi)")
    end

    if !skip_plot
        _run_cmd([
            python_cmd,
            SCRIPT_PLOT,
            "--boundary", joinpath(out_phase_dir, "boundary.csv"),
            "--cep", joinpath(out_phase_dir, "cep.csv"),
            "--spinodal", joinpath(out_phase_dir, "spinodals.csv"),
            "--crossover", joinpath(out_phase_dir, "crossover.csv"),
            "--output-dir", out_fig_dir,
            "--no-show",
        ]; title="相图绘制")
    end

    manifest = joinpath(out_root, "manifest.txt")
    open(manifest, "w") do io
        println(io, "tag=", tag)
        println(io, "profile=", profile)
        println(io, "xi_values=", join(xi_values, ","))
        println(io, "scan_dir=", out_scan_dir)
        println(io, "phase_dir=", out_phase_dir)
        println(io, "fig_dir=", out_fig_dir)
        println(io, "generated_at=", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    end

    println("\n完成。清单文件：", manifest)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
