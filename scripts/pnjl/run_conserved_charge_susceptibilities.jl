#!/usr/bin/env julia
#=
PNJL 守恒荷广义磁化率与累积量统一入口脚本

目标：为用户提供单一稳定入口，统一输出
- chi_BQS(...)
- cumulant_BQS(...)
- baryon_Ssigma(...)
- baryon_kappa_sigma2(...)

当前定位：
- `scripts/pnjl/run_*.jl` 作为稳定用户入口
- `scripts/dev/` 保留开发期导出/一次性工具
- `scripts/analysis/` 保留后处理分析脚本
=#

const HELP_TEXT = """
PNJL 守恒荷广义磁化率与累积量统一入口

用法：
    julia --project=. scripts/pnjl/run_conserved_charge_susceptibilities.jl [options]

常用选项：
    --T=150            温度 (MeV)
    --muB=0            重子化学势 (MeV)
    --muQ=0            电荷化学势 (MeV)
    --muS=0            奇异化学势 (MeV)
    --V=64             体积 (fm^3)，用于累积量
    --scan=none        扫描轴: none|T|muB|muQ|muS
    --scan_min=0       扫描起点 (MeV)
    --scan_max=0       扫描终点 (MeV)
    --scan_step=10     扫描步长 (MeV)
    --xi=0.0           各向异性参数
    --p_num=24         动量积分节点数
    --t_num=8          角度积分节点数
    --output=PATH      可选：将结果写入 CSV 文件
    --help             显示帮助

默认输出：
    - chi2_B, chi2_Q, chi2_S
    - chi11_BQ, chi11_BS, chi11_QS
    - chi1_B, chi3_B, chi4_B
    - C200, C020, C002, C110, C101, C011
    - baryon Ssigma, baryon kappa_sigma2
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using Main.Constants_PNJL: ħc_MeV_fm
const PNJL = Main.Models.pnjl_module()

Base.@kwdef struct SusceptibilityCLIConfig
    T_MeV::Float64 = 150.0
    muB_MeV::Float64 = 0.0
    muQ_MeV::Float64 = 0.0
    muS_MeV::Float64 = 0.0
    V_fm3::Float64 = 64.0
    scan_axis::Symbol = :none
    scan_min_MeV::Float64 = 0.0
    scan_max_MeV::Float64 = 0.0
    scan_step_MeV::Float64 = 10.0
    xi::Float64 = 0.0
    p_num::Int = 24
    t_num::Int = 8
    output::Union{Nothing, String} = nothing
end

function parse_args(args::Vector{String})
    T_MeV = 150.0
    muB_MeV = 0.0
    muQ_MeV = 0.0
    muS_MeV = 0.0
    V_fm3 = 64.0
    scan_axis = :none
    scan_min_MeV = 0.0
    scan_max_MeV = 0.0
    scan_step_MeV = 10.0
    xi = 0.0
    p_num = 24
    t_num = 8
    output = nothing

    for arg in args
        if arg == "--help" || arg == "-h"
            println(HELP_TEXT)
            exit(0)
        elseif startswith(arg, "--T=")
            T_MeV = parse(Float64, arg[5:end])
        elseif startswith(arg, "--muB=")
            muB_MeV = parse(Float64, arg[7:end])
        elseif startswith(arg, "--muQ=")
            muQ_MeV = parse(Float64, arg[7:end])
        elseif startswith(arg, "--muS=")
            muS_MeV = parse(Float64, arg[7:end])
        elseif startswith(arg, "--V=")
            V_fm3 = parse(Float64, arg[5:end])
        elseif startswith(arg, "--scan=")
            value = Symbol(arg[8:end])
            value in (:none, :T, :muB, :muQ, :muS) || error("invalid --scan value: $(value)")
            scan_axis = value
        elseif startswith(arg, "--scan_min=")
            scan_min_MeV = parse(Float64, arg[12:end])
        elseif startswith(arg, "--scan_max=")
            scan_max_MeV = parse(Float64, arg[12:end])
        elseif startswith(arg, "--scan_step=")
            scan_step_MeV = parse(Float64, arg[13:end])
        elseif startswith(arg, "--xi=")
            xi = parse(Float64, arg[6:end])
        elseif startswith(arg, "--p_num=")
            p_num = parse(Int, arg[9:end])
        elseif startswith(arg, "--t_num=")
            t_num = parse(Int, arg[9:end])
        elseif startswith(arg, "--output=")
            output = arg[10:end]
        else
            error("unknown argument: $arg")
        end
    end

    return SusceptibilityCLIConfig(
        T_MeV=T_MeV,
        muB_MeV=muB_MeV,
        muQ_MeV=muQ_MeV,
        muS_MeV=muS_MeV,
        V_fm3=V_fm3,
        scan_axis=scan_axis,
        scan_min_MeV=scan_min_MeV,
        scan_max_MeV=scan_max_MeV,
        scan_step_MeV=scan_step_MeV,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        output=output,
    )
end

function default_output_path(cfg::SusceptibilityCLIConfig)
    base = joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl", "susceptibilities")
    if cfg.scan_axis === :none
        return joinpath(base, "single_point.csv")
    end
    return joinpath(base, "scan", "$(cfg.scan_axis)_scan.csv")
end

function compute_summary(cfg::SusceptibilityCLIConfig)
    T_fm = cfg.T_MeV / ħc_MeV_fm
    muB_fm = cfg.muB_MeV / ħc_MeV_fm
    muQ_fm = cfg.muQ_MeV / ħc_MeV_fm
    muS_fm = cfg.muS_MeV / ħc_MeV_fm

    kwargs = (; xi=cfg.xi, p_num=cfg.p_num, t_num=cfg.t_num)

    chi1B = PNJL.chi1_B(T_fm, muB_fm; kwargs...)
    chi2B = PNJL.chi2_B(T_fm, muB_fm; kwargs...)
    chi3B = PNJL.chi3_B(T_fm, muB_fm; kwargs...)
    chi4B = PNJL.chi4_B(T_fm, muB_fm; kwargs...)

    chi2Q = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 2, 0), kwargs...)
    chi2S = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 0, 2), kwargs...)
    chi11BQ = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), kwargs...)
    chi11BS = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 0, 1), kwargs...)
    chi11QS = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 1, 1), kwargs...)

    c200 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, cfg.V_fm3; orders=(2, 0, 0), kwargs...)
    c020 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, cfg.V_fm3; orders=(0, 2, 0), kwargs...)
    c002 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, cfg.V_fm3; orders=(0, 0, 2), kwargs...)
    c110 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, cfg.V_fm3; orders=(1, 1, 0), kwargs...)
    c101 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, cfg.V_fm3; orders=(1, 0, 1), kwargs...)
    c011 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, cfg.V_fm3; orders=(0, 1, 1), kwargs...)

    return (
        T_MeV=cfg.T_MeV,
        muB_MeV=cfg.muB_MeV,
        muQ_MeV=cfg.muQ_MeV,
        muS_MeV=cfg.muS_MeV,
        V_fm3=cfg.V_fm3,
        chi1_B=chi1B,
        chi2_B=chi2B,
        chi3_B=chi3B,
        chi4_B=chi4B,
        chi2_Q=chi2Q,
        chi2_S=chi2S,
        chi11_BQ=chi11BQ,
        chi11_BS=chi11BS,
        chi11_QS=chi11QS,
        C200=c200,
        C020=c020,
        C002=c002,
        C110=c110,
        C101=c101,
        C011=c011,
        Ssigma=PNJL.baryon_Ssigma(T_fm, muB_fm; kwargs...),
        kappa_sigma2=PNJL.baryon_kappa_sigma2(T_fm, muB_fm; kwargs...),
    )
end

function print_summary(result)
    println("=" ^ 68)
    println("PNJL conserved-charge susceptibilities")
    println("=" ^ 68)
    @printf("T = %.3f MeV, muB = %.3f MeV, muQ = %.3f MeV, muS = %.3f MeV\n", result.T_MeV, result.muB_MeV, result.muQ_MeV, result.muS_MeV)
    @printf("V = %.3f fm^3\n", result.V_fm3)
    println()
    @printf("chi1_B   = %.10e\n", result.chi1_B)
    @printf("chi2_B   = %.10e\n", result.chi2_B)
    @printf("chi3_B   = %.10e\n", result.chi3_B)
    @printf("chi4_B   = %.10e\n", result.chi4_B)
    println()
    @printf("chi2_Q   = %.10e\n", result.chi2_Q)
    @printf("chi2_S   = %.10e\n", result.chi2_S)
    @printf("chi11_BQ = %.10e\n", result.chi11_BQ)
    @printf("chi11_BS = %.10e\n", result.chi11_BS)
    @printf("chi11_QS = %.10e\n", result.chi11_QS)
    println()
    @printf("C200     = %.10e\n", result.C200)
    @printf("C020     = %.10e\n", result.C020)
    @printf("C002     = %.10e\n", result.C002)
    @printf("C110     = %.10e\n", result.C110)
    @printf("C101     = %.10e\n", result.C101)
    @printf("C011     = %.10e\n", result.C011)
    println()
    @printf("Ssigma        = %.10e\n", result.Ssigma)
    @printf("kappa_sigma2  = %.10e\n", result.kappa_sigma2)
end

function maybe_write_csv(result, output_path::String)
    mkpath(dirname(output_path))
    names = propertynames(result)
    open(output_path, "w") do io
        println(io, join(string.(names), ','))
        println(io, join((@sprintf("%.16e", Float64(getproperty(result, k))) for k in names), ','))
    end
end

function write_scan_csv(results, output_path::String)
    mkpath(dirname(output_path))
    names = propertynames(first(results))
    open(output_path, "w") do io
        println(io, join(string.(names), ','))
        for result in results
            row = String[]
            for k in names
                value = getproperty(result, k)
                if value isa Real
                    push!(row, @sprintf("%.16e", Float64(value)))
                else
                    push!(row, string(value))
                end
            end
            println(io, join(row, ','))
        end
    end
end

function with_scan_value(cfg::SusceptibilityCLIConfig, axis::Symbol, value::Float64)
    axis === :T && return SusceptibilityCLIConfig(
        T_MeV=value,
        muB_MeV=cfg.muB_MeV,
        muQ_MeV=cfg.muQ_MeV,
        muS_MeV=cfg.muS_MeV,
        V_fm3=cfg.V_fm3,
        scan_axis=cfg.scan_axis,
        scan_min_MeV=cfg.scan_min_MeV,
        scan_max_MeV=cfg.scan_max_MeV,
        scan_step_MeV=cfg.scan_step_MeV,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        output=cfg.output,
    )
    axis === :muB && return SusceptibilityCLIConfig(
        T_MeV=cfg.T_MeV,
        muB_MeV=value,
        muQ_MeV=cfg.muQ_MeV,
        muS_MeV=cfg.muS_MeV,
        V_fm3=cfg.V_fm3,
        scan_axis=cfg.scan_axis,
        scan_min_MeV=cfg.scan_min_MeV,
        scan_max_MeV=cfg.scan_max_MeV,
        scan_step_MeV=cfg.scan_step_MeV,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        output=cfg.output,
    )
    axis === :muQ && return SusceptibilityCLIConfig(
        T_MeV=cfg.T_MeV,
        muB_MeV=cfg.muB_MeV,
        muQ_MeV=value,
        muS_MeV=cfg.muS_MeV,
        V_fm3=cfg.V_fm3,
        scan_axis=cfg.scan_axis,
        scan_min_MeV=cfg.scan_min_MeV,
        scan_max_MeV=cfg.scan_max_MeV,
        scan_step_MeV=cfg.scan_step_MeV,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        output=cfg.output,
    )
    axis === :muS && return SusceptibilityCLIConfig(
        T_MeV=cfg.T_MeV,
        muB_MeV=cfg.muB_MeV,
        muQ_MeV=cfg.muQ_MeV,
        muS_MeV=value,
        V_fm3=cfg.V_fm3,
        scan_axis=cfg.scan_axis,
        scan_min_MeV=cfg.scan_min_MeV,
        scan_max_MeV=cfg.scan_max_MeV,
        scan_step_MeV=cfg.scan_step_MeV,
        xi=cfg.xi,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        output=cfg.output,
    )
    throw(ArgumentError("unsupported scan axis: $(axis)"))
end

function run_scan(cfg::SusceptibilityCLIConfig)
    cfg.scan_axis !== :none || throw(ArgumentError("run_scan requires scan_axis != :none"))
    cfg.scan_step_MeV > 0 || throw(ArgumentError("scan_step must be positive"))
    cfg.scan_max_MeV >= cfg.scan_min_MeV || throw(ArgumentError("scan_max must be >= scan_min"))

    values = collect(cfg.scan_min_MeV:cfg.scan_step_MeV:cfg.scan_max_MeV)
    isempty(values) && throw(ArgumentError("empty scan grid"))

    results = NamedTuple[]
    for value in values
        point_cfg = with_scan_value(cfg, cfg.scan_axis, value)
        push!(results, compute_summary(point_cfg))
    end
    return results
end

function main(args=ARGS)
    cfg = parse_args(args)
    output_path = cfg.output === nothing ? default_output_path(cfg) : cfg.output

    if cfg.scan_axis === :none
        result = compute_summary(cfg)
        print_summary(result)
        if cfg.output !== nothing
            maybe_write_csv(result, output_path)
            println()
            println("wrote CSV: $(output_path)")
        end
    else
        results = run_scan(cfg)
        println("=" ^ 68)
        println("PNJL conserved-charge susceptibility scan")
        println("=" ^ 68)
        println("scan axis : $(cfg.scan_axis)")
        println("scan grid : $(cfg.scan_min_MeV):$(cfg.scan_step_MeV):$(cfg.scan_max_MeV) MeV")
        println("points    : $(length(results))")
        println("fixed     : T=$(cfg.T_MeV) MeV, muB=$(cfg.muB_MeV) MeV, muQ=$(cfg.muQ_MeV) MeV, muS=$(cfg.muS_MeV) MeV, V=$(cfg.V_fm3) fm^3")
        write_scan_csv(results, output_path)
        println("wrote CSV: $(output_path)")
    end
end

main()