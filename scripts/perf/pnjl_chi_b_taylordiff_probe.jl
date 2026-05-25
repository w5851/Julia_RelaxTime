using Statistics

const SCRIPT_START = time()
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS_PATH = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
const CONSTANTS_PNJL_PATH = joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl")
const CEP_REFERENCE_PATH = joinpath(PROJECT_ROOT, "data", "reference", "pnjl", "cep.csv")

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, CONSTANTS_PNJL_PATH)
end
if !isdefined(Main, :Models)
    Base.include(Main, MODELS_PATH)
end

using .Models
using Main.Constants_PNJL: ħc_MeV_fm

const TD = Models.PNJLChiBTaylorDiff

function _arg_value(name::String, default::String)
    prefix = "--" * name * "="
    for arg in ARGS
        startswith(arg, prefix) && return arg[length(prefix) + 1:end]
    end
    return default
end

function _load_xi0_cep(path::String)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    idx_xi = findfirst(==("xi"), header)
    idx_T = findfirst(==("T_CEP_MeV"), header)
    idx_muB = findfirst(==("muB_CEP_MeV"), header)
    for line in lines[2:end]
        cols = split(strip(line), ',')
        parse(Float64, cols[idx_xi]) == 0.0 || continue
        return (
            T_CEP_MeV=parse(Float64, cols[idx_T]),
            muB_CEP_MeV=parse(Float64, cols[idx_muB]),
        )
    end
    error("xi=0 CEP reference row not found: $path")
end

function _rss_bytes()
    try
        return Sys.maxrss()
    catch
        return -1
    end
end

max_order = parse(Int, _arg_value("max_order", "10"))
p_num = parse(Int, _arg_value("p_num", "8"))
t_num = parse(Int, _arg_value("t_num", "4"))
repeats = parse(Int, _arg_value("repeats", "7"))
linear_solve = Symbol(_arg_value("linear_solve", "auto"))
series_residual_tol = parse(Float64, _arg_value("series_residual_tol", "1e-6"))

cep = _load_xi0_cep(CEP_REFERENCE_PATH)
T_fm = cep.T_CEP_MeV / ħc_MeV_fm
muB_fm = cep.muB_CEP_MeV / ħc_MeV_fm
T_next_fm = (cep.T_CEP_MeV + 0.5) / ħc_MeV_fm

kwargs = (;
    max_order=max_order,
    xi=0.0,
    p_num=p_num,
    t_num=t_num,
    linear_solve=linear_solve,
    series_residual_tol=series_residual_tol,
)

GC.gc()
first = @timed TD.chi_B_taylordiff_all(T_fm, muB_fm; kwargs...)
GC.gc()
next_point = @timed TD.chi_B_taylordiff_all(T_next_fm, muB_fm; kwargs...)

warm_times = Float64[]
warm_allocs = Int[]
for _ in 1:repeats
    GC.gc()
    t = @timed TD.chi_B_taylordiff_all(T_fm, muB_fm; kwargs...)
    push!(warm_times, t.time)
    push!(warm_allocs, t.bytes)
end

println("PNJL chi_B TaylorDiff performance probe")
println("T_CEP_MeV=$(cep.T_CEP_MeV)")
println("muB_CEP_MeV=$(cep.muB_CEP_MeV)")
println("max_order=$(max_order)")
println("p_num=$(p_num)")
println("t_num=$(t_num)")
println("linear_solve=$(linear_solve)")
println("first_call_s=$(first.time)")
println("first_alloc_bytes=$(first.bytes)")
println("same_process_next_point_s=$(next_point.time)")
println("same_process_next_point_alloc_bytes=$(next_point.bytes)")
println("warm_median_s=$(median(warm_times))")
println("warm_min_s=$(minimum(warm_times))")
println("warm_alloc_median_bytes=$(median(warm_allocs))")
println("peak_rss_bytes=$(_rss_bytes())")
println("cold_wall_s=$(time() - SCRIPT_START)")
println("values=$(first.value)")
