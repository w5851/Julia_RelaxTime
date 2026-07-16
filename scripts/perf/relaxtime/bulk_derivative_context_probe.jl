using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MODELS_PATH = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")

if !isdefined(Main, :Models)
    Base.include(Main, MODELS_PATH)
end

using .Models

const THERMO = Models.ThermoDerivatives

function _arg_value(name::String, default::String)
    prefix = "--" * name * "="
    for arg in ARGS
        startswith(arg, prefix) && return arg[length(prefix) + 1:end]
    end
    return default
end

T_fm = parse(Float64, _arg_value("T_fm", "0.5"))
mu_fm = parse(Float64, _arg_value("mu_fm", "1.5"))
xi = parse(Float64, _arg_value("xi", "0.0"))
p_num = parse(Int, _arg_value("p_num", "12"))
t_num = parse(Int, _arg_value("t_num", "6"))
repeats = parse(Int, _arg_value("repeats", "5"))

model = Models.create_model(:PNJL)
equilibrium = Models.solve(
    model,
    Models.FixedMu(),
    T_fm,
    mu_fm;
    xi=xi,
    p_num=p_num,
    t_num=t_num,
)

function independent_series_path()
    pressure = THERMO._pressure_derivatives_order2_taylordiff(
        model,
        T_fm,
        mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=nothing,
        linear_solve=:auto,
        series_residual_tol=1e-7,
    )
    masses = THERMO._mass_derivatives_taylordiff(
        T_fm,
        mu_fm;
        order=1,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=nothing,
        linear_solve=:auto,
        series_residual_tol=1e-7,
    )
    return pressure, masses
end

function shared_context_path()
    return Models.bulk_viscosity_coefficients(
        T_fm,
        mu_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        base_state=equilibrium.x_state,
        base_masses=equilibrium.masses,
        base_mu_vec=equilibrium.mu_vec,
        base_state_source=:performance_probe_equilibrium,
    )
end

independent_series_path()
shared_context_path()

independent_times = Float64[]
independent_bytes = Int[]
shared_times = Float64[]
shared_bytes = Int[]
shared_result_ref = Ref{Any}(nothing)
for _ in 1:repeats
    GC.gc()
    independent = @timed independent_series_path()
    push!(independent_times, independent.time)
    push!(independent_bytes, independent.bytes)

    GC.gc()
    shared = @timed shared_context_path()
    push!(shared_times, shared.time)
    push!(shared_bytes, shared.bytes)
    shared_result_ref[] = shared.value
end

shared_result = shared_result_ref[]

println("PNJL bulk derivative context performance probe")
println("T_fm=$(T_fm)")
println("mu_fm=$(mu_fm)")
println("xi=$(xi)")
println("p_num=$(p_num)")
println("t_num=$(t_num)")
println("repeats=$(repeats)")
println("independent_series_primal_solve_count=5")
println("independent_series_jacobian_build_count=5")
println("shared_context_primal_solve_count=$(shared_result.primal_solve_count)")
println("shared_context_jacobian_factorization_count=$(shared_result.jacobian_factorization_count)")
println("shared_context_derivative_series_count=$(shared_result.derivative_series_count)")
println("independent_median_s=$(median(independent_times))")
println("shared_context_median_s=$(median(shared_times))")
println("median_speedup=$(median(independent_times) / median(shared_times))")
println("independent_median_alloc_bytes=$(median(independent_bytes))")
println("shared_context_median_alloc_bytes=$(median(shared_bytes))")
