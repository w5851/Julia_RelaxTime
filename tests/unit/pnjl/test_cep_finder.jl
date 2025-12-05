using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.PhaseTransition
using .PNJL.CEPFinder

function synthetic_curves()
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    grouped[60.0] = [(μ, 0.02 * μ + 0.4 * sin(μ / 40)) for μ in -80.0:20.0:80.0]
    grouped[80.0] = [(μ, 0.01 * μ + 0.6 * sin(μ / 30)) for μ in -80.0:20.0:80.0]
    grouped[100.0] = [(μ, 0.03 * μ) for μ in -80.0:20.0:80.0]
    return grouped
end

function s_shape_curve()
    mu = Float64[-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    rho = Float64[0.0, 0.6, 1.0, 0.5, 0.2, 0.8]
    return (mu, rho)
end

function monotonic_curve()
    mu = Float64[-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    rho = Float64[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    return (mu, rho)
end

@testset "build_curves" begin
    grouped = Dict(70.0 => [(0.0, 0.0), (10.0, 0.1), (20.0, 0.2)])
    curves = build_curves(grouped)
    @test haskey(curves, 70.0)
    mu, rho = curves[70.0]
    @test length(mu) == 3 == length(rho)
end

@testset "find_cep auto fetch" begin
    curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}(
        80.0 => s_shape_curve(),
        120.0 => monotonic_curve(),
    )
    invoked = Float64[]
    header = "T_MeV,rho,xi,mu_u_MeV,mu_d_MeV,mu_s_MeV,mu_avg_MeV,pressure_fm4,entropy_fm3,energy_fm4,phi_u,phi_d,phi_s,Phi1,Phi2,iterations,residual_norm,converged,message"
    fake_runner = function (; T_values, rho_values, xi_values, output_path, seed_path, overwrite, resume, p_num, t_num)
        push!(invoked, T_values[1])
        mkpath(dirname(output_path))
        open(output_path, "w") do io
            println(io, header)
            for (idx, rho) in enumerate(rho_values)
                mu = 250.0 + rho + idx
                values = (
                    T_values[1],
                    rho,
                    xi_values[1],
                    mu,
                    mu,
                    mu,
                    mu,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    5,
                    1e-6,
                    "true",
                    "",
                )
                println(io, join(values, ','))
            end
        end
    end
    raw_path = tempname() * "_raw.csv"
    processed_path = tempname() * "_curves.csv"
    result = find_cep(curves; temp_tol=1.0, auto_fetch=(
        xi = 0.0,
        runner = fake_runner,
        rho_values = [0.0, 0.5],
        seed_path = "",
        output_path = raw_path,
        processed_path = processed_path,
        adaptive = false,
    ))
    @test result.has_cep
    @test any(t -> isapprox(t, 100.0; atol=1e-8), invoked)
    @test haskey(curves, invoked[end])
    @test isfile(raw_path)
    @test isfile(processed_path)
    raw_lines = readlines(raw_path)
    proc_lines = readlines(processed_path)
    @test startswith(raw_lines[1], "T_MeV")
    @test length(proc_lines) > 1
    rm(raw_path; force=true)
    rm(processed_path; force=true)
end

@testset "find_cep" begin
    grouped = synthetic_curves()
    curves = build_curves(grouped)
    result = find_cep(curves; auto_fetch=nothing)
    @test result.has_cep
    @test result.T_cep_MeV !== nothing
    @test result.bracket == (80.0, 100.0)
end
