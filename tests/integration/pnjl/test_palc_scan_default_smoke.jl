using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "Qualified density scans use PALC and preserve first-order structure" begin
    mktempdir() do directory
        grid = collect(0.2:0.15:3.5)
        rows = Dict{Float64,Any}()
        reference = Dict{Float64,Any}()
        paths = Ref(0)
        file = joinpath(directory, "auto.csv")
        common = (; T_values=[120.0], rho_values=grid, xi_values=[0.0],
            p_num=24, t_num=8, iterations=300, reverse_rho=true,
            constraint_mode=:fixed_rho, overwrite=true, resume=false)
        auto = Models.run_trho_scan(; common..., output_path=file,
            progress_cb=(point, result) -> (rows[point.rho] = result),
            path_diagnostics_cb=(point, result) -> (paths[] += 1))
        single = Models.run_trho_scan(; common..., continuation=:single_pass,
            output_path=joinpath(directory, "single.csv"),
            progress_cb=(point, result) -> (reference[point.rho] = result))
        @test auto.success == single.success == length(grid)
        @test paths[] == 1
        @test all(line -> occursin("continuation=palc", line), readlines(file)[2:end])
        for rho in grid
            @test rows[rho].residual_norm < 1e-6
            @test abs(rows[rho].rho_norm - rho) < 1e-6
            @test maximum(abs.(rows[rho].solution .- reference[rho].solution)) < 1e-6
            @test abs(rows[rho].pressure - reference[rho].pressure) < 1e-6
        end
        mu = [sum(rows[rho].mu_vec)/3 for rho in grid]
        @test any(<(0), diff(mu))
        shape = Models.detect_s_shape(mu, grid)
        @test shape.has_s_shape
        boundary = Models.maxwell_construction(mu, grid; spinodal_hint=shape)
        @test boundary.converged
        @test boundary.rho_hadron < boundary.rho_quark

        # A partial file retains its existing row and is filled in scan order.
        partial = joinpath(directory, "resume.csv")
        write(partial, join(readlines(file)[1:2], '\n') * "\n")
        resumed = Models.run_trho_scan(; common..., output_path=partial, overwrite=false, resume=true)
        @test resumed.skipped == 1
        @test resumed.success == length(grid)-1
        @test read(partial) == read(file)
    end
end
