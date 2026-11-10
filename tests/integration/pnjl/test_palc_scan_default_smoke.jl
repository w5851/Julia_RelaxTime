using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "Qualified chemical-potential scans use PALC and preserve the density jump" begin
    mktempdir() do directory
        grid = collect(0.0:15.0:450.0)
        rows = Dict{Float64,Any}()
        reference = Dict{Float64,Any}()
        paths = Models.PathSolveResult[]
        callback_order = Float64[]
        file = joinpath(directory, "auto_tmu.csv")
        config = Models.TmuScanConfig(; T_values=[120.0], mu_values=grid, xi_values=[0.0],
            p_num=24, t_num=8, overwrite=true, resume=false, nlsolve_kwargs=(; iterations=300))
        auto = Models.run_tmu_scan(config; output_path=file,
            progress_cb=(point, result) -> begin
                push!(callback_order, point.mu)
                rows[point.mu] = result
            end,
            path_diagnostics_cb=(point, result) -> push!(paths, result))
        single = Models.run_tmu_scan(config; continuation=:single_pass,
            output_path=joinpath(directory, "single_tmu.csv"),
            progress_cb=(point, result) -> (reference[point.mu] = result))
        @test auto.success == single.success == length(grid)
        @test auto.failure == single.failure == 0
        @test length(paths) == 1
        @test only(paths).diagnostics.continuation_backend === :bifurcationkit
        @test callback_order == grid
        lines = readlines(file)
        @test all(line -> occursin("continuation=palc", line), lines[2:end])
        @test [parse(Float64, split(line, ',')[2]) for line in lines[2:end]] == grid

        model = Models.create_model(:PNJL)
        selected = Models.selected_path_results(model, only(paths); iterations=300)
        hbarc = Main.Constants_PNJL.ħc_MeV_fm
        rho0 = Main.Constants_PNJL.ρ0_inv_fm3
        for (index, mu) in enumerate(grid)
            @test rows[mu].residual_norm < 1e-6
            @test maximum(abs.(rows[mu].mu_vec .- mu / hbarc)) < 1e-10
            @test maximum(abs.(rows[mu].solution .- reference[mu].solution)) < 1e-6
            @test abs(rows[mu].pressure - reference[mu].pressure) < 1e-6
            @test abs(rows[mu].rho_norm - reference[mu].rho_norm) < 1e-7
            @test rows[mu].rho_norm ≈ selected[index].rho_norm * rho0
        end
        @test abs(rows[0.0].rho_norm) < 1e-8
        # At T=120 MeV the pressure-selected branch changes in this mu interval.
        @test rows[315.0].rho_norm - rows[300.0].rho_norm > rho0

        partial = joinpath(directory, "resume_tmu.csv")
        write(partial, join(lines[1:4], '\n') * "\n")
        resumed = Models.run_tmu_scan(config; output_path=partial, overwrite=false, resume=true)
        @test resumed.skipped == 3
        @test resumed.success == length(grid)-3
        @test resumed.failure == 0
        @test read(partial) == read(file)
    end
end

@testset "Qualified density scans use PALC and preserve first-order structure" begin
    mktempdir() do directory
        grid = vcat(0.0, collect(0.2:0.15:3.5))
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
        mu_q_MeV = [sum(rows[rho].mu_vec)/3 * Main.Constants_PNJL.ħc_MeV_fm for rho in grid]
        @test any(<(0), diff(mu_q_MeV))
        shape = Models.detect_s_shape(mu_q_MeV, grid)
        @test shape.has_s_shape
        boundary = Models.maxwell_construction(mu_q_MeV, grid; spinodal_hint=shape)
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
