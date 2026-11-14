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
        @test only(paths).diagnostics.anchor_strategy === :endpoint_anchor
        @test length(only(paths).anchors) == 1
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
            for value in (rows[mu], reference[mu], selected[index])
                c = Models.solver_convergence_view(value)
                @test c.solver_converged isa Bool
                @test c.residual_ok === true
                @test c.physical_valid === true
                @test c.source in (:path_nlsolve, :path_newton)
                @test c.residual_norm_max == 1e-6
            end
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

@testset "Automatic single anchor retains the three roots of a first-order slice" begin
    model = Models.create_model(:PNJL)
    hbarc = Main.Constants_PNJL.ħc_MeV_fm
    path = Models.FixedMuPath(120/hbarc,collect(296.0:0.5:318.0)./hbarc;p_num=24,t_num=8)
    result = Models.solve_path(model,path;continuation_strategy=Models.PALCContinuation())
    @test result.diagnostics.anchor_strategy === :endpoint_anchor
    @test result.diagnostics.anchor_solve_count == 1
    @test result.diagnostics.failure_count == 0
    target = 308/hbarc
    roots = [point for branch in result.branches for point in branch.points
        if point.param_value == target && point.converged]
    @test length(roots) == 3
    catalog = Models.solve_path(model,path;continuation_strategy=Models.PALCContinuation(),
        anchor_strategy=Models.SeedPoolAnchor())
    reference = sort([point.rho_norm for branch in catalog.branches for point in branch.points
        if point.param_value == target && point.converged])
    @test sort([point.rho_norm for point in roots]) ≈ reference atol=1e-6
    @test all(point -> point.residual_norm < 1e-6,roots)
end

@testset "Qualified density scans use natural continuation and preserve first-order structure" begin
    mktempdir() do directory
        grid = vcat(0.0, collect(0.2:0.15:3.5))
        rows = Dict{Float64,Any}()
        reference = Dict{Float64,Any}()
        paths = Models.PathSolveResult[]
        file = joinpath(directory, "auto.csv")
        common = (; T_values=[120.0], rho_values=grid, xi_values=[0.0],
            p_num=24, t_num=8, iterations=300, reverse_rho=true,
            constraint_mode=:fixed_rho, overwrite=true, resume=false)
        auto = Models.run_trho_scan(; common..., output_path=file,
            progress_cb=(point, result) -> (rows[point.rho] = result),
            path_diagnostics_cb=(point, result) -> push!(paths,result))
        single = Models.run_trho_scan(; common..., continuation=:single_pass,
            output_path=joinpath(directory, "single.csv"),
            progress_cb=(point, result) -> (reference[point.rho] = result))
        @test auto.success == single.success == length(grid)
        @test length(paths) == 1
        @test only(paths).diagnostics.continuation_backend === :natural_parameter
        @test only(paths).diagnostics.anchor_solve_count == 1
        @test all(line -> occursin("continuation=natural", line), readlines(file)[2:end])
        for rho in grid
            @test rows[rho].residual_norm < 1e-6
            @test abs(rows[rho].rho_norm - rho) < 1e-6
            @test maximum(abs.(rows[rho].solution .- reference[rho].solution)) < 1e-6
            @test abs(rows[rho].pressure - reference[rho].pressure) < 1e-6
            @test Models.solver_convergence_view(rows[rho]).residual_ok === true
            @test Models.solver_convergence_view(rows[rho]).source in (:path_nlsolve, :path_newton)
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

        phase = Models.run_phase_pipeline(:PNJL;mode=:production,T_grid=[120.0],
            rho_grid=grid,xi=0.0,output_dir=joinpath(directory,"phase"),continuation=:auto,
            seed_policy=:candidates,memoize_uniform=true,p_num=24,t_num=8,iterations=300,
            cep_max_refine_level=0,rho_geometry_convergence=false,promote_reference=false)
        @test phase.config_snapshot["continuation"] == "auto"
        phase_csv = joinpath(directory,"phase","trho_scan.csv")
        @test all(line -> occursin("continuation=natural",line),readlines(phase_csv)[2:end])
        phase_curves = Models.load_curves_from_trho_csv(phase_csv;xi=0.0,min_points=3)
        phase_mu,phase_rho = phase_curves[120.0]
        phase_boundary = Models.maxwell_construction(phase_mu,phase_rho;
            spinodal_hint=Models.detect_s_shape(phase_mu,phase_rho))
        @test phase_boundary.converged
        @test phase_boundary.mu_transition ≈ boundary.mu_transition atol=1e-6
    end
end
