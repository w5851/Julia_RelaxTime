using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "Qualified constraint paths default to PALC with complete joint residuals" begin
    model = Models.create_model(:PNJL)
    hbarc = Main.Constants_PNJL.ħc_MeV_fm
    for (mode, temperatures) in (
        (Models.FixedSigma(12.0759), collect(100.0:2.0:140.0)),
        (Models.FixedMuBConservedCharges(900/hbarc), collect(100.0:3.0:160.0)),
    )
        path = Models.ConstraintTPath(mode, reverse(temperatures)./hbarc; p_num=24, t_num=8)
        result = Models.solve_path(model, path; iterations=300)
        @test result.diagnostics.continuation_backend === :bifurcationkit
        @test result.diagnostics.failure_count == 0
        @test result.diagnostics.anchor_strategy === :endpoint_anchor
        @test result.diagnostics.anchor_solve_count == 1
        @test length(result.anchors) == 1
        values = Models.selected_path_results(model, result)
        options = Models._path_options(model, path)
        residual = Models._path_residual(model, path, options)
        for (temperature, value) in zip(path.T_values, values)
            @test value !== nothing
            value === nothing && continue
            @test value.converged
            @test length(value.solution) == 8
            @test maximum(abs, residual(value.solution, temperature)) < 1e-6
            @test value.masses[3] > value.masses[1] > 0
            c = Models.solver_convergence_view(value)
            @test c.solver_converged isa Bool
            @test c.residual_ok === true
            @test c.physical_valid === true
            @test c.source in (:path_nlsolve, :path_newton)
        end
        for i in (1,11,21)
            legacy = Models.solve(model, mode, path.T_values[i];
                seed_guess=values[i].solution, p_num=24, t_num=8, iterations=300)
            @test legacy.converged
            @test maximum(abs.(legacy.solution .- values[i].solution)) < 1e-4
            @test abs(legacy.pressure - values[i].pressure) < 1e-6
        end
    end
end

@testset "Isentropic meson entry consumes the PALC equilibrium in input order" begin
    mktempdir() do directory
        grid = collect(140.0:-2.0:100.0)
        paths = Any[]
        result = Models.run_isentropic_meson_mass_scan(; T_MeV_values=grid,
            sigma_target=12.0759, traversal=:T_descending, mesons=(), p_num=24, t_num=8,
            output_path=joinpath(directory, "isentropic.csv"),
            path_diagnostics_cb=(indices, path) -> push!(paths, (indices, path)))
        @test result.palc_paths == 1
        @test only(paths)[1] == 1:21
        selected = Models.selected_path_results(Models.create_model(:PNJL), only(paths)[2])
        lines = readlines(result.output_path)
        header = split(first(lines), ',')
        columns = Dict(name=>i for (i,name) in enumerate(header))
        rows = [split(line, ','; keepempty=true) for line in lines[2:end]]
        @test length(rows) == length(grid)
        @test [parse(Float64,row[columns["T_MeV"]]) for row in rows] == grid
        for (row, value) in zip(rows, selected)
            @test strip(row[columns["message"]], '"') == ""
            @test parse(Bool,row[columns["equilibrium_converged"]])
            @test parse(Float64,row[columns["equilibrium_residual_norm"]]) < 1e-6
            @test parse(Float64,row[columns["Phi"]]) ≈ value.x_state[4] atol=1e-6
            @test parse(Float64,row[columns["muB_MeV"]]) ≈
                3*value.mu_vec[1]*Main.Constants_PNJL.ħc_MeV_fm atol=1e-5
        end
    end
end
