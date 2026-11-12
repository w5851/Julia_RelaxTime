using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

function assert_same_path_roots(a, b)
    @test a.diagnostics.failure_count == b.diagnostics.failure_count == 0
    for value in Models._path_values(a.path)
        left = [p for branch in a.branches for p in branch.points if p.param_value == value && p.converged]
        right = [p for branch in b.branches for p in branch.points if p.param_value == value && p.converged]
        @test length(left) == length(right)
        for p in left
            matches = filter(q -> maximum(abs.(p.solution .- q.solution)) <= 1e-4, right)
            @test length(matches) == 1
            isempty(matches) && continue
            @test abs(p.pressure-first(matches).pressure) < 1e-6
            @test first(matches).residual_norm < 1e-6
        end
    end
end

@testset "Temperature reuse preserves every Tmu root and scan ordering" begin
    mktempdir() do directory
        common = (; T_values=[130.0,120.0,110.0], mu_values=collect(296.0:0.5:320.0),
            xi_values=[0.0,0.4], p_num=24, t_num=8, iterations=300,
            overwrite=true, resume=false, continuation=:palc)
        reference = Models.PathSolveResult[]
        warm = Models.PathSolveResult[]
        order = Tuple{Float64,Float64,Float64}[]
        Models.run_tmu_scan(;common...,reuse_temperature_anchor=false,
            output_path=joinpath(directory,"reference.csv"),
            path_diagnostics_cb=(_,result)->push!(reference,result))
        result = Models.run_tmu_scan(;common...,reuse_temperature_anchor=true,
            palc_options=(;jacobian=:forwarddiff,target_predictor=:cubic),
            reuse_corrector=true,corrector_method=:newton,
            output_path=joinpath(directory,"warm.csv"),
            progress_cb=(point,_)->push!(order,(point.T,point.mu,point.xi)),
            path_diagnostics_cb=(_,result)->push!(warm,result))
        @test result.failure == 0
        @test result.success == length(order)
        @test order == [(T,mu,xi) for xi in common.xi_values for T in common.T_values for mu in common.mu_values]
        for (index,(a,b)) in enumerate(zip(reference,warm))
            assert_same_path_roots(a,b)
            snapshot = b.branches[1].diagnostics.work
            @test snapshot.reused_anchor_solves == (index in (1,4) ? 0 : 1)
        end
    end
end

@testset "Density session reuses an endpoint and recalculates the new temperature" begin
    hbarc = Main.Constants_PNJL.ħc_MeV_fm
    model = Models.create_model(:PNJL)
    grid = vcat(0.0, collect(0.2:0.15:3.5))
    session = Models.TrhoScan.new_rho_point_session(;continuation=:natural,
        reuse_temperature_anchor=true,iterations=300)
    for (index,T) in enumerate((130.0,120.0,110.0))
        result = Models.TrhoScan.prepare_rho_path!(session,T,0.0,grid)
        reference = Models.solve_path(model,Models.FixedRhoPath(T/hbarc,grid;p_num=24,t_num=8);
            continuation_strategy=Models.NaturalContinuation(),iterations=300)
        assert_same_path_roots(reference,result)
        @test result.branches[1].diagnostics.work.reused_anchor_solves == (index == 1 ? 0 : 1)
        for rho in grid
            value = session.prepared[(T,0.0,rho)]
            @test value !== nothing && value.converged
            @test abs(value.rho_norm-rho) < 1e-6
        end
    end
end
