using Test
if !isdefined(Main,:Models)
    include(joinpath(@__DIR__,"..","..","..","src","models","Models.jl"))
end

@testset "FixedMu and solve_gap enforce the requested iteration budget" begin
    model=Main.Models.create_model(:PNJL)
    seed=collect(Main.Models.HADRON_SEED_5)
    options=(seed_guess=seed,seed_candidates=[seed],p_num=12,t_num=6,residual_norm_max=1e-8)
    @test_throws ArgumentError Main.Models.solve_constraint(model,Main.Models.FixedMu(),1.1;
        μ_fm=0.0,options...,iterations=0)
    @test_throws ArgumentError Main.Models.solve_gap(model,1.1,0.0;
        initial_guess=seed,p_num=12,t_num=6,iterations=0)
    @test_throws ErrorException Main.Models.solve_gap(model,1.1,0.0;
        initial_guess=seed,p_num=12,t_num=6,iterations=1)
    limited=Main.Models.solve_constraint(model,Main.Models.FixedMu(),1.1;
        μ_fm=0.0,options...,iterations=1)
    @test !limited.converged
    converged=Main.Models.solve_constraint(model,Main.Models.FixedMu(),1.1;
        μ_fm=0.0,options...,iterations=40)
    @test converged.converged
    @test converged.residual_norm<=1e-8
end

@testset "FixedMu forwards an explicit local seed pool" begin
    model=Main.Models.create_model(:PNJL)
    T,mu=1.1,0.0
    original=Main.Models.solve(model,Main.Models.FixedMu(),T,mu;p_num=12,t_num=6)
    @test original.converged
    seed=collect(original.x_state)
    local_result=Main.Models.solve_constraint(model,Main.Models.FixedMu(),T;μ_fm=mu,
        seed_guess=seed,seed_candidates=[seed],p_num=12,t_num=6,residual_norm_max=1e-8)
    @test local_result.converged
    @test local_result.candidate_count==1
    @test local_result.residual_norm<=1e-8
    @test all(isapprox.(local_result.masses,original.masses;rtol=1e-7,atol=1e-10))
    perturbed=Main.Models.solve_constraint(model,Main.Models.FixedMu(),T;μ_fm=mu,
        seed_guess=seed,p_num=12,t_num=6,residual_norm_max=1e-8)
    @test perturbed.converged
    @test perturbed.candidate_count>local_result.candidate_count
    @test all(isapprox.(local_result.masses,perturbed.masses;rtol=1e-7,atol=1e-10))
end
