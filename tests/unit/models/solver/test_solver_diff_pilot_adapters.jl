using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver diff pilot adapters contract" begin
    T_fm = 100.0 / 197.327
    mu_fm = 0.0
    model = Models.create_model(:PNJL)
    mode = Models.FixedMu()
    result = Models.solve(model, mode, T_fm, mu_fm; p_num=8, t_num=4, residual_norm_max=1e-6)

    @testset "mu_fm and μ_fm alias paths" begin
        ctx_mu = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )
        @test ctx_mu !== nothing

        ctx_mu_alias = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, μ_fm=mu_fm, xi=0.0),
        )
        @test ctx_mu_alias !== nothing

        out_mu = Models.eval_pilot_derivatives(
            ctx_mu;
            targets=[:pressure],
            params=[:T_fm],
        )
        out_mu_alias = Models.eval_pilot_derivatives(
            ctx_mu_alias;
            targets=[:pressure],
            params=[:T_fm],
        )

        @test size(out_mu.jacobian) == size(out_mu_alias.jacobian)
        @test Set(keys(out_mu.by_name)) == Set(keys(out_mu_alias.by_name))
        @test isapprox(out_mu.jacobian[1, 1], out_mu_alias.jacobian[1, 1])
    end

    @testset "eval_pilot_derivatives returns jacobian and by_name" begin
        ctx = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )

        out = Models.eval_pilot_derivatives(
            ctx;
            targets=[:pressure, :entropy],
            params=[:T_fm, :mu_fm],
        )

        @test hasproperty(out, :jacobian)
        @test hasproperty(out, :by_name)
        @test size(out.jacobian) == (2, 2)
        @test haskey(out.by_name, Symbol("pressure__dT_fm"))
        @test haskey(out.by_name, Symbol("entropy__dmu_fm"))
    end

    @testset "unknown target throws" begin
        ctx = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )

        @test_throws ArgumentError Models.eval_pilot_derivatives(
            ctx;
            targets=[:unknown_target],
            params=[:T_fm],
        )
    end

    @testset "invalid param throws" begin
        ctx = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )

        @test_throws ArgumentError Models.eval_pilot_derivatives(
            ctx;
            targets=[:pressure],
            params=[:bad_param],
        )
    end

    @testset "duplicate targets throw" begin
        ctx = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )

        err = try
            Models.eval_pilot_derivatives(
                ctx;
                target_names=[:pressure, :pressure],
                param_names=[:T_fm],
            )
            nothing
        catch ex
            ex
        end

        @test err isa ArgumentError
        @test occursin("target_names must not contain duplicates", sprint(showerror, err))
    end

    @testset "params alias equivalence on params keyword path" begin
        ctx = Models.build_pilot_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )

        out_ascii = Models.eval_pilot_derivatives(
            ctx;
            targets=[:pressure],
            params=[:mu_fm],
        )
        out_alias = Models.eval_pilot_derivatives(
            ctx;
            targets=[:pressure],
            params=[Symbol("μ_fm")],
        )

        @test size(out_ascii.jacobian) == size(out_alias.jacobian)
        @test Set(keys(out_ascii.by_name)) == Set(keys(out_alias.by_name))
        @test isapprox(out_ascii.jacobian[1, 1], out_alias.jacobian[1, 1])
        @test haskey(out_alias.by_name, Symbol("pressure__dmu_fm"))
    end
end
