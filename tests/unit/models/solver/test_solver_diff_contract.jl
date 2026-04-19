using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver diff contract scaffold" begin
    T_fm = 100.0 / 197.327
    mu_fm = 0.0
    model = Models.create_model(:PNJL)
    mode = Models.FixedMu()
    result = Models.solve(model, mode, T_fm, mu_fm; p_num=8, t_num=4, residual_norm_max=1e-6)

    @testset "context and paramspec contract" begin
        ctx = Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )
        @test ctx isa Models.ThermoDiffContext

        params = Models.ParamSpec([:xi, :T_fm])
        @test params.names == [:T_fm, :xi]

        @test_throws ArgumentError Models.ParamSpec([:bad_param])
        @test_throws ArgumentError Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=NaN, mu_fm=mu_fm, xi=0.0),
        )

        ctx_mu_alias = Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, μ_fm=mu_fm, xi=0.0),
        )
        @test hasproperty(ctx_mu_alias.theta, :mu_fm)
        @test ctx_mu_alias.theta.mu_fm == mu_fm
        @test_throws ArgumentError Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, μ_fm=mu_fm, xi=0.0),
        )
    end

    @testset "target registry and unimplemented path" begin
        target = Models.diff_target(:pressure)
        @test target isa Models.DiffTarget
        @test_throws ArgumentError Models.diff_target(:unknown_target)

        ctx = Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )
        params = Models.ParamSpec([:T_fm])

        @test_throws ErrorException Models.jacobian(ctx, target, params)
    end

    @testset "jacobian shape convention" begin
        backend = function (ctx, target_or_targets, params)
            _ = ctx
            _ = params
            if target_or_targets isa AbstractVector
                return fill(2.5, length(target_or_targets))
            end
            return 2.5
        end
        ctx = Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
            jacobian_backend=backend,
        )
        target = Models.DiffTarget(:stub_target, c -> 0.0)
        params = Models.ParamSpec([:T_fm])

        J = Models.jacobian(ctx, target, params)
        @test size(J) == (1, 1)
        @test J[1, 1] == 2.5

        J2 = Models.jacobian(ctx, [target, target], params)
        @test size(J2) == (2, 1)

        params2 = Models.ParamSpec([:T_fm, :mu_fm])
        backend2 = function (ctx, target_or_targets, params)
            _ = ctx
            n_targets = target_or_targets isa AbstractVector ? length(target_or_targets) : 1
            return fill(1.25, n_targets, length(params.names))
        end
        ctx2 = Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
            jacobian_backend=backend2,
        )
        J3 = Models.jacobian(ctx2, target, params2)
        @test size(J3) == (1, 2)
        J4 = Models.jacobian(ctx2, [target, target], params2)
        @test size(J4) == (2, 2)

        ctx_default = Models.build_thermo_diff_context(
            result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )
        err_msg = try
            Models.jacobian(ctx_default, [target, target], params)
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("targets", err_msg)
    end
end
