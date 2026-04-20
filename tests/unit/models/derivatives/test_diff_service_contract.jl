using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "diff service contract" begin
    T_fm = 100.0 / 197.327
    mu_fm = 0.0
    model = Models.create_model(:PNJL)
    mode = Models.FixedMu()
    result = Models.solve(model, mode, T_fm, mu_fm; p_num=8, t_num=4, residual_norm_max=1e-6)

    @test isdefined(Models, :build_diff_service_context)
    @test isdefined(Models, :eval_diff_service_jacobian)

    backend = function (ctx, target_or_targets, params)
        _ = ctx
        _ = params
        if target_or_targets isa AbstractVector
            return fill(2.0, length(target_or_targets))
        end
        return 2.0
    end

    ctx = Models.build_diff_service_context(
        result;
        mode=mode,
        model=model,
        theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        jacobian_backend=backend,
    )

    payload = Models.eval_diff_service_jacobian(
        ctx;
        target_names=[:pressure],
        param_names=[:T_fm],
    )

    @test size(payload.jacobian) == (1, 1)
    @test payload.jacobian[1, 1] == 2.0
    @test payload.target_names == [:pressure]
    @test payload.param_names == [:T_fm]

    @testset "failure paths" begin
        @test_throws ArgumentError Models.eval_diff_service_jacobian(
            ctx;
            target_names=Symbol[],
            param_names=[:T_fm],
        )

        @test_throws ArgumentError Models.eval_diff_service_jacobian(
            ctx;
            target_names=[:pressure, :pressure],
            param_names=[:T_fm],
        )

        err = try
            Models.eval_diff_service_jacobian(
                ctx;
                target_names=[:pressure],
                param_names=[:μ_fm, :mu_fm],
            )
            nothing
        catch e
            e
        end

        @test err isa ArgumentError
        @test occursin("duplicated", sprint(showerror, err))
    end
end
