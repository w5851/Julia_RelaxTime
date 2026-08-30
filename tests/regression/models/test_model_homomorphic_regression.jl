using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "model homomorphic regression" begin
    kinds = (:NJL, :NJL2, :PNJL, :PNJLMagnetic, :RPNJL, :Rotation, :GasLiquid)

    for kind in kinds
        model = kind === :PNJLMagnetic ?
            Models.create_model(kind; eB_fm2=0.1) : Models.create_model(kind)
        dim = Models.gap_state_dim(model)
        @test dim >= 2

        schema = Models.schema_for_model(kind)
        @test Models.state_dim(schema) == dim

        if kind === :Rotation || kind === :GasLiquid
            # 这两个模型工作流当前是专用入口，回归阶段先验证入口存在。
            if kind === :Rotation
                @test isdefined(Models, :solve_rotation_point)
            else
                @test isdefined(Models, :solve_gas_liquid_point)
            end
            continue
        end

        if kind === :PNJLMagnetic
            # Magnetic equilibrium is a dedicated FixedMu route; the shared
            # scan/constraint regression must not launch a full magnetic solve.
            @test isdefined(Models, :run_magnetic_scan)
            @test_throws ArgumentError Models.run_tmu_scan(
                T_values=[150.0], mu_values=[0.0], xi_values=[0.0],
                output_path=joinpath(mktempdir(), "magnetic_tmu.csv"),
                model_kind=:PNJLMagnetic,
            )
            @test_throws ArgumentError Models.run_trho_scan(
                T_values=[150.0], rho_values=[0.2], xi_values=[0.0],
                output_path=joinpath(mktempdir(), "magnetic_trho.csv"),
                model_kind=:PNJLMagnetic,
            )
            continue
        end

        mu_vec = if kind === :NJL2
            [0.0, 0.0, 0.0]
        else
            [0.0, 0.0, 0.0]
        end
        st = Models.solve_gap(model, 0.5, mu_vec; p_num=8, t_num=4)
        @test st isa Models.MeanFieldState
        @test length(Models.state_vector(st)) == 5
    end
end
