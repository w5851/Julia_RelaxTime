# ScanConfig 模块单元测试

using Test

const PROJECT_ROOT_SCF = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _SCAN_CONFIG_PATH = normpath(joinpath(PROJECT_ROOT_SCF, "src", "models", "scans", "ScanConfig.jl"))
if !isdefined(Main, :ScanConfig)
    Base.include(Main, _SCAN_CONFIG_PATH)
end

using Main.ScanConfig: TmuScanConfig, TrhoScanConfig, FreezeoutScanConfig, scan_kwargs

@testset "ScanConfig" begin
    @testset "TmuScanConfig 默认构建" begin
        cfg = TmuScanConfig(
            T_values=[0.5, 0.6, 0.7],
            mu_values=[0.0, 0.1],
            output_path=tempname()
        )
        @test cfg.T_values == [0.5, 0.6, 0.7]
        @test cfg.mu_values == [0.0, 0.1]
    end

    @testset "TrhoScanConfig 默认构建" begin
        cfg = TrhoScanConfig(
            T_values=[0.5],
            rho_values=[0.0, 0.5],
            output_path=tempname()
        )
        @test cfg.T_values == [0.5]
    end

    @testset "scan_kwargs 过滤 nothing 字段" begin
        cfg = TmuScanConfig(
            T_values=[0.5],
            mu_values=[0.0],
            output_path=tempname()
        )
        kw = scan_kwargs(cfg)
        @test kw isa NamedTuple
    end

    @testset "scan config 保留显式热积分策略" begin
        cfg = TrhoScanConfig(
            T_values=[150.0],
            rho_values=[0.1],
            output_path=tempname(),
            thermo_quadrature_policy=:rs_reduced_adaptive,
            thermo_quadrature_rtol=1e-7,
            thermo_quadrature_atol=1e-9,
            thermo_quadrature_maxevals=12345,
        )
        kw = scan_kwargs(cfg)
        @test kw.thermo_quadrature_policy === :rs_reduced_adaptive
        @test kw.thermo_quadrature_rtol == 1e-7
        @test kw.thermo_quadrature_atol == 1e-9
        @test kw.thermo_quadrature_maxevals == 12345
    end

    @testset "FreezeoutScanConfig 默认构建" begin
        cfg = FreezeoutScanConfig(
            sqrt_s_NN_values=[7.7, 11.5],
            xi_values=[0.0],
            output_path=tempname(),
            profile_name="default",
            path_profile_name="baseline_freezeout",
        )
        @test cfg.sqrt_s_NN_values == [7.7, 11.5]
        @test cfg.profile_name == "default"
        @test cfg.path_profile_name == "baseline_freezeout"
        kw = scan_kwargs(cfg)
        @test haskey(kw, :sqrt_s_NN_values)
        @test haskey(kw, :profile_name)
        @test haskey(kw, :path_profile_name)
    end
end
