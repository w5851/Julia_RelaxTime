# ScanConfig 模块单元测试

using Test

const PROJECT_ROOT_SCF = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _SCAN_CONFIG_PATH = normpath(joinpath(PROJECT_ROOT_SCF, "src", "models", "scans", "ScanConfig.jl"))
if !isdefined(Main, :ScanConfig)
    Base.include(Main, _SCAN_CONFIG_PATH)
end

using Main.ScanConfig: TmuScanConfig, TrhoScanConfig, scan_kwargs

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
end
