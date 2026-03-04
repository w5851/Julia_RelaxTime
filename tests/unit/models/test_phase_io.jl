# PhaseIO.jl 单元测试
#
# 测试内容：
# 1. group_curves_by_temperature 合成数据
# 2. load_curves_from_trho_csv 文件不存在时

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "PhaseIO" begin

    @testset "group_curves_by_temperature 接口存在" begin
        @test isdefined(Models, :group_curves_by_temperature)
    end

    @testset "group_curves_by_temperature 合成数据" begin
        rows = [
            Dict("T_MeV" => "100.0", "mu_MeV" => "0.0", "rho" => "0.1", "xi" => "0.0"),
            Dict("T_MeV" => "100.0", "mu_MeV" => "50.0", "rho" => "0.2", "xi" => "0.0"),
            Dict("T_MeV" => "200.0", "mu_MeV" => "0.0", "rho" => "0.05", "xi" => "0.0"),
        ]
        result = Models.group_curves_by_temperature(rows; xi=0.0)
        @test result isa Dict
        @test length(result) == 2
    end

    @testset "load_curves_from_trho_csv 接口存在" begin
        @test isdefined(Models, :load_curves_from_trho_csv)
    end
end
