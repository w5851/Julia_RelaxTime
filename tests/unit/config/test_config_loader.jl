# ConfigLoader 模块单元测试

using Test
using TOML

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "config", "ConfigLoader.jl"))
if !isdefined(Main, :ConfigLoader)
    Base.include(Main, _CONFIG_LOADER_PATH)
end

using Main.ConfigLoader: deep_merge, load_config

@testset "ConfigLoader" begin
    @testset "deep_merge 浅层合并" begin
        base = Dict{String,Any}("a" => 1, "b" => 2)
        override = Dict{String,Any}("b" => 3, "c" => 4)
        merged = deep_merge(base, override)
        @test merged["a"] == 1
        @test merged["b"] == 3
        @test merged["c"] == 4
    end

    @testset "deep_merge 深层嵌套" begin
        base = Dict{String,Any}("x" => Dict{String,Any}("a" => 1, "b" => 2))
        override = Dict{String,Any}("x" => Dict{String,Any}("b" => 3))
        merged = deep_merge(base, override)
        @test merged["x"]["a"] == 1
        @test merged["x"]["b"] == 3
    end

    @testset "load_config 加载 models 目录默认配置" begin
        config_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "pnjl"))
        config, profile, path = load_config(config_dir, Dict{String,Any}())
        @test config isa Dict
        @test profile isa String
        @test path !== nothing
    end
end
