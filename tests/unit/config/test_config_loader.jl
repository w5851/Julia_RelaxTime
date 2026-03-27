# ConfigLoader 模块单元测试

using Test
using TOML

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "config", "ConfigLoader.jl"))
if !isdefined(Main, :ConfigLoader)
    Base.include(Main, _CONFIG_LOADER_PATH)
end

using Main.ConfigLoader: config_loader_cache_stats, deep_merge, load_config, reset_config_loader_cache!

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

    @testset "test_config_loader_cache_hit_miss_reset" begin
        reset_config_loader_cache!()
        dir = mktempdir()
        try
            write(joinpath(dir, "default.toml"), "a = 1\n")
            default_config = Dict{String, Any}("root" => Dict{String, Any}("x" => 1))

            _ = load_config(dir, default_config; profile="default")
            stats1 = config_loader_cache_stats()
            @test stats1.miss == 1
            @test stats1.hit == 0

            _ = load_config(dir, default_config; profile="default")
            stats2 = config_loader_cache_stats()
            @test stats2.miss == 1
            @test stats2.hit == 1

            reset_config_loader_cache!()
            stats3 = config_loader_cache_stats()
            @test stats3.miss == 0
            @test stats3.hit == 0
            @test stats3.entries == 0
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "test_config_loader_cache_invalidate_on_profile_content_change" begin
        reset_config_loader_cache!()
        dir = mktempdir()
        try
            cfg_path = joinpath(dir, "default.toml")
            write(cfg_path, "a = 1\n")
            default_config = Dict{String, Any}()

            first = load_config(dir, default_config; profile="default")
            @test first.config["a"] == 1

            write(cfg_path, "a = 2\n")
            second = load_config(dir, default_config; profile="default")
            @test second.config["a"] == 2

            stats = config_loader_cache_stats()
            @test stats.miss == 2
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "test_config_loader_cache_returns_defensive_copy" begin
        reset_config_loader_cache!()
        dir = mktempdir()
        try
            write(joinpath(dir, "default.toml"), "[section]\nvalue = 3\n")
            default_config = Dict{String, Any}("section" => Dict{String, Any}("other" => 1))

            result1 = load_config(dir, default_config; profile="default")
            result1.config["section"]["value"] = 999

            result2 = load_config(dir, default_config; profile="default")
            @test result2.config["section"]["value"] == 3
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "test_config_loader_fingerprint_key_sort_and_normalized_dir" begin
        reset_config_loader_cache!()
        dir = mktempdir()
        try
            write(joinpath(dir, "default.toml"), "x = 1\n")
            dir_variant = joinpath(dir, ".", "sub", "..")
            config_a = Dict{String, Any}("b" => 2, "a" => 1)
            config_b = Dict{String, Any}("a" => 1, "b" => 2)

            _ = load_config(dir_variant, config_a; profile="default")
            _ = load_config(dir, config_b; profile="default")

            stats = config_loader_cache_stats()
            @test stats.miss == 1
            @test stats.hit == 1
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "test_config_loader_fingerprint_float_negzero_and_naninf_policy" begin
        reset_config_loader_cache!()
        dir = mktempdir()
        try
            write(joinpath(dir, "default.toml"), "x = 1\n")

            _ = load_config(dir, Dict{String, Any}("x" => -0.0); profile="default")
            _ = load_config(dir, Dict{String, Any}("x" => 0.0); profile="default")
            stats = config_loader_cache_stats()
            @test stats.miss == 1
            @test stats.hit == 1

            @test_throws ArgumentError load_config(dir, Dict{String, Any}("x" => NaN); profile="default")
            @test_throws ArgumentError load_config(dir, Dict{String, Any}("x" => Inf); profile="default")
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "test_config_loader_cache_concurrent_reads_no_crash_no_dirty_write" begin
        reset_config_loader_cache!()
        dir = mktempdir()
        try
            write(joinpath(dir, "default.toml"), "[solver]\nmax_iter = 128\n")
            default_config = Dict{String, Any}("solver" => Dict{String, Any}("tol" => 1.0e-6))

            tasks = [
                Threads.@spawn begin
                    res = load_config(dir, default_config; profile="default")
                    res.config["solver"]["max_iter"]
                end for _ in 1:16
            ]
            values = fetch.(tasks)
            @test all(==(128), values)

            probe = load_config(dir, default_config; profile="default")
            @test probe.config["solver"]["max_iter"] == 128
        finally
            rm(dir; recursive=true, force=true)
        end
    end
end
