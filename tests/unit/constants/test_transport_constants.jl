# TransportConstants 单元测试
#
# 测试内容：
# 1. SCATTERING_MESON_MAP 结构完整性
# 2. SCATTERING_PROCESS_KEYS 枚举
# 3. 散射类型/通道结构正确性

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
using .Constants_PNJL: SCATTERING_MESON_MAP, SCATTERING_PROCESS_KEYS

# ============================================================================
# SCATTERING_MESON_MAP 结构测试
# ============================================================================

@testset "TransportConstants" begin

    @testset "SCATTERING_MESON_MAP 非空" begin
        @test !isempty(SCATTERING_MESON_MAP)
        @test length(SCATTERING_MESON_MAP) >= 17  # 4 qq + 4 q̄q̄ + 9 qq̄
    end

    @testset "每个过程有 :type 和 :channels" begin
        for (proc, info) in SCATTERING_MESON_MAP
            @test haskey(info, :type)
            @test haskey(info, :channels)
            @test info[:type] in (:qq, :qqbar)
        end
    end

    @testset "qq 过程有 t/u 道" begin
        qq_procs = filter(kv -> kv.second[:type] === :qq, SCATTERING_MESON_MAP)
        @test !isempty(qq_procs)
        for (proc, info) in qq_procs
            ch = info[:channels]
            @test haskey(ch, :t)
            @test haskey(ch, :u)
            @test !haskey(ch, :s)
        end
    end

    @testset "qqbar 过程有 t/s 道" begin
        qqbar_procs = filter(kv -> kv.second[:type] === :qqbar, SCATTERING_MESON_MAP)
        @test !isempty(qqbar_procs)
        for (proc, info) in qqbar_procs
            ch = info[:channels]
            @test haskey(ch, :t)
            @test haskey(ch, :s)
            @test !haskey(ch, :u)
        end
    end

    @testset "通道结构字段完整" begin
        for (proc, info) in SCATTERING_MESON_MAP
            for (ch_name, ch_data) in info[:channels]
                @test haskey(ch_data, :simple)
                @test haskey(ch_data, :mixed_P)
                @test haskey(ch_data, :mixed_S)
                @test ch_data[:simple] isa Vector{Symbol}
                @test ch_data[:mixed_P] isa Bool
                @test ch_data[:mixed_S] isa Bool
            end
        end
    end

    @testset "介子种类合法" begin
        valid_mesons = Set([:pi, :K, :sigma_pi, :sigma_K])
        for (proc, info) in SCATTERING_MESON_MAP
            for (ch_name, ch_data) in info[:channels]
                for m in ch_data[:simple]
                    @test m in valid_mesons
                end
            end
        end
    end

    @testset "SCATTERING_PROCESS_KEYS 一致性" begin
        @test Set(SCATTERING_PROCESS_KEYS) == Set(keys(SCATTERING_MESON_MAP))
        @test length(SCATTERING_PROCESS_KEYS) == length(SCATTERING_MESON_MAP)
    end

    @testset "对称散射过程互补存在" begin
        # ud_to_ud 的电荷共轭 ubardbar_to_ubardbar 也应存在
        @test haskey(SCATTERING_MESON_MAP, :ud_to_ud)
        @test haskey(SCATTERING_MESON_MAP, :ubardbar_to_ubardbar)
        # uubar_to_uubar (q-qbar 自身散射)
        @test haskey(SCATTERING_MESON_MAP, :uubar_to_uubar)
    end
end
