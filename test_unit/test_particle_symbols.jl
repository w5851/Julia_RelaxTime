"""
测试 ParticleSymbols 工具模块

测试内容：
1. extract_flavor - 味标识解析
2. normalize_particle_symbol - 符号归一化
3. parse_scattering_process - 散射过程解析
4. parse_particle_pair - 粒子对解析
5. get_mass - 质量查询
6. get_chemical_potential - 化学势查询
7. get_wavefunction - 波函数查询
"""

using Test

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../src/utils"))

include("../src/utils/ParticleSymbols.jl")

using .ParticleSymbols

@testset "ParticleSymbols 工具模块测试" begin
    
    # ========== 测试1：味标识解析 ==========
    @testset "extract_flavor" begin
        println("\n" * "="^70)
        println("测试1：extract_flavor - 味标识解析")
        println("="^70)
        
        # 正粒子
        flavor, is_anti = extract_flavor(:u)
        @test flavor == :u
        @test is_anti == false
        println("  :u → (:u, false) ✓")
        
        flavor, is_anti = extract_flavor(:d)
        @test flavor == :d
        @test is_anti == false
        println("  :d → (:d, false) ✓")
        
        flavor, is_anti = extract_flavor(:s)
        @test flavor == :s
        @test is_anti == false
        println("  :s → (:s, false) ✓")
        
        # 反粒子
        flavor, is_anti = extract_flavor(:ubar)
        @test flavor == :u
        @test is_anti == true
        println("  :ubar → (:u, true) ✓")
        
        flavor, is_anti = extract_flavor(:dbar)
        @test flavor == :d
        @test is_anti == true
        println("  :dbar → (:d, true) ✓")
        
        flavor, is_anti = extract_flavor(:sbar)
        @test flavor == :s
        @test is_anti == true
        println("  :sbar → (:s, true) ✓")
    end
    
    # ========== 测试2：符号归一化 ==========
    @testset "normalize_particle_symbol" begin
        println("\n" * "="^70)
        println("测试2：normalize_particle_symbol - 符号归一化")
        println("="^70)
        
        # Unicode 反粒子符号
        @test normalize_particle_symbol("ū") == :ubar
        println("  'ū' → :ubar ✓")
        
        @test normalize_particle_symbol("đ") == :dbar
        println("  'đ' → :dbar ✓")
        
        @test normalize_particle_symbol("s̄") == :sbar
        println("  's̄' → :sbar ✓")
        
        # 标准符号
        @test normalize_particle_symbol("u") == :u
        println("  'u' → :u ✓")
        
        @test normalize_particle_symbol("ubar") == :ubar
        println("  'ubar' → :ubar ✓")
    end
    
    # ========== 测试3：散射过程解析 ==========
    @testset "parse_scattering_process" begin
        println("\n" * "="^70)
        println("测试3：parse_scattering_process - 散射过程解析")
        println("="^70)
        
        # qq 散射
        i, j, c, d = parse_scattering_process(:uu_to_uu)
        @test (i, j, c, d) == (:u, :u, :u, :u)
        println("  :uu_to_uu → (:u, :u, :u, :u) ✓")
        
        i, j, c, d = parse_scattering_process(:ud_to_ud)
        @test (i, j, c, d) == (:u, :d, :u, :d)
        println("  :ud_to_ud → (:u, :d, :u, :d) ✓")
        
        i, j, c, d = parse_scattering_process(:us_to_us)
        @test (i, j, c, d) == (:u, :s, :u, :s)
        println("  :us_to_us → (:u, :s, :u, :s) ✓")
        
        # qqbar 散射
        i, j, c, d = parse_scattering_process(:udbar_to_udbar)
        @test (i, j, c, d) == (:u, :dbar, :u, :dbar)
        println("  :udbar_to_udbar → (:u, :dbar, :u, :dbar) ✓")
        
        i, j, c, d = parse_scattering_process(:uubar_to_uubar)
        @test (i, j, c, d) == (:u, :ubar, :u, :ubar)
        println("  :uubar_to_uubar → (:u, :ubar, :u, :ubar) ✓")
        
        i, j, c, d = parse_scattering_process(:uubar_to_ssbar)
        @test (i, j, c, d) == (:u, :ubar, :s, :sbar)
        println("  :uubar_to_ssbar → (:u, :ubar, :s, :sbar) ✓")
    end
    
    # ========== 测试4：粒子对解析 ==========
    @testset "parse_particle_pair" begin
        println("\n" * "="^70)
        println("测试4：parse_particle_pair - 粒子对解析")
        println("="^70)
        
        # 标准格式
        p1, p2 = parse_particle_pair("uu")
        @test (p1, p2) == (:u, :u)
        println("  'uu' → (:u, :u) ✓")
        
        p1, p2 = parse_particle_pair("ud")
        @test (p1, p2) == (:u, :d)
        println("  'ud' → (:u, :d) ✓")
        
        p1, p2 = parse_particle_pair("udbar")
        @test (p1, p2) == (:u, :dbar)
        println("  'udbar' → (:u, :dbar) ✓")
        
        p1, p2 = parse_particle_pair("ubarū")
        @test (p1, p2) == (:ubar, :ubar)
        println("  'ubarsbar' → (:ubar, :sbar) ✓")
    end
    
    # ========== 测试5：参数查询 ==========
    @testset "参数查询" begin
        println("\n" * "="^70)
        println("测试5：参数查询 (质量、化学势、波函数)")
        println("="^70)
        
        # 创建测试参数
        quark_params = (
            m = (u = 1.52, d = 1.52, s = 2.50),
            μ = (u = 0.3, d = 0.3, s = 0.0)
        )
        
        # 测试5.1: get_mass
        println("\n质量查询:")
        m_u = get_mass(:u, quark_params)
        @test m_u ≈ 1.52
        println("  get_mass(:u) = $m_u fm⁻¹ ✓")
        
        m_ubar = get_mass(:ubar, quark_params)
        @test m_ubar ≈ 1.52  # 反粒子质量相同
        println("  get_mass(:ubar) = $m_ubar fm⁻¹ (反粒子质量相同) ✓")
        
        m_s = get_mass(:s, quark_params)
        @test m_s ≈ 2.50
        println("  get_mass(:s) = $m_s fm⁻¹ ✓")
        
        # 测试5.2: get_chemical_potential
        println("\n化学势查询:")
        μ_u = get_chemical_potential(:u, quark_params)
        @test μ_u ≈ 0.3
        println("  get_chemical_potential(:u) = $μ_u fm⁻¹ ✓")
        
        μ_ubar = get_chemical_potential(:ubar, quark_params)
        @test μ_ubar ≈ -0.3  # 反粒子化学势取负
        println("  get_chemical_potential(:ubar) = $μ_ubar fm⁻¹ (反粒子取负) ✓")
        
        μ_s = get_chemical_potential(:s, quark_params)
        @test μ_s ≈ 0.0
        println("  get_chemical_potential(:s) = $μ_s fm⁻¹ ✓")
        
        # 测试5.3: get_wavefunction
        println("\n波函数查询:")
        ψ_u = get_wavefunction(:u)
        @test ψ_u == [1.0, 0.0, 0.0]
        println("  get_wavefunction(:u) = [1.0, 0.0, 0.0] (列向量) ✓")
        
        ψbar_u = get_wavefunction(:ubar)
        @test ψbar_u == [1.0 0.0 0.0]
        println("  get_wavefunction(:ubar) = [1.0 0.0 0.0] (行向量) ✓")
        
        ψ_s = get_wavefunction(:s)
        @test ψ_s == [0.0, 0.0, 1.0]
        println("  get_wavefunction(:s) = [0.0, 0.0, 1.0] ✓")
    end
    
    # ========== 测试6：错误处理 ==========
    @testset "错误处理" begin
        println("\n" * "="^70)
        println("测试6：错误处理")
        println("="^70)
        
        # 无效格式
        @test_throws ErrorException parse_scattering_process(:invalid_format)
        println("  无效过程格式: 正确抛出错误 ✓")
        
        @test_throws ErrorException parse_particle_pair("xyz")
        println("  无效粒子符号: 正确抛出错误 ✓")
        
        quark_params = (
            m = (u = 1.52, d = 1.52, s = 2.50),
            μ = (u = 0.3, d = 0.3, s = 0.0)
        )
        
        @test_throws ErrorException get_mass(:unknown, quark_params)
        println("  未知味标识: 正确抛出错误 ✓")
    end
    
end

println("\n" * "="^70)
println("ParticleSymbols 工具模块测试完成！")
println("="^70)
