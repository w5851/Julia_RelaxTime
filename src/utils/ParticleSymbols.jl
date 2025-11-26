"""
# ParticleSymbols.jl

统一的粒子符号处理工具模块，提供符号解析、参数查询等通用功能。

## 功能概述

1. **味标识解析**: 从符号中提取味类型和反粒子信息
2. **过程解析**: 解析散射过程符号 (如 :uu_to_uu → (:u, :u, :u, :u))
3. **参数查询**: 根据味标识获取质量、化学势、波函数等

## 设计原则

- 统一命名规范: :u, :d, :s (正粒子), :ubar, :dbar, :sbar (反粒子)
- 支持多种输入格式: "ū" → :ubar, "đ" → :dbar, "s̄" → :sbar
- 提供清晰的错误信息
- 避免代码重复

## 使用示例

```julia
using .ParticleSymbols

# 味标识解析
flavor, is_antiparticle = extract_flavor(:ubar)  # (:u, true)

# 过程解析
i, j, c, d = parse_scattering_process(:uu_to_uu)  # (:u, :u, :u, :u)

# 参数查询
m = get_mass(:u, quark_params)  # 1.52 fm⁻¹
μ = get_chemical_potential(:ubar, quark_params)  # -0.3 fm⁻¹ (反粒子)
ψ = get_wavefunction(:u, false)  # [1.0, 0.0, 0.0] (列向量)
```
"""
module ParticleSymbols

# 导出函数
export extract_flavor
export normalize_particle_symbol
export parse_scattering_process
export parse_particle_pair
export get_mass
export get_chemical_potential
export get_wavefunction

# 导入必要的常量
include(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
using .Constants_PNJL: ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s

# ----------------------------------------------------------------------------
# 1. 味标识解析
# ----------------------------------------------------------------------------

"""
    extract_flavor(particle::Symbol) -> (flavor::Symbol, is_antiparticle::Bool)

从粒子符号中提取味类型和反粒子标志。

# 参数

- `particle::Symbol`: 粒子符号 (:u, :d, :s, :ubar, :dbar, :sbar)

# 返回

元组 `(flavor, is_antiparticle)`:
- `flavor::Symbol`: 味类型 (:u, :d, :s)
- `is_antiparticle::Bool`: 是否为反粒子

# 示例

```julia
extract_flavor(:u)     # (:u, false)
extract_flavor(:ubar)  # (:u, true)
extract_flavor(:dbar)  # (:d, true)
```
"""
function extract_flavor(particle::Symbol)
    s = string(particle)
    
    if endswith(s, "bar")
        # 反粒子: 移除 "bar" 后缀
        flavor = Symbol(s[1:end-3])
        return (flavor, true)
    else
        # 正粒子
        return (particle, false)
    end
end

"""
    normalize_particle_symbol(s::AbstractString) -> Symbol

将各种粒子符号格式归一化为标准符号。

# 支持的输入格式

- 正粒子: "u", "d", "s" → :u, :d, :s
- 反粒子: "ū", "đ", "s̄" → :ubar, :dbar, :sbar
- 反粒子: "ubar", "dbar", "sbar" → :ubar, :dbar, :sbar

# 参数

- `s::AbstractString`: 粒子符号字符串

# 返回

标准化的粒子符号 (Symbol)

# 示例

```julia
normalize_particle_symbol("ū")     # :ubar
normalize_particle_symbol("ubar")  # :ubar
normalize_particle_symbol("u")     # :u
```
"""
function normalize_particle_symbol(s::AbstractString)
    # 替换 Unicode 反粒子标记为 "bar" 后缀
    s = replace(s, "ū" => "ubar")
    s = replace(s, "đ" => "dbar")
    s = replace(s, "s̄" => "sbar")
    
    return Symbol(s)
end

# ----------------------------------------------------------------------------
# 2. 散射过程解析
# ----------------------------------------------------------------------------

"""
    parse_scattering_process(process::Symbol) -> (i::Symbol, j::Symbol, c::Symbol, d::Symbol)

解析散射过程符号，提取初末态粒子。

# 参数

- `process::Symbol`: 散射过程符号 (格式: :ab_to_cd)

# 返回

元组 `(particle_i, particle_j, particle_c, particle_d)`:
- 初态粒子: i, j
- 末态粒子: c, d

# 支持的格式

- :uu_to_uu → (:u, :u, :u, :u)
- :udbar_to_udbar → (:u, :dbar, :u, :dbar)
- :uubar_to_ssbar → (:u, :ubar, :s, :sbar)

# 示例

```julia
i, j, c, d = parse_scattering_process(:uu_to_uu)
# 返回: (:u, :u, :u, :u)

i, j, c, d = parse_scattering_process(:udbar_to_udbar)
# 返回: (:u, :dbar, :u, :dbar)
```

# 错误处理

- 如果格式不符合 "ab_to_cd"，抛出错误
- 如果粒子符号无法解析，抛出错误
"""
function parse_scattering_process(process::Symbol)
    # 转换为字符串并分割 "_to_"
    process_str = string(process)
    parts = split(process_str, "_to_")
    
    if length(parts) != 2
        error("Invalid process format: $process (expected format: 'ab_to_cd')")
    end
    
    initial_str = String(parts[1])
    final_str = String(parts[2])
    
    # 解析初态和末态粒子对
    particle_i, particle_j = parse_particle_pair(initial_str)
    particle_c, particle_d = parse_particle_pair(final_str)
    
    return (particle_i, particle_j, particle_c, particle_d)
end

"""
    parse_particle_pair(pair_str::AbstractString) -> (p1::Symbol, p2::Symbol)

从字符串中解析粒子对。

# 参数

- `pair_str::AbstractString`: 粒子对字符串 (如 "uu", "udbar", "uū")

# 返回

元组 `(particle1, particle2)`

# 示例

```julia
parse_particle_pair("uu")     # (:u, :u)
parse_particle_pair("udbar")  # (:u, :dbar)
parse_particle_pair("uū")     # (:u, :ubar)
```
"""
function parse_particle_pair(pair_str::AbstractString)
    # 归一化 Unicode 反粒子符号
    pair_str = replace(pair_str, "ū" => "ubar")
    pair_str = replace(pair_str, "đ" => "dbar")
    pair_str = replace(pair_str, "s̄" => "sbar")
    
    # 解析粒子序列
    particles = Symbol[]
    
    i = 1
    while i <= length(pair_str)
        # 检查是否是反粒子 (xxxbar)
        if i + 3 <= length(pair_str) && pair_str[i:i+3] == "ubar"
            push!(particles, :ubar)
            i += 4
        elseif i + 3 <= length(pair_str) && pair_str[i:i+3] == "dbar"
            push!(particles, :dbar)
            i += 4
        elseif i + 3 <= length(pair_str) && pair_str[i:i+3] == "sbar"
            push!(particles, :sbar)
            i += 4
        elseif pair_str[i] in ['u', 'd', 's']
            # 正粒子
            push!(particles, Symbol(pair_str[i]))
            i += 1
        else
            error("Unknown particle symbol at position $i in '$pair_str'")
        end
    end
    
    if length(particles) != 2
        error("Expected 2 particles in '$pair_str', found $(length(particles))")
    end
    
    return (particles[1], particles[2])
end

# ----------------------------------------------------------------------------
# 3. 参数查询
# ----------------------------------------------------------------------------

"""
    get_mass(particle::Symbol, quark_params::NamedTuple) -> Float64

获取指定粒子的质量。

# 参数

- `particle::Symbol`: 粒子符号 (:u, :d, :s, :ubar, :dbar, :sbar)
- `quark_params::NamedTuple`: 夸克参数 (需包含 `m` 字段)

# 返回

质量 [fm⁻¹]

# 说明

反粒子与正粒子质量相同。

# 示例

```julia
m_u = get_mass(:u, quark_params)      # 1.52 fm⁻¹
m_ubar = get_mass(:ubar, quark_params)  # 1.52 fm⁻¹ (相同)
```
"""
function get_mass(particle::Symbol, quark_params::NamedTuple)
    flavor, _ = extract_flavor(particle)
    
    if flavor == :u
        return quark_params.m.u
    elseif flavor == :d
        return quark_params.m.d
    elseif flavor == :s
        return quark_params.m.s
    else
        error("Unknown flavor: $flavor")
    end
end

"""
    get_chemical_potential(particle::Symbol, quark_params::NamedTuple) -> Float64

获取指定粒子的化学势。

# 参数

- `particle::Symbol`: 粒子符号 (:u, :d, :s, :ubar, :dbar, :sbar)
- `quark_params::NamedTuple`: 夸克参数 (需包含 `μ` 字段)

# 返回

化学势 [fm⁻¹]

# 说明

反粒子化学势为负: μ(q̄) = -μ(q)

# 示例

```julia
μ_u = get_chemical_potential(:u, quark_params)      # +0.3 fm⁻¹
μ_ubar = get_chemical_potential(:ubar, quark_params)  # -0.3 fm⁻¹
```
"""
function get_chemical_potential(particle::Symbol, quark_params::NamedTuple)
    flavor, is_antiparticle = extract_flavor(particle)
    
    # 获取正粒子化学势
    if flavor == :u
        μ = quark_params.μ.u
    elseif flavor == :d
        μ = quark_params.μ.d
    elseif flavor == :s
        μ = quark_params.μ.s
    else
        error("Unknown flavor: $flavor")
    end
    
    # 反粒子化学势取负
    return is_antiparticle ? -μ : μ
end

"""
    get_wavefunction(particle::Symbol) -> Union{Vector{Float64}, Matrix{Float64}}

获取指定粒子的味空间波函数。

# 参数

- `particle::Symbol`: 粒子符号 (:u, :d, :s, :ubar, :dbar, :sbar)

# 返回

- 正粒子: 列向量 (3×1 Vector)
- 反粒子: 行向量 (1×3 Matrix)

# 物理意义

- ψ_u = [1, 0, 0]: u夸克在味空间的表示
- ψbar_u = [1 0 0]: u反夸克在味空间的表示

# 示例

```julia
ψ_u = get_wavefunction(:u)      # [1.0, 0.0, 0.0]
ψbar_u = get_wavefunction(:ubar)  # [1.0 0.0 0.0]
```
"""
function get_wavefunction(particle::Symbol)
    flavor, is_antiparticle = extract_flavor(particle)
    
    if is_antiparticle
        # 反粒子: 行向量
        if flavor == :u
            return ψbar_u
        elseif flavor == :d
            return ψbar_d
        elseif flavor == :s
            return ψbar_s
        else
            error("Unknown flavor: $flavor")
        end
    else
        # 正粒子: 列向量
        if flavor == :u
            return ψ_u
        elseif flavor == :d
            return ψ_d
        elseif flavor == :s
            return ψ_s
        else
            error("Unknown flavor: $flavor")
        end
    end
end

end  # module ParticleSymbols
