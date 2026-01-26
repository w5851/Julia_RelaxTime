"""
# ParticleSymbols Module

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

## Dual Interface Pattern

This module supports **both struct and NamedTuple parameters** for parameter lookup functions:

```julia
# Using structs (recommended)
using Main.ParameterTypes: QuarkParams

q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
m = get_mass(:u, q)  # 1.52 fm⁻¹
μ = get_chemical_potential(:ubar, q)  # +0.3 fm⁻¹

# Using NamedTuples (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
m = get_mass(:u, q_nt)  # 1.52 fm⁻¹
μ = get_chemical_potential(:ubar, q_nt)  # +0.3 fm⁻¹
```

Both produce identical results. Internal normalization ensures type stability.

## 使用示例

```julia
using .ParticleSymbols

# 味标识解析
flavor, is_antiparticle = extract_flavor(:ubar)  # (:u, true)

# 过程解析
i, j, c, d = parse_scattering_process(:uu_to_uu)  # (:u, :u, :u, :u)

# 参数查询（支持结构体和NamedTuple）
m = get_mass(:u, quark_params)  # 1.52 fm⁻¹
μ = get_chemical_potential(:ubar, quark_params)  # +0.3 fm⁻¹ (同号)
ψ = get_wavefunction(:u, false)  # [1.0, 0.0, 0.0] (列向量)
```
"""
module ParticleSymbols

# 导出函数
export extract_flavor
export extract_quark_flavor
export is_antiquark
export normalize_particle_symbol
export parse_scattering_process
export parse_particle_pair
export parse_particle_pair_str
export parse_scattering_process_flavors
export parse_scattering_process_flavor_codes
export FLAVOR_U, FLAVOR_D, FLAVOR_S
export get_quark_masses_for_process
export get_mass
export get_chemical_potential
export get_wavefunction

# 导入必要的常量
include(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
using .Constants_PNJL: ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s, SCATTERING_MESON_MAP

# 导入参数类型
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "..", "ParameterTypes.jl"))
end
using Main.ParameterTypes: QuarkParams, as_namedtuple

# ----------------------------------------------------------------------------
# 0. Normalization Helper
# ----------------------------------------------------------------------------

"""
    _nt_quark(q) -> NamedTuple

Internal normalization helper that converts QuarkParams struct to NamedTuple.
If input is already a NamedTuple, returns it unchanged.

This ensures consistent internal representation regardless of input type.
"""
@inline _nt_quark(q) = q isa QuarkParams ? as_namedtuple(q) : q

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
    # Hot-path: avoid string allocations for common particles.
    if particle === :u || particle === :d || particle === :s
        return (particle, false)
    elseif particle === :ubar
        return (:u, true)
    elseif particle === :dbar
        return (:d, true)
    elseif particle === :sbar
        return (:s, true)
    end

    # Fallback for uncommon/extended symbols.
    s = string(particle)
    if endswith(s, "bar")
        flavor = Symbol(s[1:end-3])
        return (flavor, true)
    end
    return (particle, false)
end

"""
    extract_quark_flavor(particle::Symbol) -> Symbol

从粒子符号中提取夸克味（忽略反粒子 `bar` 标记）。

等价于 `first(extract_flavor(particle))`，但提供单独命名用于可读性。
"""
@inline function extract_quark_flavor(particle::Symbol)::Symbol
    flavor, _ = extract_flavor(particle)
    return flavor
end

"""
    is_antiquark(particle::Symbol) -> Bool

判断给定粒子符号是否为反夸克（`:ubar/:dbar/:sbar`）。
"""
@inline function is_antiquark(particle::Symbol)::Bool
    return particle === :ubar || particle === :dbar || particle === :sbar
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
    cached = get(_PARSED_SCATTERING_PROCESS, process, nothing)
    if cached !== nothing
        return cached
    end
    return _parse_scattering_process_uncached(process)
end

"""
    parse_scattering_process_flavors(process::Symbol) -> (f1, f2, f3, f4)

解析散射过程并返回四个外腿的“夸克味” `(f1,f2,f3,f4)`，忽略反粒子 `bar` 标记。

示例：
`parse_scattering_process_flavors(:uubar_to_ssbar) == (:u,:u,:s,:s)`
"""
@inline function parse_scattering_process_flavors(process::Symbol)::NTuple{4, Symbol}
    cached = get(_PARSED_SCATTERING_PROCESS_FLAVORS, process, nothing)
    if cached !== nothing
        return cached
    end
    q1, q2, q3, q4 = parse_scattering_process(process)
    return (
        extract_quark_flavor(q1),
        extract_quark_flavor(q2),
        extract_quark_flavor(q3),
        extract_quark_flavor(q4),
    )
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
    return parse_particle_pair_str(pair_str)
end

"""
    parse_particle_pair_str(pair_str::AbstractString) -> (p1::Symbol, p2::Symbol)

将形如 `"ud"`, `"uubar"`, `"ubarsbar"`, `"dū"` 的二粒子字符串解析为两个粒子 `Symbol`。

支持 token：
- 夸克：`u`, `d`, `s`
- 反夸克：`ubar`, `dbar`, `sbar`
- Unicode overbar：`ū`, `đ`, `s̄`（会先归一化为 `*bar`）
"""
function parse_particle_pair_str(pair_str::AbstractString)::Tuple{Symbol, Symbol}
    s = replace(String(pair_str), "ū" => "ubar", "đ" => "dbar", "s̄" => "sbar")

    local p1::Symbol
    local p2::Symbol
    n = 0
    i = firstindex(s)
    while i <= lastindex(s)
        rest = SubString(s, i)
        if startswith(rest, "ubar")
            n += 1
            if n == 1
                p1 = :ubar
            elseif n == 2
                p2 = :ubar
            else
                error("Invalid particle pair '$s' (expected 2 particles)")
            end
            i = nextind(s, i, 4)
        elseif startswith(rest, "dbar")
            n += 1
            if n == 1
                p1 = :dbar
            elseif n == 2
                p2 = :dbar
            else
                error("Invalid particle pair '$s' (expected 2 particles)")
            end
            i = nextind(s, i, 4)
        elseif startswith(rest, "sbar")
            n += 1
            if n == 1
                p1 = :sbar
            elseif n == 2
                p2 = :sbar
            else
                error("Invalid particle pair '$s' (expected 2 particles)")
            end
            i = nextind(s, i, 4)
        else
            c = s[i]
            if c == 'u'
                n += 1
                if n == 1
                    p1 = :u
                elseif n == 2
                    p2 = :u
                else
                    error("Invalid particle pair '$s' (expected 2 particles)")
                end
            elseif c == 'd'
                n += 1
                if n == 1
                    p1 = :d
                elseif n == 2
                    p2 = :d
                else
                    error("Invalid particle pair '$s' (expected 2 particles)")
                end
            elseif c == 's'
                n += 1
                if n == 1
                    p1 = :s
                elseif n == 2
                    p2 = :s
                else
                    error("Invalid particle pair '$s' (expected 2 particles)")
                end
            else
                error("Invalid particle token at position $i in '$s'")
            end
            i = nextind(s, i)
        end
    end
    n == 2 || error("Invalid particle pair '$s' (expected 2 particles)")
    return (p1, p2)
end

# ----------------------------------------------------------------------------
# 2.x 散射过程解析缓存
# ----------------------------------------------------------------------------

function _parse_scattering_process_uncached(process::Symbol)::NTuple{4, Symbol}
    s = string(process)
    parts = split(s, "_to_")
    length(parts) == 2 || error("Invalid process format: $process (expected format: 'ab_to_cd')")
    q1, q2 = parse_particle_pair_str(parts[1])
    q3, q4 = parse_particle_pair_str(parts[2])
    return (q1, q2, q3, q4)
end

const _PARSED_SCATTERING_PROCESS = let d = Dict{Symbol, NTuple{4, Symbol}}()
    for p in keys(SCATTERING_MESON_MAP)
        d[p] = _parse_scattering_process_uncached(p)
    end
    d
end

const _PARSED_SCATTERING_PROCESS_FLAVORS = let d = Dict{Symbol, NTuple{4, Symbol}}()
    for (p, (q1, q2, q3, q4)) in _PARSED_SCATTERING_PROCESS
        d[p] = (
            extract_quark_flavor(q1),
            extract_quark_flavor(q2),
            extract_quark_flavor(q3),
            extract_quark_flavor(q4),
        )
    end
    d
end

# ----------------------------------------------------------------------------
# 2.x 夸克味编码（性能敏感路径）
# ----------------------------------------------------------------------------

"""
夸克味编码（isbits）：用于在热路径中避免 `Symbol` 元组参与 hash/逃逸。

- `FLAVOR_U == 0x01` 表示 `u`
- `FLAVOR_D == 0x02` 表示 `d`
- `FLAVOR_S == 0x03` 表示 `s`
"""
const FLAVOR_U::UInt8 = 0x01
const FLAVOR_D::UInt8 = 0x02
const FLAVOR_S::UInt8 = 0x03

@inline function _hex_u8(x::UInt8)::String
    return "0x" * lpad(string(x, base=16), 2, '0')
end

@inline function _unknown_flavor_code_error(code::UInt8)
    error(
        "Unknown flavor code: $(_hex_u8(code)) (expected $(_hex_u8(FLAVOR_U))=u, $(_hex_u8(FLAVOR_D))=d, $(_hex_u8(FLAVOR_S))=s)"
    )
end

@inline function _flavor_code(flavor::Symbol)::UInt8
    if flavor === :u
        return FLAVOR_U
    elseif flavor === :d
        return FLAVOR_D
    elseif flavor === :s
        return FLAVOR_S
    else
        error("Unknown quark flavor: $flavor (expected :u, :d, or :s)")
    end
end

const _PARSED_SCATTERING_PROCESS_FLAVOR_CODES = let d = Dict{Symbol, NTuple{4, UInt8}}()
    for (p, (f1, f2, f3, f4)) in _PARSED_SCATTERING_PROCESS_FLAVORS
        d[p] = (_flavor_code(f1), _flavor_code(f2), _flavor_code(f3), _flavor_code(f4))
    end
    d
end

@inline function _parse_scattering_process_flavor_codes(process::Symbol)::NTuple{4, UInt8}
    cached = get(_PARSED_SCATTERING_PROCESS_FLAVOR_CODES, process, (0x00, 0x00, 0x00, 0x00))
    if cached[1] != 0x00
        return cached
    end
    f1, f2, f3, f4 = parse_scattering_process_flavors(process)
    return (_flavor_code(f1), _flavor_code(f2), _flavor_code(f3), _flavor_code(f4))
end

"""
    parse_scattering_process_flavor_codes(process::Symbol) -> (c1,c2,c3,c4)

解析散射过程并返回四个外腿的“夸克味代码” `(c1,c2,c3,c4)`，其中

- `FLAVOR_U == 0x01` 表示 `u`
- `FLAVOR_D == 0x02` 表示 `d`
- `FLAVOR_S == 0x03` 表示 `s`

该接口用于性能敏感路径：返回值为 isbits 的 `NTuple{4,UInt8}`，可避免 `Symbol` 元组在微基准/热循环中发生逃逸分配。
"""
@inline function parse_scattering_process_flavor_codes(process::Symbol)::NTuple{4, UInt8}
    return _parse_scattering_process_flavor_codes(process)
end

# ----------------------------------------------------------------------------
# 3. 参数查询
# ----------------------------------------------------------------------------

"""
    get_mass(particle::Symbol, quark_params) -> Float64

获取指定粒子的质量。

# 参数

- `particle::Symbol`: 粒子符号 (:u, :d, :s, :ubar, :dbar, :sbar)
- `quark_params`: 夸克参数，可以是 `QuarkParams` 结构体或包含 `m` 字段的 `NamedTuple`

# 返回

质量 [fm⁻¹]

# 说明

反粒子与正粒子质量相同。

# 示例

```julia
# Using QuarkParams struct (recommended)
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
m_u = get_mass(:u, q)      # 1.52 fm⁻¹
m_ubar = get_mass(:ubar, q)  # 1.52 fm⁻¹ (相同)

# Using NamedTuple (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
m_u = get_mass(:u, q_nt)      # 1.52 fm⁻¹
```
"""
function get_mass(particle::Symbol, quark_params::Union{NamedTuple, QuarkParams})
    quark_params = _nt_quark(quark_params)
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

@inline function _mass_for_quark_flavor(flavor::Symbol, quark_params)::Float64
    quark_params = _nt_quark(quark_params)
    if flavor === :u
        return quark_params.m.u
    elseif flavor === :d
        return quark_params.m.d
    elseif flavor === :s
        return quark_params.m.s
    else
        error("Unknown flavor: $flavor")
    end
end

@inline function _mass_for_flavor_code(code::UInt8, quark_params)::Float64
    quark_params = _nt_quark(quark_params)
    if code == FLAVOR_U
        return quark_params.m.u
    elseif code == FLAVOR_D
        return quark_params.m.d
    elseif code == FLAVOR_S
        return quark_params.m.s
    else
        _unknown_flavor_code_error(code)
    end
end

"""
    get_quark_masses_for_process(process::Symbol, quark_params) -> NTuple{4, Float64}

根据散射过程符号返回四个外腿粒子的质量 `(m1, m2, m3, m4)`（单位：fm⁻¹）。

约定：
- 夸克与反夸克质量相同（因此会先用 `extract_quark_flavor` 忽略 `bar` 标记）。
- `process` 的解析由 `parse_scattering_process` 完成（对常见过程有缓存）。

该函数刻意只做“符号→质量索引”的工具工作，不涉及任何动力学或散射道（s/t/u）逻辑。
"""
function get_quark_masses_for_process(process::Symbol, quark_params::Union{NamedTuple, QuarkParams})::NTuple{4, Float64}
    quark_params = _nt_quark(quark_params)
    c1, c2, c3, c4 = _parse_scattering_process_flavor_codes(process)
    return (
        _mass_for_flavor_code(c1, quark_params),
        _mass_for_flavor_code(c2, quark_params),
        _mass_for_flavor_code(c3, quark_params),
        _mass_for_flavor_code(c4, quark_params),
    )
end

"""
    get_chemical_potential(particle::Symbol, quark_params) -> Float64

获取指定粒子的化学势。

# 参数

- `particle::Symbol`: 粒子符号 (:u, :d, :s, :ubar, :dbar, :sbar)
- `quark_params`: 夸克参数，可以是 `QuarkParams` 结构体或包含 `μ` 字段的 `NamedTuple`

# 返回

化学势 [fm⁻¹]

# 说明

本仓库约定：该函数 **始终返回“夸克化学势” μ_q（同号）**，不对反夸克取负。
反夸克占据数的差异应由分布函数接口（如 `antiquark_distribution*` 或 `sign_flag`）处理。

# 示例

```julia
# Using QuarkParams struct (recommended)
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
μ_u = get_chemical_potential(:u, q)      # +0.3 fm⁻¹
μ_ubar = get_chemical_potential(:ubar, q)  # +0.3 fm⁻¹

# Using NamedTuple (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
μ_u = get_chemical_potential(:u, q_nt)      # +0.3 fm⁻¹
μ_ubar = get_chemical_potential(:ubar, q_nt)  # +0.3 fm⁻¹
```
"""
function get_chemical_potential(particle::Symbol, quark_params::Union{NamedTuple, QuarkParams})
    quark_params = _nt_quark(quark_params)
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
    
    # 约定：始终返回“夸克化学势” μ_q（同号），不对反夸克取负。
    return μ
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
