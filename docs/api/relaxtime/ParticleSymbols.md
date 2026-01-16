# ParticleSymbols 模块 API 文档

## 模块定位

`ParticleSymbols` 是目前 **relaxtime 链路专用** 的“粒子/过程符号解析工具”。它负责把形如 `:uu_to_uu`、`"uubar"`、`"dū"` 的过程/粒子对表示解析为统一的粒子 `Symbol` 表达（`:u/:d/:s/:ubar/:dbar/:sbar`）。

之所以放在 `docs/api/relaxtime/`：当前只有 relaxtime 依赖它，并且其缓存策略（基于 `SCATTERING_MESON_MAP` 的预填充）也与 relaxtime 的过程集合强相关。

## 设计目标

- 统一入口：避免在多个模块里重复写 `split/replace/collect` 解析逻辑。
- 性能友好：热点路径尽量 **零/低分配**（避免构造 `Vector{Symbol}`、避免频繁 `string(...)`）。
- 兼容输入：支持 ASCII `ubar/dbar/sbar` 以及 Unicode overbar（`ū/đ/s̄`）。

## 核心 API

### `parse_scattering_process(process)`

```julia
parse_scattering_process(process::Symbol) -> (q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol)
```

解析散射过程 key，例如：

- `:uu_to_uu` → `(:u, :u, :u, :u)`
- `:uubar_to_ssbar` → `(:u, :ubar, :s, :sbar)`

**缓存策略**

- 对 `keys(SCATTERING_MESON_MAP)` 中出现的过程进行预解析并缓存。
- 其它过程走 uncached 路径（仍支持解析，但会有 `split` 等分配）。

### `parse_particle_pair_str(pair_str)` / `parse_particle_pair(pair_str)`

```julia
parse_particle_pair_str(pair_str::AbstractString) -> (p1::Symbol, p2::Symbol)
parse_particle_pair(pair_str::AbstractString) -> (p1::Symbol, p2::Symbol)
```

解析二粒子 token 字符串：

- `"ud"` → `(:u, :d)`
- `"uubar"` → `(:u, :ubar)`
- `"dū"` → `(:d, :ubar)`（会先归一化）

`parse_particle_pair` 目前是对 `parse_particle_pair_str` 的兼容 wrapper；建议新代码直接使用 `parse_particle_pair_str`。

### `extract_flavor(particle)` / `extract_quark_flavor(particle)`

```julia
extract_flavor(particle::Symbol) -> (flavor::Symbol, is_antiparticle::Bool)
extract_quark_flavor(particle::Symbol) -> Symbol
```

- `extract_flavor(:ubar)` → `(:u, true)`
- `extract_quark_flavor(:ubar)` → `:u`

内部对 `:u/:d/:s/:ubar/:dbar/:sbar` 做了 `Symbol` 快速分支，以减少热点路径字符串开销。

### `is_antiquark(particle)`

```julia
is_antiquark(particle::Symbol) -> Bool
```

判断粒子符号是否为反夸克（`:ubar/:dbar/:sbar`）。

说明：该函数只做“符号层”的判断，不涉及化学势正负约定（μ 的约定见下文）。

### `parse_scattering_process_flavors(process)`

```julia
parse_scattering_process_flavors(process::Symbol) -> (f1, f2, f3, f4)
```

返回四个外腿粒子的“夸克味”元组，忽略 `bar` 标记：

- `:uubar_to_ssbar` → `(:u, :u, :s, :s)`

**缓存策略**

- 对 `keys(SCATTERING_MESON_MAP)` 中出现的过程预先缓存其 flavor 元组（避免热点路径重复调用 `extract_*`）。
- 其它过程仍可解析，但会走非缓存路径。

### `get_quark_masses_for_process(process, quark_params)`

```julia
get_quark_masses_for_process(process::Symbol, quark_params::NamedTuple) -> (m1, m2, m3, m4)
```

给定散射过程 `process`（例如 `:uubar_to_ssbar`），返回四个外腿粒子的质量。

- 该函数会先用 `parse_scattering_process` 得到 `(q1,q2,q3,q4)`，再用 `extract_quark_flavor` 忽略 `bar` 标记。
- 约定：夸克与反夸克质量相同。

## 使用位置（relaxtime）

- `TotalPropagator`：用于 `parse_scattering_process` 与 `extract_quark_flavor`（质量/味因子/混合介子味结构）。
- `TotalCrossSection`：用于从 `process::Symbol` 解析四个粒子（从而确定末态统计因子等）。

## 注意事项

- 本模块的 `get_chemical_potential` 约定为：**始终返回“夸克化学势” μ_q（同号）**，反夸克的负号由分布函数（如 `antiquark_distribution*`）处理。
