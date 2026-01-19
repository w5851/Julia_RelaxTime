# MesonMass

介子质量计算接口，基于极点方程 $1 - 4K\Pi(p_0)=0$，并支持 $p_0=M+i\Gamma/2$ 的宽度耦合。

## 模块
- `src/relaxtime/MesonMass.jl`

## 主要函数

### `ensure_quark_params_has_A(quark_params, thermo_params; p_nodes=32)`
为 `quark_params` 补齐 `A` 字段（若缺失）。

**参数**
- `quark_params::NamedTuple`：包含 `m`、`μ` 字段
- `thermo_params::NamedTuple`：包含 `T`、`Φ`、`Φbar`
- `p_nodes::Int`：A 积分节点数（默认 32）

**返回**
- 新的 `quark_params`（包含 `A` 字段）

### `meson_mass_equation(meson, k0, gamma, k_norm, quark_params, thermo_params, K_coeffs)`
返回复数残差 $f=1-4K\Pi$。

**参数**
- `meson::Symbol`：`:pi`, `:K`, `:sigma_pi`, `:sigma_K`
- `k0::Float64`：能量实部
- `gamma::Float64`：宽度
- `k_norm::Float64`：三动量模长
- `quark_params::NamedTuple`
- `thermo_params::NamedTuple`
- `K_coeffs::NamedTuple`

**返回**
- `ComplexF64`

### `solve_meson_mass(meson, quark_params, thermo_params; K_coeffs=nothing, k_norm=0.0, initial_mass=nothing, initial_gamma=0.0, nlsolve_kwargs...)`
求解介子质量与宽度。

**返回**
- `NamedTuple`：`(mass, gamma, converged, residual_norm, solution)`

## 使用示例
```julia
using .MesonMass: solve_meson_mass

quark_params = (m=(u=0.3, d=0.3, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

res = solve_meson_mass(:pi, quark_params, thermo_params)
@show res.mass res.gamma res.converged
```