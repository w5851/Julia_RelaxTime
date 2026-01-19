# MottTransition

Mott 相变判据与差值计算接口。

## 模块
- `src/relaxtime/MottTransition.jl`

## 主要函数

### `mott_threshold_mass(meson, quark_params)`
返回介子对应的组分夸克阈值质量。

**参数**
- `meson::Symbol`：`:pi`, `:K`, `:sigma_pi`, `:sigma_K`
- `quark_params::NamedTuple`

**返回**
- `Float64`

### `mott_gap(meson, meson_mass, quark_params)`
返回 $\Delta = M_{meson} - (M_{q_1}+M_{q_2})$。

**参数**
- `meson::Symbol`
- `meson_mass::Float64`
- `quark_params::NamedTuple`

**返回**
- `Float64`

### `is_mott_point(meson, meson_mass, quark_params; atol=1e-6)`
判断是否满足 Mott 相变条件。

**返回**
- `Bool`

## 使用示例
```julia
using .MottTransition: mott_gap, is_mott_point

quark_params = (m=(u=0.3, d=0.3, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
gap = mott_gap(:K, 0.8, quark_params)
@show gap is_mott_point(:K, 0.8, quark_params)
```