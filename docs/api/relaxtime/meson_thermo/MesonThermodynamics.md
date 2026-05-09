# MesonThermodynamics

## 定位

`src/relaxtime/MesonThermodynamics.jl` 提供介子热力学的 pressure 核。

它是 `MesonDensity` 的 pressure 对应层，当前只做最小补位，不改写已有 meson density / meson mass 主链。

## 已实现能力

- `stable_meson_pressure`
- `stable_meson_pressure_summary`
- `strict_bw_meson_pressure`
- `strict_bw_meson_pressure_summary`
- `phase_shift_meson_pressure`
- `phase_shift_meson_pressure_summary`

## 物理口径

### Stable

稳定粒子极限采用标准玻色 pressure 积分：

```math
P_M = d_M \int \frac{dq\,q^2}{2\pi^2}
T \log(1 - e^{-(E_M-\mu_M)/T})^{-1}.
```

### Reduced strict BW

第一版只支持 reduced strict BW：

- 读取 workflow 已产出的 `mass/gamma`
- 不在 pressure 层重求 `q` 依赖复极点

### Phase-shift / generalized BU reference

第一版复用现有：

- propagator
- polarization
- weighted phase

并支持两种 scheme：

- `:current`
- `:gbu_reference`

当前已显式拆出：

- `QP`（timelike / `omega >= q`）部分
- `LD`（Landau damping / `omega < q`）部分

并支持最小治理参数：

- `ld_cutoff`
- `ld_cutoff_mode`
- `ld_threshold_mode`

## 输出摘要

summary helper 统一返回：

- `P_pi`
- `P_K`
- `P_meson`
- `P_K_over_P_pi`

phase-shift 摘要还会携带：

- `scheme`
- `P_pi_qp`, `P_pi_ld`
- `P_K_qp`, `P_K_ld`
- `P_meson_qp`, `P_meson_ld`
- `ld_cutoff`
- `ld_cutoff_mode`
- `ld_threshold_mode`
- `pi_pressure`
- `k_pressure`

## 当前边界

- 只承诺 `π/K`
- 尚未覆盖 scalar partner / pseudoscalar nonet 扩张
- 尚未提供 validation / regression baseline 资产
