# Magnetic PNJL 参数说明

本文档给出外磁场 PNJL 的参数口径、默认值与建议范围。

## 一、运行参数

### 1) 状态向量 `x_state`

- 类型：`SVector{5}`
- 顺序：`(φ_u, φ_d, φ_s, Φ, Φbar)`
- 单位：`φ_i` 为 `fm⁻³` 口径下的平均场变量（沿用现有 PNJL 约定），`Φ/Φbar` 无量纲。

### 2) 化学势向量 `mu_vec`

- 类型：`SVector{3}`
- 顺序：`(μ_u, μ_d, μ_s)`
- 单位：`fm⁻¹`

### 3) 温度 `T_fm`

- 单位：`fm⁻¹`
- 从 MeV 换算：`T_fm = T_MeV / ħc_MeV_fm`

### 4) 外磁场 `eB_fm2`

- 单位：`fm⁻²`
- 从 MeV² 换算：`eB_fm2 = eB_MeV2 / ħc_MeV_fm^2`

## 二、`MagneticConfig` 字段

### `eB_fm2::Float64`

外磁场强度。`|eB_fm2| <= 1e-14` 时自动退化到零磁场 PNJL 路径。

### `n_max::Union{Int,Nothing}`

Landau 层截断。

- `nothing`：自动估计（推荐）
- 指定 `Int`：手动固定层数

### `p_num::Int`

`p_z` 积分节点数。

- smoke 推荐：`24`
- 默认生产：`96`
- nightly 可按资源提高

### `pz_max::Float64`

`p_z` 积分上界（`fm⁻¹`）。

- `0.0`：自动取默认上界
- 手动建议：`10~30`（视点位与精度要求）

### `cutoff_N::Int`

平滑截断指数。

- 默认：`10`
- 建议保持不变，除非做方法学对比

### `imc::MagneticIMCParams`

`G(B)` 的 IMC 参数。

## 三、IMC 参数 `MagneticIMCParams`

`G(B) = G0 * (1 + aζ² + bζ³) / (1 + cζ² + dζ⁴)`，其中 `ζ = eB / Λ_QCD²`。

默认值：

- `a = 0.108805`
- `b = -1.0133e-4`
- `c = 0.02228`
- `d = 1.84558e-4`
- `Λ_QCD_MeV = 300.0`

## 四、收敛判据（`n_max`）

统一使用：

- `magnetic_nmax_convergence_report(...; delta_n=6, rtol=3e-2)`

含义：

- 基于 `n_base` 与 `n_base+6` 的 `Ω` 相对差
- 判定阈值：`rel_diff <= 3e-2`

## 五、基线与容差

### smoke 基线

- 文件：`tests/baselines/pnjl/baseline_pnjl_magnetic_fixedpoints_smoke_v1.csv`
- 测试：`tests/unit/pnjl/test_magnetic_fixedpoint_baseline_smoke.jl`
- 容差：`rtol=8e-2`, `atol=1e-6`

### nightly 基线

- 文件：`tests/baselines/pnjl/baseline_pnjl_magnetic_fixedpoints_nightly_v1.csv`
- 门禁脚本：`scripts/pnjl/run_magnetic_fixedpoint_regression.jl`
- 容差：`rtol=5e-2`, `atol=1e-7`

## 六、配置文件口径

磁场配置模板：`config/models/pnjl/magnetic_default.toml`

字段映射：

- `[magnetic]`
  - `eB_fm2`
  - `n_max`
  - `p_num`
  - `pz_max`
  - `cutoff_N`
- `[magnetic.imc]`
  - `a`, `b`, `c`, `d`, `Lambda_QCD_MeV`

## 七、推荐使用顺序

1. 用 `default_magnetic_config` 跑单点。
2. 用 `magnetic_nmax_convergence_report` 检查截断收敛。
3. 用 `run_magnetic_eb_scan.jl` 做固定 `(T,μ)` 的 eB 扫描。
4. 修改实现后先过 smoke，再跑 nightly 回归。
