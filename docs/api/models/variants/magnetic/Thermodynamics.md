# Magnetic 热力学接口

本页说明 magnetic 主题中真正面向业务计算的热力学主接口。实现位于 [src/models/pnjl_physics/core/MagneticThermodynamics.jl](../../../../src/models/pnjl_physics/core/MagneticThermodynamics.jl#L1)。

## 主要导出

- `coupling_GB`
- `calculate_magnetic_omega_components`
- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`
- `calculate_magnetic_rho`
- `calculate_magnetic_number_densities`
- `magnetic_nmax_convergence_report`

## `coupling_GB`

计算 IMC 耦合 `G(B)`。它主要服务于热力学组件计算，而不是独立业务工作流。

## `calculate_magnetic_omega_components`

返回 magnetic 主题里最完整的热力学拆分结果，包括：

- `chi`
- `poly`
- `vac`
- `therm`
- `masses`
- `omega`
- `n_max`
- `G_B`

这是理解磁场下 `omega` 组装边界的核心入口。

## `calculate_magnetic_omega` / `calculate_magnetic_pressure`

这两个是最常见的业务层接口：

- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`

若你只需要最终热力学量，优先使用它们，而不是直接调用 Landau 低层函数。

`calculate_magnetic_pressure=-Omega` 只表示固定外部磁场背景下的标量物质压力；当前
没有 Maxwell 自能、磁化强度或横向/纵向压力张量输出。

## `calculate_magnetic_rho`

通过数值导数路径计算 flavor 相关密度。它更偏进阶入口，使用时应注意步长与数值敏感区的影响。

## `calculate_magnetic_number_densities`

返回：

- `quark`
- `net`
- `antiquark = nothing`（磁场路线不单独输出夸克/反夸克两支）
- `baryon`

其中 `quark` 是为兼容历史脚本保留的字段名，正式语义与 `net` 相同，表示
`q - qbar` 净密度；`baryon` 由三味净密度求和得到。该返回值不满足普通 PNJL
`(quark, antiquark)` 的独立分量合同，不应直接传入需要独立反夸克占据数的输运路径。
这是受支持正 `eB` 的 magnetic core 语义；输入低于
`MAGNETIC_EB_MIN_FM2` 会明确拒绝，不会转回普通 PNJL。

因此 `PNJLMagneticModel` 会将通用
`supports_number_densities` capability 标为 `false`，避免该结果被普通输运路径误接；
这不撤销专用 `calculate_magnetic_number_densities` API。
适合磁场固定点、扫描脚本或回归基线中直接消费。

## `magnetic_nmax_convergence_report`

这是 magnetic 主题必须显式强调的收敛治理入口。它比较 `n_base` 与 `n_base + delta_n` 下的 `omega` 差异，并给出：

- `converged`
- `rtol`
- `rel_diff`
- `n_base`
- `n_probe`
- `omega_base`
- `omega_probe`

建议在以下场景优先调用：

- 新参数区首次跑点
- 调整 `eB_fm2`、`p_num` 或 `pz_max` 后
- 更新实现后做 baseline 或 regression 检查前
