# Rotation 变体 API

`RotationModel` 通过 `Models` 统一接口提供旋转背景的模型适配，`solve_rotation_point` 完成单点求解与可观测量整理。当前为最小实现，具体物理范围与公式见 [Rotation-PNJL 公式](../../../../reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md)。

## 调用与单位

在仓库根目录运行：

```julia
include("src/models/Models.jl")
model = Models.create_model(:Rotation)
result = Models.solve_rotation_point(0.80, 0.20; omega=0.05)
```

`T_fm`、`mu_fm` 与 `omega` 使用 `fm⁻¹`。返回字段 `pressure`、`energy`、`omega_potential` 使用 `fm⁻⁴`，`rho` 与 `entropy` 使用 `fm⁻³`；`rho` 是工作流将夸克净密度除以 3 后的重子密度。实现见 [RotationWorkflow](../../../../../src/models/variants/rotation/workflows/RotationWorkflow.jl)。

## 实现职责

- `RotationModel` 对接模型接口，包括 `solve_gap`、`omega_components` 和 `number_densities`。
- `RotationWorkflow` 将求解与热力学输出组织为单点结果。
- `src/models/variants/rotation/core/RotationThermo.jl` 承载数值核，由模型适配层与工作流复用。

调用方通过 `Models` 使用这些能力。核心函数供实现层复用，不作为单点计算的默认入口。

[自动生成的公开导出索引](generated/Exports.md)列出本主题的完整导出表面。
