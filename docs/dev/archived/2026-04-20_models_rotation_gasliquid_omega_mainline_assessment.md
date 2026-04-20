---
title: Rotation / GasLiquid 回归通用 Ω 主线评估（#95 本地草案）
archived: true
original: docs/dev/active/2026-04-20_models_rotation_gasliquid_omega_mainline_assessment.md
archived_date: 2026-04-20
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Rotation / GasLiquid 回归通用 Ω 主线评估（#95 本地草案）

> 日期：2026-04-20
> 对应 issue：#95
> 目的：从“公式项级语义”而非“当前工程实现形态”判断 Rotation / GasLiquid 能否回归 `omega.jl` 通用组装主线。
> 状态：主线回归已完成并验证，见 commits `d7adcec` / `d96fd64`。

---

## 1) 判定口径

- 通用主线定义：`omega_components(model, ...) = chi + poly + vac + therm`，见 `src/models/omega.jl:21`。
- 判定分级：
  - `可直接接入`：当前实现已能无语义损失映射到四项。
  - `条件接入`：公式可映射，但现实现含占位/工程近似，需先收敛前置项。
  - `暂不接入`：当前实现与公式关键定义存在未收敛差异，直接回归会掩盖物理口径风险。

---

## 2) Rotation 模型评估

### 2.1 公式侧项级结构

- 公式给出：
  - `chi = G <qq>^2`（`docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md:37`）
  - `poly = U(Phi, PhiBar; T)`（`docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md:70`）
  - `therm = 旋转模和 + 双对数核`（`docs/reference/formula/models/rotation/Rotation_PNJL_CoreEquations.md:39`）
  - `vac` 在文档最小口径中未单独解析，可并入 therm/vac 的工程分拆。

### 2.2 代码侧项级映射

- 当前实现已经返回四项结构：
  - `rotation_omega_components`：`src/models/variants/rotation/RotationModel.jl:96`
  - `chi = G*phi^2`：`src/models/variants/rotation/core/RotationThermo.jl:154`
  - `poly` 多项式势：`src/models/variants/rotation/core/RotationThermo.jl:106`
  - `therm` 由旋转核积分构成：`src/models/variants/rotation/core/RotationThermo.jl:123`, `src/models/variants/rotation/core/RotationThermo.jl:156`

### 2.3 风险与占位

- `RotationModel` 顶层接口仍保留占位函数（不走核心积分核）：
  - `calculate_chiral = 0.05*phi^2`：`src/models/variants/rotation/RotationModel.jl:64`
  - `vacuum_contribution = -0.04*masses[1]`：`src/models/variants/rotation/RotationModel.jl:71`
- 但其 `omega_components` 已直接委托 `RotationThermo.omega_components`，因此主 Ω 计算口径并不依赖上述占位函数。

### 2.4 结论

- **判定：`条件接入`（偏“接近可直接接入”）**。
- 理由：公式项级与通用四项分解可一一映射；主要阻碍在于顶层占位 API 尚未清理，易引发“不同入口不同口径”。

### 2.5 最小改造路径

1. 让 `calculate_chiral/vacuum_contribution/thermal_contribution` 与 `RotationThermo.omega_components` 同源，消除占位常数。
2. 明确 `vac` 在 Rotation 口径中的分拆约定（是否独立项或并入热核），并写入模型注释与文档。
3. 增加一组固定点验证：直接四项求和与 `omega_components.omega` 数值一致（rtol/atol 明确）。

---

## 3) GasLiquid 模型评估

### 3.1 公式侧项级结构

- 文档给出 RMF 口径：
  - `Omega_RMF = U(sigma) - vector fields - thermal log integral`（`docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md:101`）
  - 与通用分解可映射为：
    - `chi`: 标量势（含 sigma / delta）
    - `vac`: 向量场静态项（omega/rho）
    - `therm`: 费米双对数积分
    - `poly`: 0（非 Polyakov 模型）

### 3.2 代码侧项级映射

- `GasLiquidThermodynamics.omega_components` 已按上述思路分拆：
  - `chi`：`src/models/variants/gas_liquid/core/Thermodynamics.jl:56`
  - `vac`：`src/models/variants/gas_liquid/core/Thermodynamics.jl:57`
  - `therm`：`src/models/variants/gas_liquid/core/Thermodynamics.jl:58`
  - `poly = 0`：`src/models/variants/gas_liquid/core/Thermodynamics.jl:59`
- `GasLiquidModel.omega_components` 也已直接委托上述实现：`src/models/variants/gas_liquid/GasLiquidModel.jl:86`。

### 3.3 风险与占位

- 顶层适配层仍有与四项分拆不一致的占位实现：
  - `calculate_chiral = 0`：`src/models/variants/gas_liquid/GasLiquidModel.jl:54`
  - `vacuum_contribution = 0`：`src/models/variants/gas_liquid/GasLiquidModel.jl:64`
  - `thermal_contribution` 通过 `sigma=delta=0` 构造临时状态，仅返回 `-pressure`：`src/models/variants/gas_liquid/GasLiquidModel.jl:69`
- 公式文档也已明确“当前实现包含占位”风险：`docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md:243`。

### 3.4 结论

- **判定：`条件接入`（当前更接近“暂不建议切主入口”）**。
- 理由：底层 `omega_components` 的项级语义基本可映射，但顶层通用 term 函数仍存在占位/简化，若强行回归到“只走通用 term 组装”，可能放大口径不一致风险。

### 3.5 最小改造路径

1. 统一 `GasLiquidModel` 顶层 term 函数与 `GasLiquidThermodynamics.omega_components` 的来源，禁止混用占位返回。
2. 固化 `mu_p/mu_n` 与 `rho03` 的职责边界（质量分裂 vs 化学势分裂），避免把约束口径差异掺入 Ω 分项。
3. 为一个固定参数集建立 `P=-Omega` 与分项守恒回归点（至少 smoke+core 两档）。

---

## 4) 统一结论（#95 首轮）

- Rotation：`条件接入`（去占位后可基本转为直接接入）。
- GasLiquid：`条件接入`（需先清理顶层占位并锁定参数/约束口径）。
- 两者都不是“公式天然不可映射”；阻碍主要在“当前实现层次存在多入口占位与语义分叉”。

---

## 5) 建议验证清单（用于 #95 后续）

- [ ] V1：`omega_components.omega` 与 `chi+vac+therm+poly` 点值一致性（Rotation/GasLiquid 各 3 个点）
- [ ] V2：`P=-Omega` 点值一致性（FixedMu）
- [ ] V3：变体模型入口与通用入口的结果一致性（同一初值、同一积分参数）
- [ ] V4：将验证脚本接入 `tests/integration/` 的 smoke profile
