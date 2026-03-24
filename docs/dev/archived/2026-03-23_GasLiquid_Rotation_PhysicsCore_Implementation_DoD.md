---
title: Gas-Liquid + Rotation 物理内核实现任务单（本 PR）
archived: true
original: docs/dev/active/2026-03-23_GasLiquid_Rotation_PhysicsCore_Implementation_DoD.md
archived_date: 2026-03-24
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Gas-Liquid + Rotation 物理内核实现任务单（本 PR）

更新日期：2026-03-23

> 目标确认：结合仓库通用能力与两份公式文档，完成 `gas_liquid` 与 `rotation` 的**真实物理内核**实现（替换当前占位核），并保持 `Models` 统一入口、配置驱动、测试分层和文档一致性。

---

## 1. 范围与非范围

### 1.1 本 PR 范围（必须完成）

- [x] 用文献对齐公式替换 `gas_liquid` 占位内核：
  - `Omega_RMF`、`rho_s/rho_B`、平均场闭合方程、`P/rho/s/epsilon`。
- [x] 用文献对齐公式替换 `rotation` 占位内核：
  - `epsilon_n`、`n` 模求和、Bessel 权重、PNJL 双对数热项、`Omega_rot` 与 gap 驻点。
- [x] 建立配置驱动参数层：
  - `config/models/gas_liquid/` 与 `config/models/rotation/`。
- [x] 接入 `Models` 统一入口与 workflow，保持 API 契约不破坏。
- [x] 增补 unit + integration smoke 测试覆盖，并通过核心门禁。

### 1.2 本 PR 非范围（明确不做）

- [ ] 不在本轮引入 GPU/MPI。
- [ ] 不在本轮实现高阶导数全套（如所有 `dP/domega^n`）。
- [ ] 不在本轮做大规模性能极限优化（仅做必要数值稳定与明显低效修正）。

---

## 2. 现状差距（基于当前代码）

### 2.1 Gas-Liquid

- [x] `solve_equilibrium` 仍是 `tanh` 目标函数驱动，不是 RMF 自洽方程。
- [x] `pressure_density_entropy_energy` 为代数占位表达，不是 `Omega_RMF` 导出。
- [x] 参数不完整：当前只有 `g_sigma/g_delta` 等少量项，缺 `g_omega/g_rho` 与介子质量集。

### 2.2 Rotation

- [x] `solve_gap` 仍是代数占位残差，不是 `∂Omega/∂(phi,Phi,PhiBar)=0`。
- [x] 热力学核缺少 `n` 求和与 Bessel 权重积分结构。
- [x] Polyakov 势与文献参数化未真正进入内核计算。

---

## 3. 实施任务分解（可勾选）

## 3.1 Gas-Liquid 物理核落地

- [x] **参数结构重建**
  - 文件：`src/models/variants/gas_liquid/core/EquationSet.jl`
  - 补齐 `m_N, m_sigma, m_omega, m_rho(,m_delta), f_sigma, f_omega, f_rho, b, c`。
  - 明确 `f_x = g_x / m_x` 到 `g_x` 的还原路径与单位换算。

- [x] **密度与热核实现**
  - 文件：`src/models/variants/gas_liquid/core/Thermodynamics.jl`
  - 实现 `f_i/fbar_i`, `rho_i`, `rho_s,i`, `Omega_RMF`, `P/rho/s/epsilon`。
  - 支持 `rho_B = rho_p + rho_n` 主口径，兼容层处理历史 `/2`。

- [x] **自洽求解器替换**
  - 文件：`src/models/variants/gas_liquid/core/EquationSet.jl`
  - 用 `R_sigma/R_omega/R_rho(/R_delta)` 替代占位残差。
  - 保留数值收敛诊断字段（残差范数、收敛标志）。

- [x] **模型壳对齐**
  - 文件：`src/models/variants/gas_liquid/GasLiquidModel.jl`
  - 清理占位 `vacuum_contribution/thermal_contribution` 路径，改由统一 `omega_components` 真实返回。

## 3.2 Rotation 物理核落地

- [x] **参数与势函数实现**
  - 文件：`src/models/variants/rotation/core/RotationThermo.jl`
  - 落地文献参数：`m, Lambda, G, a0..a3, b3, b4, T0, r, n_cut`。
  - 实现 Polyakov 多项式势 `U(Phi,PhiBar;T)`。

- [x] **巨势积分核实现**
  - 文件：`src/models/variants/rotation/core/RotationThermo.jl`
  - 实现 `epsilon_n`、`W_n=J_n^2+J_{n+1}^2`、`n` 模求和与 `p_t/p_z` 双积分（含 cutoff）。

- [x] **gap 方程真实化**
  - 文件：`src/models/variants/rotation/RotationModel.jl`
  - 用 `∂Omega/∂phi=0`, `∂Omega/∂Phi=0`, `∂Omega/∂PhiBar=0` 替代占位 `tanh` 残差。

- [x] **workflow 输出一致化**
  - 文件：`src/models/variants/rotation/workflows/RotationWorkflow.jl`
  - 输出字段保持现契约，同时来源改为真实内核值。

## 3.3 配置与文档对齐

- [x] **参数配置文件新增**
  - 目录：`config/models/gas_liquid/`, `config/models/rotation/`
  - 明确每个键：符号来源、单位、默认值、换算关系。

- [ ] **文档/API 同步**
  - 公式文档只做最小必要更新（若实现口径有细节差异）。
  - 稳定入口更新 `docs/api/`（如新增/调整对外字段或配置入口）。

## 3.4 测试与门禁

- [x] **unit：内核级验证**
  - Gas-Liquid：密度、残差收敛、热力学恒等式检查。
  - Rotation：`omega->0` 连续性、`omega*r<1` 约束、`n_cut` 收敛性 smoke。

- [x] **integration：workflow smoke**
  - 保持 `solve_gas_liquid_point` 与 `solve_rotation_point` 可用，结果有限且结构稳定。

- [x] **必跑门禁**
  - `UNIT_PROFILE=smoke`
  - `INTEGRATION_PROFILE=smoke`
  - `scripts/dev/check_pnjl_migration_guard.jl`
  - `scripts/dev/check_docs_consistency.jl`

---

## 4. DoD（Definition of Done）

仅当以下全部满足，才算本 PR 完成：

- [x] `gas_liquid` 与 `rotation` 不再依赖占位代数公式（代码层可追溯到文献方程结构）。
- [x] 参数由配置驱动，且单位/符号映射在文档中可追溯。
- [x] `Models` 统一入口与现有 workflow 契约不破坏（兼容现有调用）。
- [x] unit + integration smoke 全通过，迁移守卫与文档一致性检查通过。
- [x] 开发文档更新并可归档（包含验证命令与结果摘要）。

---

## 5. 风险与回退

- [ ] 风险：旋转双积分 + 模求和带来计算耗时增加。
  - 预案：先实现可配置 `n_cut` 与积分网格，测试使用小网格 smoke。
- [ ] 风险：Gas-Liquid 参数映射导致收敛敏感。
  - 预案：保留初值策略与收敛诊断，必要时增加 continuation seed。
- [ ] 风险：新内核改变现有测试数值。
  - 预案：分层更新测试断言（结构断言与物理趋势断言分离）。

---

## 6. 执行顺序建议

- [x] 先落 Gas-Liquid 参数与内核（复杂度可控，先收敛一套配置链路）。
- [x] 再落 Rotation 内核与 gap 驻点（积分核更复杂，放第二阶段）。
- [x] 最后做统一测试收口与文档/API 对齐。

---

## 7. 验证记录（本次执行）

- [x] `julia --project=. -e 'include("tests/unit/models/test_gas_liquid_model.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_rotation_model.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_gas_liquid_workflow.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_rotation_workflow.jl")'`
- [x] `julia --project=. -e 'include("tests/integration/models/test_gas_liquid_workflow_smoke.jl")'`
- [x] `julia --project=. -e 'include("tests/integration/models/test_rotation_workflow_smoke.jl")'`
- [x] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_pnjl_migration_guard.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
