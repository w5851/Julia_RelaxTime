---
name: baseline-regression-governance
description: 管理数值基线回归的全流程规范（基线生成、存储、测试校验、CI门禁与变更准入），适用于固定点/采样点回归、防止版本更新引入数值漂移。关键词：baseline, regression, fixed points, rtol, atol, smoke, nightly, PNJL, transport
---

# baseline-regression-governance

## 1. 技能目标

本技能用于在仓库内建立并维护“数值基线回归治理”机制，确保：

1. 关键功能有固定点/采样点基线；
2. 版本更新不会静默改变关键数值结果；
3. 基线变更有明确理由与可追溯证据；
4. smoke 与 nightly 的回归分层清晰。

---

## 2. 何时使用

当出现以下任一场景时，应启用本技能：

- 新增模型功能（如 PNJL 求解、输运系数、相图扫描）；
- 重构/优化数值内核（积分器、求解器、缓存、并行）；
- 调整物理参数口径、单位换算、默认配置；
- 需要将“文档约定”落地为“自动化回归门禁”。

---

## 3. 目录与文件规范（必须遵守）

### 3.1 基线文件目录

- 统一放在：`tests/baselines/<domain>/`
- 示例：
	- `tests/baselines/relaxtime/baseline_transport_fixedpoints_v1.csv`

### 3.2 校验测试目录

- 放在对应测试域：`tests/unit/<domain>/`
- 示例：
	- `tests/unit/relaxtime/test_transport_legacy_models_bridge_smoke.jl`

### 3.3 基线导出脚本目录

- 放在：`scripts/dev/`
- 示例：
	- `scripts/dev/export_transport_fixedpoint_baseline.jl`

### 3.4 命名规则

- 基线文件：`baseline_<feature>_<scope>_vN.csv`
	- 例：`baseline_transport_fixedpoints_v1.csv`
- 导出脚本：`export_<feature>_baseline.jl`
	- 例：`export_transport_fixedpoint_baseline.jl`
- 测试文件：`test_<feature>_baseline_smoke.jl`（或并入现有 smoke）

---

## 4. 基线列设计规范

基线 CSV 至少应包含：

1. 点位标识列（如 `T, mu, xi`）；
2. 被比较量（如 `eta, sigma, zeta`）；
3. 可稳定复现，不依赖运行时随机性。

建议：

- 小数格式固定（如 `%.6f` / `%.16e`）；
- 列顺序固定，避免无意义 diff；
- 避免写入时间戳、路径等易变字段。

---

## 5. 容差与校验规范

### 5.1 容差来源

- 容差必须在任务文档或验收矩阵中声明（例如 `rtol=8e-2`, `atol=1e-6`）；
- 测试实现必须与文档口径一致。

### 5.2 校验方式

- 使用 `isapprox(actual, expected; rtol=..., atol=...)`；
- 每个固定点至少比较核心输出量（例如 `eta/sigma/zeta`）；
- 同时保留“实现间一致性”校验（如 `legacy` vs `models`）与“基线校验”两条线。

---

## 6. CI 分层策略

### 6.1 smoke（PR 必过）

- 放入关键、小规模、快速固定点基线回归；
- 推荐并入现有 `unit-smoke`，避免过多碎片 workflow。

### 6.2 nightly（全量/重负载）

- 大网格、更多点位、更重计算放 nightly；
- 与 smoke 保持“同口径不同规模”。

---

## 7. 基线变更准入规则（必须）

任何基线文件改动都必须同时提交：

1. 变更原因（算法修复、物理口径修正、参数更新）；
2. 影响范围（哪些模块、哪些指标）；
3. 新旧差异摘要（至少固定点级别）；
4. 验证命令与结果（本地或 CI）；
5. 对应文档更新（README/STATUS/active 任务单）。

禁止：

- 仅更新基线文件，不更新说明与测试；
- 未说明原因的大幅数值漂移直接“刷基线”。

---

## 8. 执行流程（标准作业）

1. 确定固定点与比较量；
2. 编写/更新导出脚本到 `scripts/dev/`；
3. 生成基线到 `tests/baselines/<domain>/`；
4. 在 `tests/unit/<domain>/` 写入基线校验；
5. 跑定向测试确认通过；
6. 接入（或复用）smoke CI；
7. 在 active 任务单记录证据与命令。

---

## 9. 当前项目落地示例（已存在）

- 基线文件：`tests/baselines/relaxtime/baseline_transport_fixedpoints_v1.csv`
- 导出脚本：`scripts/dev/export_transport_fixedpoint_baseline.jl`
- 校验测试：`tests/unit/relaxtime/test_transport_legacy_models_bridge_smoke.jl`
- 定向测试命令：
	- `UNIT_FILES='relaxtime/test_transport_legacy_models_bridge_smoke.jl'`
	- `julia --project=. tests/unit/runtests.jl`

---

## 10. 面向未来扩展（PNJL 等）

新增其他功能基线时，按同一模式扩展：

- 目录：`tests/baselines/pnjl/`、`tests/baselines/<new-domain>/`
- 文件：`baseline_pnjl_<feature>_v1.csv`
- 脚本：`scripts/dev/export_pnjl_<feature>_baseline.jl`
- 测试：`tests/unit/pnjl/test_<feature>_baseline_smoke.jl`

核心原则：

- 先有文档口径，再有自动化测试；
- 基线是“质量契约”，不是“临时输出”。