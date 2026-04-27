---
description: 数值基线回归治理（基线生成/存储/校验/CI门禁/变更准入）
---

# 基线回归治理

## 目录规范
- 基线文件：`tests/baselines/<domain>/baseline_<feature>_<scope>_vN.csv`
- 导出脚本：`scripts/dev/export_<feature>_baseline.jl`
- 校验测试：`tests/unit/<domain>/test_<feature>_baseline_smoke.jl`

## CSV 列设计
- 至少含：点位标识列（T, mu, xi）、被比较量（eta, sigma, zeta）。
- 小数格式固定（`%.6f` / `%.16e`）；列顺序固定；不写时间戳/路径。

## 容差与校验
- 容差在任务文档中声明（`rtol` / `atol`），测试与文档一致。
- 使用 `isapprox(actual, expected; rtol=..., atol=...)`。
- 同时保留"实现间一致性"与"基线校验"两条线。

## CI 分层
- **smoke（PR 必过）**：关键小规模固定点，并入 unit-smoke。
- **nightly**：大网格、更多点位放 nightly，与 smoke 同口径不同规模。

## 基线变更准入（必须同时提交）
1. 变更原因（算法修复/物理口径修正/参数更新）
2. 影响范围
3. 新旧差异摘要（固定点级别）
4. 验证命令与结果
5. 对应文档更新

**禁止：** 仅更新基线文件不更新说明/测试；未说明原因的大幅数值漂移直接"刷基线"。
