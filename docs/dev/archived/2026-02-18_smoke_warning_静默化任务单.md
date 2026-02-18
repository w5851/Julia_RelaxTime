---
title: smoke warning 静默化任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_18_smoke_warning_静默化任务单.md
archived_date: 2026-02-18
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# smoke warning 静默化任务单（可勾选）

## 文档目的

在不改变功能行为的前提下，收敛 `UNIT_PROFILE=smoke` 输出中的可避免 warning 噪声，提升 CI 日志可读性。

---

## 任务清单

### W1 重复常量定义 warning 收敛

- [x] 将 models smoke 测试中的 `const _MODELS_ENTRY` 改为局部 `_models_entry`
- [x] 保持 `Models.jl` 按需 include 逻辑不变

### W2 模块重复 include warning 收敛

- [x] 为 `TransportWorkflow`/`RelaxationTime` 相关 smoke 测试增加 `isdefined(Main, ...)` include guard
- [x] 为 `GaussLegendre` / `Constants_PNJL` 在 smoke 入口测试中增加 include guard

### W3 运行验证

- [x] 全量执行 `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
- [x] 验证结果保持全绿

---

## 结果摘要

- 最终 smoke：`690 passed, 0 failed, 0 errored`。
- 已消除此前主要噪声：
  - `WARNING: redefinition of constant Main._MODELS_ENTRY`
  - `WARNING: replacing module TransportWorkflow`
  - `WARNING: replacing module Constants_PNJL`
  - `WARNING: replacing module GaussLegendre`
- 仍保留 1 条领域性提示：
  - `quark_params missing :A under ξ≠0; auto-building anisotropic A via A_aniso`
  - 该提示来自功能路径本身（非重复加载噪声），用于提示自动补全 `A` 字段。

---

## 执行记录

- 2026-02-18：完成 warning 来源定位与分组修复。
- 2026-02-18：完成 smoke 全量复跑并确认全绿。
- 2026-02-18：任务单完成，准备归档。
