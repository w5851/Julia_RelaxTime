---
title: DEFAULT_SKIP 清零转正记录（P0-1）
archived: true
original: docs/dev/active/2026-02-19_DEFAULT_SKIP转正记录.md
archived_date: 2026-02-19
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# DEFAULT_SKIP 清零转正记录（P0-1）

## 目标

为 `tests/unit/runtests.jl` 中历史 `DEFAULT_SKIP` 的 4 个用例建立逐项转正记录，并完成清零。

## 逐项结果

### 1) `test_core_integrals.jl`

- 现状：可执行通过（`29 passed, 1 broken`）。
- 失败根因：无致命失败；文件内存在 `@test_skip`（`Integrals.calculate_log_sum_derivatives` 未实现路径），属于“功能可选分支未实现”而非测试框架失败。
- 修复策略：从 `DEFAULT_SKIP` 转正进入 full profile；保留文件内 `@test_skip` 作为局部能力标记，后续单独跟踪该 API 实现。
- 回归点：`UNIT_FILES='pnjl/test_core_integrals.jl'` 可稳定通过。

### 2) `test_implicit_jacobian.jl`

- 现状：可执行通过（`12 passed`）。
- 失败根因：历史状态已过期；当前接口与实现可用。
- 修复策略：从 `DEFAULT_SKIP` 转正。
- 回归点：`UNIT_FILES='pnjl/test_implicit_jacobian.jl'` 可稳定通过。

### 3) `test_scans.jl`

- 现状：可执行通过（`33 passed`）。
- 失败根因：历史状态已过期；当前扫描路径可用。
- 修复策略：从 `DEFAULT_SKIP` 转正。
- 回归点：`UNIT_FILES='pnjl/test_scans.jl'` 可稳定通过。

### 4) `test_differential_cross_section.jl`

- 现状：可执行通过（`22 passed`）。
- 失败根因：历史“deprecated API”备注已过期；当前用例与实现一致。
- 修复策略：从 `DEFAULT_SKIP` 转正。
- 回归点：`UNIT_FILES='relaxtime/test_differential_cross_section.jl'` 可稳定通过。

## 执行变更

- `tests/unit/runtests.jl`：`DEFAULT_SKIP` 清零。
- `config/ci/unit_skip_policy.toml`：`max_skip=0`、`required_entries=[]`、phase 升级为 `phase-2`。

## 验证命令与结果

- `julia --project=. scripts/dev/check_unit_skip_policy.jl`
  - 结果：`[unit-skip-policy] OK`，`max_skip = 0`
- `$env:UNIT_PROFILE='full'; julia --project=. tests/unit/runtests.jl`
  - 结果：`100 passed, 1 broken, 0 failed`

## 备注

- `broken` 来自 `test_core_integrals.jl` 内部 `@test_skip`，不属于 `DEFAULT_SKIP`。
- 后续若需临时豁免测试，必须同步更新 `config/ci/unit_skip_policy.toml` 并附原因。
