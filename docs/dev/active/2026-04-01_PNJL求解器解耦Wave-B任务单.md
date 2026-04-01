# PNJL求解器解耦 Wave-B 任务单

> 执行主线说明：本任务单用于 Wave-B 的唯一执行主线（勾选、证据、验收）。
> 设计与实现参考：
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-b-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-b-implementation-plan.md`

## 1. 目标

- [ ] 将 Wave-A 已冻结契约推进到扫描/单点共用执行路径。
- [ ] 准备兼容入口 deprecation 节奏（不在本波次强删）。
- [ ] 固化 fallback/selection 诊断稳定性回归。

## 2. 范围

### 2.1 本期范围

- [ ] scan 侧候选治理接线收敛。
- [ ] 兼容入口统一路由与迁移状态补齐。
- [ ] fallback 原因标签规范化。

### 2.2 非范围

- [ ] 不移除 `solve_fixed*` 兼容入口。
- [ ] 不替换求解后端。

## 3. 任务分解

### B1：scan 侧治理并轨

- [ ] 为 scan 治理并轨补失败回归测试。
- [ ] 让 `TmuScan/TrhoScan` 走共享 governance 接口。
- [ ] 验证代表点无非预期漂移。

### B2：兼容路由 deprecation-ready

- [ ] 为统一入口/兼容入口等价性补失败 smoke。
- [ ] 完善迁移状态与查询接口。
- [ ] 保持外部签名兼容。

### B3：fallback 诊断稳定性

- [ ] 为 fallback reason 标签补失败回归。
- [ ] 统一 reason tag 与输出顺序。
- [ ] 验证诊断输出可稳定复现。

### B4：文档治理同步

- [ ] 更新迁移映射状态与删除门槛。
- [ ] 回填验证证据到本任务单。
- [ ] 通过 docs/governance 检查。

## 4. 验收标准

- [ ] scan 与单点使用统一治理口径。
- [ ] compatibility 路由行为可追溯且保持非破坏。
- [ ] fallback 诊断回归稳定。
- [ ] smoke + docs/governance 通过。

## 5. 验证命令

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`

## 6. DoD

- [ ] 任务项与验收项全部勾选。
- [ ] 关键验证命令可复现通过。
- [ ] 迁移状态文档同步更新。

## 7. 执行记录

- [x] 2026-04-01：基于 Wave-A 完成状态创建 Wave-B 任务单。
