# Solver Phase 3 能力边界固化与治理门禁任务单（稳定契约与防回耦）

## 1. 目标

- 在 Phase 1（主链收敛）与 Phase 2（职责拆分）基础上，完成 solver 能力边界固化。
- 建立稳定对外契约（Result/Diagnostic/Error taxonomy），杜绝上层对 solver 私有实现细节的穿透依赖。
- 建立“防回耦”治理门禁：新增功能可演进，但依赖方向与契约不再漂移。

### 1.1 前置依赖

- 必须在以下文档 DoD 达成后执行：
  - `docs/dev/active/2026-04-08_solver_phase1主链收敛任务单.md`
  - `docs/dev/active/2026-04-08_solver_phase2职责拆分降耦任务单.md`

## 2. 范围与非范围

### 2.1 本期范围（Phase 3）

- [x] 发布 solver 稳定契约（类型 + 序列化视图 + 兼容策略）。
- [x] 建立 phase/scans/workflows 适配层，仅消费稳定契约。
- [x] 增加耦合守卫测试与治理脚本检查点。
- [x] 清理过渡期字段与旁路参数（在兼容期策略允许范围内）。

### 2.2 非范围（本期不做）

- [ ] 不在本期改动物理模型公式、扫描策略阈值和科学结论口径。
- [ ] 不做大规模 UI/前端重构（仅做契约消费适配）。
- [ ] 不做跨仓库接口治理（仅限当前仓库 `src/models` 域）。

## 3. 稳定契约定义（对外）

### 3.1 结果契约（Result Contract）

- [x] 统一 `SolverResult`/raw result 的稳定字段集合（必须字段、可选字段、弃用字段）。
- [x] 明确字段语义与单位（例如 `rho_norm`, `residual_norm`, `xi`, `pressure`）。
- [x] 定义版本化策略：`contract_version`（建议从 `v1` 起）。

### 3.2 诊断契约（Diagnostic Contract）

- [x] 统一 summary/full 语义与字段最小集合。
- [x] 诊断字段分层：
  - 稳定公共字段（上层可依赖）
  - 内部调试字段（禁止上层依赖）
- [x] 提供兼容输出视图（NamedTuple）与类型化对象双通道。

### 3.3 错误契约（Error Taxonomy）

- [x] 定义错误类别（参数错误、约束失败、求解失败、中断等）。
- [x] 统一错误码/错误消息模板，避免上层按字符串关键词做脆弱判断。

## 4. 文件级任务分解（可执行）

### Batch-A：契约发布与文档化

- [x] A1 新增/更新 `docs/api/` 中 solver 契约文档（Result/Diagnostic/Error）。
- [x] A2 在 `src/models/solver` 增加契约导出入口（建议集中在 `Solver.jl` 或新契约模块）。
- [x] A3 更新 `src/models/Models.jl` 导出清单，避免导出私有实现符号。

### Batch-B：上层适配层收敛

- [x] B1 `src/models/phase/PMPhaseDiagnostic.jl`
  - 通过适配器读取稳定诊断字段，不再直接枚举 solver 私有字段。
- [x] B2 `src/models/scans/TmuScan.jl`
  - 统一从稳定结果契约读取成功判据与摘要字段。
- [x] B3 `src/models/scans/TrhoScan.jl`
  - 同步适配稳定契约，移除对临时兼容字段的硬依赖。
- [x] B4 `src/models/workflows/*`
  - 仅消费稳定 solver 结果，不拼装 solver 内部 kwargs。

### Batch-C：防回耦门禁

- [x] C1 新增契约稳定性测试（API 层）
  - 校验关键字段存在性、类型、单位语义。
- [x] C2 新增数值回归点（FixedMu/FixedRho/FixedAsymmetricRho）
  - 防止架构重构引入系统性数值漂移。
- [x] C3 新增依赖约束检查
  - 禁止 `phase/scans/workflows` 直接引用 solver 私有符号（如下划线函数、内部字段）。
- [x] C4 文档治理联动
  - `docs/dev/active`、`docs/api` 与导出接口保持一致。

### Batch-D：兼容收口与清理

- [x] D1 标记并移除到期过渡字段（按兼容策略逐步）。
- [x] D2 清理不再使用的 legacy kwargs 路径。
- [x] D3 更新迁移说明，明确“何时删除、替代方式、示例代码”。

## 5. 逐函数/模块改动清单（函数名 + 预期 diff 类型）

### 5.1 `src/models/solver/Solver.jl`

- [x] `SolverResult` 与 `_coerce_solver_result(...)`
  - diff 类型：`契约固化`
  - 加入版本化/稳定字段约束（必要时补适配层）。

- [x] `solve_constraint(...)` 系列入口
  - diff 类型：`兼容收口`
  - 移除过渡参数读取，统一走稳定参数集合。

### 5.2 `src/models/solver/SolverDiagnostics.jl` + `SolverDiagnosticsTypes.jl`

- [x] 诊断对象构造函数
  - diff 类型：`公共/私有字段分层`
  - 公共字段稳定，内部字段通过 debug 扩展通道暴露。

### 5.3 `src/models/phase/PMPhaseDiagnostic.jl`

- [x] `_pm_extract_solver_diagnostic(...)`
  - diff 类型：`适配层改造`
  - 优先读取稳定契约，保留老字段兜底（兼容窗口期）。

### 5.4 `src/models/scans/TmuScan.jl` / `TrhoScan.jl`

- [x] `_to_solver_result(...)`, `_is_success(...)` 等结果消费函数
  - diff 类型：`契约对齐`
  - 统一以稳定字段判定，不依赖临时/内部字段。

## 6. 验收标准（DoD）

- [x] solver 对外存在明确、文档化、可测试的稳定契约。
- [x] phase/scans/workflows 不再直接依赖 solver 私有实现细节。
- [x] 新增依赖门禁可拦截“反向耦合”与“私有符号穿透”。
- [x] 数值回归在既定容差内通过。

## 7. 验证计划（门禁级）

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [x] `julia --project=. scripts/dev/check_solver_contract_leakage.jl`
- [x] （新增）契约稳定性测试入口（本任务落地时补充到 `tests/` 分层目录）。

## 8. 风险与回滚点

- [ ] 风险 A：清理兼容字段导致历史调用断裂。
  - 回滚点：恢复适配层，不回滚核心契约分层。

- [ ] 风险 B：依赖门禁过严导致正常开发受阻。
  - 回滚点：先以 warning 模式运行，再切 hard fail。

- [ ] 风险 C：契约版本推进与文档不同步。
  - 回滚点：冻结版本发布，先修文档与测试。

## 9. 建议提交粒度

- [ ] Commit 1：契约文档 + 类型导出整理。
- [ ] Commit 2：phase/scans/workflows 适配层改造。
- [ ] Commit 3：契约稳定性测试 + 依赖门禁。
- [x] Commit 4：兼容字段清理与迁移说明。

## 10. 与前两阶段的衔接检查

- [x] Phase 1：主链已唯一化（无旁路分流）。
- [x] Phase 2：后处理/诊断/配置拆分已完成。
- [x] 满足以上前置后再执行 Phase 3 清理动作。
- [x] 三份 active 文档中的验证命令集保持一致（含 docs governance 检查）。

## 11. 执行记录

- [x] 2026-04-08：完成 Phase 3 任务单草案，作为三阶段计划收口文档。
- [x] 2026-04-08：完成与 Phase 1/2 文档的联动校对（前置、命令、衔接项对齐）。
- [x] 2026-04-09：完成 SolverResult / Diagnostic 稳定契约 v1 固化（版本字段、公共视图、兼容 coerce）。
- [x] 2026-04-09：完成 phase/scans 契约消费适配（PM 诊断公共视图、Tmu/Trho 结果 coerce）。
- [x] 2026-04-09：新增门禁脚本 `scripts/dev/check_solver_contract_leakage.jl` 与契约回归测试入口。
- [x] 2026-04-09：完成 Batch-D 收口（移除 `fixedmu_problem_spec_active`、`fixedmu_use_problem_spec`、`legacy_fallback_plugin` 与 legacy fallback 透传路径）。
- [x] 2026-04-09：新增 `tests/regression/models/test_solver_phase3_fixedpoint_regression.jl` 覆盖 FixedMu/FixedRho/FixedAsymmetricRho 回归点。
- [x] 2026-04-09：补齐迁移窗口期与门禁上线日期说明（见 `docs/api/models/solver/ResultDiagnosticErrorContracts.md`）。
- [x] 契约版本号：`SolverResult.contract_version=:v1`，`SolverDiagnosticSummary.diagnostic_version=:v1`。
- [x] 迁移窗口期：2026-04-09 至 2026-04-30。
- [x] 门禁上线日期：2026-04-09（`check_solver_contract_leakage.jl`）。
