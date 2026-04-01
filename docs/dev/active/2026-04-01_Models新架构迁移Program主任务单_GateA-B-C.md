# Models 新架构迁移 Program 主任务单（Gate A/B/C）

## 1. Program 背景与目标

本任务单用于串联两条已立项主线，防止执行过程漂移：

- 主线 A（Plan-B）：dimension-agnostic + schema-driven + 旧架构/compat 移除
- 主线 B（Homomorphic）：7 模型目录/API/测试骨架同构化（接口同构、内核异构）

Program 最终目标：

- [ ] 完成从旧架构到新架构的彻底迁移
- [ ] 清除旧 compat 与回退路径，不保留长期双轨
- [ ] 在 7 个模型上形成统一边界契约 + 可预测工程结构 + 稳定验证门禁

---

## 2. 关联文档（单一事实源）

### 2.1 Spec

- [ ] `docs/superpowers/specs/2026-04-01-models-dimension-agnostic-plan-b-design.md`（Gate A）
- [ ] `docs/superpowers/specs/2026-04-01-models-homomorphic-architecture-design.md`（Gate B/C）

### 2.2 Plan

- [ ] `docs/superpowers/plans/2026-04-01-models-dimension-agnostic-plan-b.md`（A 线实施）
- [ ] `docs/superpowers/plans/2026-04-01-models-homomorphic-architecture-implementation-plan.md`（B/C 线实施）
- [ ] `docs/dev/active/2026-04-01_Models固定维度耦合点基线清单_GateA.md`（A 线前置基线）

执行约束：

- [ ] 若本任务单与上述 spec/plan 冲突，以 spec/plan 为准并先回填本任务单再实施。

---

## 3. 模型覆盖范围（Program 硬门禁）

必须覆盖当前已注册 7 模型：

- [ ] `NJL`
- [ ] `NJL2`
- [ ] `PNJL`
- [ ] `PNJLMagnetic`
- [ ] `RPNJL`
- [ ] `Rotation`
- [ ] `GasLiquid`

---

## 3.5 开工前前置准备（必须先完成）

### 3.5.1 基线与清单冻结

- [ ] 固定维度耦合点清单已冻结（文件+符号+预期迁移去向）。
- [ ] 7 模型当前行为基线已记录（最小可复现实验命令+输出摘要）。
- [ ] 回归阈值已冻结（残差/关键量容忍区间，至少给出 `rtol/atol` 或范围规则）。

### 3.5.2 工程执行隔离

- [ ] 在隔离分支/工作树中执行（避免与其他开发流互相污染）。
- [ ] 数据输出目录清理策略已确定（禁止将临时产物混入提交）。
- [ ] 提交流程策略已确定（按 chunk 小步提交，避免超大混合提交）。

### 3.5.3 外部契约与回滚策略

- [ ] 外部可见接口影响清单已确认（`Models` 导出、CLI 参数、CSV/JSON 字段）。
- [ ] deprecation 文案与移除窗口已定义（避免“静默破坏”）。
- [ ] 回滚策略已定义（触发条件、回滚范围、恢复步骤）。

### 3.5.4 质量门禁映射

- [ ] 本地验证矩阵与 CI required checks 名称已建立映射。
- [ ] “停止线”条件已定义（例如：核心模型回归失败、治理检查失败即停止推进下一 Gate）。

---

## 4. Gate 串行策略（禁止并行越级）

### Gate A：Plan-B 架构收口（前置门）

入口条件：

- [ ] A 线 spec+plan 已冻结并通过评审
- [ ] 执行环境与验证命令已可复现

退出条件（全部满足）：

- [ ] 主路径去除固定 `5/3/8` 切片依赖
- [ ] schema-driven 状态/残差契约稳定
- [ ] 旧架构与 compat 在 A 线定义范围内已移除
- [ ] A 线验证矩阵通过（unit/integration/regression/governance）

### Gate B：同构化实施（依赖 A）

入口条件：

- [ ] Gate A 完成并记录证据
- [ ] 不存在未决的核心契约变更

退出条件（全部满足）：

- [ ] 7 模型最小公共接口契约一致
- [ ] 目录骨架与 workflow adapter 结构同构
- [ ] 测试模板同构（unit/integration/regression）
- [ ] 不以统一之名篡改模型内核物理语义

### Gate C：Program Final 联合收口

入口条件：

- [ ] Gate A + Gate B 均完成

退出条件（全部满足）：

- [ ] 全模型回归通过
- [ ] docs/governance 全绿
- [ ] 无 compat 残留、无旧架构回退入口
- [ ] 本任务单状态与证据完整回填

---

## 5. 漂移防护规则（执行纪律）

- [ ] 禁止跨 Gate 并行改动核心契约（`state schema / result contract / solver core`）。
- [ ] 禁止在 Gate A 未稳定时提前推进目录同构化大改。
- [ ] 禁止以临时 shim 长期替代正式迁移（每个 shim 必须有移除点）。
- [ ] 每个 Gate 完成后立即回填证据（命令、结果、问题与处置）。
- [ ] 若发现范围膨胀，先更新 spec/plan，再更新本任务单，再实施代码。

---

## 6. 验证命令与证据记录

### 6.1 覆盖矩阵基线

- [ ] `julia --project=. -e 'include("src/constants/Constants_PNJL.jl"); include("src/models/Models.jl"); println(Models.registered_model_kinds())'`

证据：

- [ ] 输出已记录（7 模型齐全）

### 6.2 通用验证矩阵（每个 Gate 收口至少执行一次）

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`

证据：

- [ ] Gate A 验证结论已回填
- [ ] Gate B 验证结论已回填
- [ ] Gate C 验证结论已回填

---

## 7. 执行日志（append-only）

> 仅追加，不回写历史。

### 2026-04-01 初始化

- [x] 建立 Program 主任务单，定义 Gate A/B/C 串行治理。
- [x] 绑定两条主线 spec+plan 作为单一事实源。
- [x] 明确最终目标：彻底迁移到新架构并移除旧 compat。

### 2026-04-01 文档复核补充

- [x] 补充“开工前前置准备”章节，新增基线冻结、执行隔离、外部契约、回滚与质量门禁映射要求。
- [x] 新增 Gate A 基线附件：`2026-04-01_Models固定维度耦合点基线清单_GateA.md`。

---

## 8. Program DoD

- [ ] Gate A/B/C 全部门禁通过
- [ ] 7 模型覆盖矩阵全绿
- [ ] 旧架构与 compat 清零
- [ ] 文档、计划、执行状态一致
- [ ] PR checks 与 review 闭环
