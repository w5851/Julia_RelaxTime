# ADR-0005: 科学计算 SOP 与文档权威治理

## 状态

已接受

## 决策结论

1. 科学计算流程采用“公共生命周期 SOP + 专题工作流 SOP”，不建设单一巨型总 SOP。
2. `docs/guides/sop/` 是面向执行与复现的 SOP 权威目录；公式、API、分析证据和任务过程继续保留在各自既有目录。
3. `config/governance/docs_authority_map.toml` 是 SOP 身份、状态、权威范围、入口和复核周期的机读单一来源。
4. 新 SOP 必须同时满足“written / configured / effective / usable”四层完成条件；仅有 Markdown 文件不代表流程已经可用。
5. 首期只覆盖公共科学计算生命周期、PNJL 相结构和 relaxtime transport，不改变任何物理公式、求解语义或数值基线。

## 背景

仓库已经存在稳定脚本白名单、测试分层、baseline 管理、端到端物理链路、运行 manifest 和局部 runbook，但这些能力分散在 `docs/guides/`、`docs/reference/`、`docs/api/`、`docs/dev/` 与脚本 README 中。

当前主要问题不是缺少文档数量，而是：

- 缺少从环境冻结、smoke、收敛性、正式运行到产物升格的统一执行闭环；
- 部分历史 notes 仍可能被误认为当前工作流；
- 文档一致性门禁可以检查已知旧路径，但不能声明某一流程只有一个权威 SOP；
- 科学计算的“诊断结果”与“正式产物”缺少跨工作流统一的升格规则。

这一方向会指导多个后续专题和维护任务，不能只保留在单个 active 任务单中。

## 决策

### 1. 文档分层

- `docs/reference/formula/`：公式、物理假设、单位与方法依据。
- `docs/api/`：稳定入口、参数、返回值与程序合同。
- `docs/guides/sop/`：从零执行、验证、复现和故障恢复。
- `docs/analysis/`：结果证据、claim boundary 与分析包。
- `docs/dev/active|backlog|archived/`：任务过程和实施历史。

SOP 通过链接消费其他层，不复制大段公式、API 全集或任务历史。

### 2. SOP 结构

所有 active SOP 采用同一组必需章节，至少覆盖：

- 目的与适用范围；
- 非适用范围；
- 权威入口；
- 物理口径、单位与参数约束；
- 输入配置及优先级；
- 环境与版本冻结；
- smoke 预检；
- 收敛性验证；
- 正式计算命令；
- 输出目录与产物合同；
- regression / validation 验收；
- 失败、断点续算与重跑；
- diagnostic 与 formal production 边界；
- 后处理与作图；
- 关联公式、API 和测试；
- 最后验证记录。

### 3. 机读治理

TOML 权威映射记录：

- SOP ID 与文件路径；
- `active / draft / deprecated` 状态；
- 唯一的 `authoritative_for` 范围；
- 稳定入口与相关配置；
- 验证命令；
- `last_verified` 与复核周期。

Julia 治理脚本检查路径存在性、权威范围唯一性、必需章节、旧入口泄漏、稳定入口白名单一致性和复核期限。

### 4. 渐进式接入

首期门禁只约束 `docs/guides/sop/` 与权威映射，不要求一次性修复全部历史 Markdown。旧文档按使用风险逐批标记为 supporting、deprecated 或 historical。

### 5. 科学计算完成语义

- **written**：SOP 与索引已经存在。
- **configured**：入口、配置、单位、输出和验证命令已登记。
- **effective**：smoke 或等价最小验证实际可执行并产生约定产物。
- **usable**：只依赖 SOP 即可判断失败、恢复、收敛和产物是否可升格。

正式产物必须先有收敛性证据；单次运行成功或图像看似合理都不足以升格。

## 理由

- 复用现有 `Models`、稳定 CLI、manifest、测试分层和治理脚本，不引入新的运行框架。
- 将易漂移的状态和入口登记在 TOML，便于 Julia 标准库直接解析并接入 CI。
- 保持公式、API、操作步骤和研究结论的职责边界，减少重复维护。
- 以高风险工作流试点，验证模板和门禁后再扩展介子热力学、介子数密度等专题。

## 备选方案

### 方案 A：仅新增一份总 SOP

优点：初始文件少。

缺点：不同工作流的参数、收敛轴和失败语义会快速膨胀，容易再次过时。

未选择原因：不适合多个物理链路长期独立演进。

### 方案 B：只依赖现有 README 和测试

优点：不新增治理机制。

缺点：无法回答从诊断到正式产物的完整流程，也不能声明权威文档唯一性。

未选择原因：不能解决当前文档分散和复现闭环缺失。

### 方案 C：为所有历史文档统一增加复杂 frontmatter

优点：元数据就地可见。

缺点：迁移面过大，并引入 Markdown 元数据重复和解析成本。

未选择原因：首期采用独立 TOML 权威映射更小、更可验证。

## 后果

### 正面影响

- 新成员和未来维护者可沿统一流程完成可靠计算。
- 稳定入口、配置、文档和验证命令可以被 CI 联合检查。
- 历史说明不再与当前主线竞争权威性。
- diagnostic 与 formal production 的边界可在不同物理专题间复用。

### 负面影响

- 每个 active SOP 需要定期复核。
- 新增工作流时需同步登记权威映射和验证命令。
- 首期只能治理新 SOP，旧文档仍需分批迁移。

### 风险与缓解

- 风险：SOP 复制公式和 API，形成新的重复事实源。
  - 缓解：模板要求链接参考层，只在 SOP 中记录运行所需的最小物理口径。
- 风险：为了通过门禁而更新日期，未实际运行验证。
  - 缓解：`last_verified` 必须与明确的验证命令和产物合同绑定。
- 风险：过早把诊断参数写成正式默认值。
  - 缓解：非默认求解/缓存策略必须显式标记 diagnostic-only，保持现有数值语义不变。

## 阶段与复核

1. MVP：公共 SOP、PNJL 相结构 SOP、transport SOP、权威映射和治理门禁。
2. 扩展：介子热力学、介子数密度和文献复现 SOP。
3. 成熟：将 SOP 与正式产物、分析包和论文证据链进一步关联。

在以下任一条件触发时复核本决策：

- 完成三个专题 SOP 后；
- 稳定 CLI 或产物 manifest 契约发生变化；
- 新增第二套相互竞争的文档元数据机制；
- SOP 门禁连续产生无法解释的误报。

## 与短期任务的关系

首期实施记录已归档至 `docs/dev/archived/2026-07-10_科学计算SOP与文档治理建设任务单.md`。后续专题必须引用本 ADR，但具体任务和验证证据留在各自 active 文档中，完成后按仓库规则归档。

## 相关决策

- [ADR-0002](0002-models-solver-contract.md)
- [ADR-0003](0003-pnjl-solver-decoupling-governance-and-ad-implicit-contract.md)
- [ADR-0004](0004-solver-three-layer-contract-and-vector-kernel.md)

## 参考资料

- `docs/guides/scripts/README.md`
- `docs/dev/testing_governance.md`
- `docs/guides/BASELINE_VERSION_MANAGEMENT.md`
- `docs/reference/formula/relaxtime/transport/Transport_EndToEnd_Pipeline.md`

---

**日期**：2026-07-10
**决策确认**：用户已确认总体落地方向
