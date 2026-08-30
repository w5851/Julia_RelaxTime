# Literature To Implementation Protocol

Status: project routing protocol.

`Julia_RelaxTime` uses literature review to support computational decisions: model formulas, numerical methods, parameter choices, validation targets, regression baselines, documentation, and paper handoff notes.

The main citation and BibTeX library remains in `D:\Desktop\paper\bib`.

## Project Agents

- `relax-literature-search-strategist`: targeted searches for formulas, parameterizations, algorithms, validation points, and reproducibility signals.
- `relax-method-reviewer`: method assumptions, equation correctness, units, reproducibility, implementation risk, and validation needs.
- `relax-evidence-synthesizer`: method comparison, engineering decision paths, validation plans, and documentation implications.
- `relax-gap-analyst`: gaps between literature and implementation, turned into experiments, regression tests, docs tasks, or paper handoffs.

## Boundary With paper

- Do not edit `D:\Desktop\paper\bib` from this project.
- Do not create a competing master bibliography in this repository.
- When citation metadata, duplicate cleanup, or citekey normalization is needed, hand off to the `paper-citation-curator` in `D:\Desktop\paper`.
- Keep local literature notes focused on implementation relevance, not complete bibliography governance.

## Evidence-To-Code Rules

- Literature evidence alone is not enough to change source code, tests, or baselines.
- Every implementation recommendation should state assumptions, equations, units, parameter definitions, validation targets, and regression risk.
- Use repository test layering: unit, integration, regression, validation, and benchmark as appropriate.
- Do not loosen tolerances to match literature without explaining numerical and physical meaning.
- Treat archived development docs as historical evidence unless explicitly requested for provenance.

## Formula-route closure gate

向新领域拓展、增加新可观测量、引入外部作用量/传播子/响应函数，或把诊断路线
提议为 production candidate 时，必须先建立一个公式路线闭合包。该门禁由
`docs/guides/sop/formula_route_closure.md` 执行，机读登记位于
`config/governance/formula_route_closure.toml`，检查命令为
`julia --project=. scripts/dev/check_formula_route_closure.jl`。

### 必需内容

闭合包至少需要把下面九项放在同一条可审阅链上：

1. **问题和边界**：目标 observable、模型区间、明确的非目标和 diagnostic/formal
   边界；
2. **规范记号**：所有符号、粒子/反粒子约定、单位、化学势和简并度；
3. **模型起点与推导**：从微观拉氏量/哈密顿量开始，逐步写出平均场、二次作用量、
   响应/传播子、统计权重和最终输出；
4. **外部来源矩阵**：每条外部公式都记录 DOI/arXiv、方程号或页码、版本，以及
   与项目符号、归一化、单位和正则化的转换；
5. **选择与替代项**：列出相互冲突的文献约定、项目选择和选择理由，禁止在代码
   中用未登记默认值静默解决冲突；
6. **公式 → 代码 → 测试映射**：每个关键公式对应源码入口、配置、测试层和已知
   缺口；
7. **不变量与验收**：维度、对称性、极限、守恒量，以及数值收敛/回归/外部验证
   的具体判据；
8. **生产状态**：`draft`、`candidate`、`production_authorized` 或 `deprecated`，
   并明确谁可以进行人工升格；
9. **冻结审查快照**：基线 SHA、依赖环境、日期、命令和保留的失败/诊断证据。

“外部公式已经发表”不等于“项目路线已闭合”。如果来源和项目的极化函数、
传播子分母、相位分支或自然单位不同，必须先写出代数转换并测试不变量；不能
把不同论文的公式片段直接拼接后再以数值调参掩盖差异。

### 风险等级和状态转换

- **轻量**：只涉及单位/符号或已有公式的无语义文档修订；仍需纯契约检查。
- **标准**：新传播子、响应函数、统计权重或新的诊断 observable；至少完成公式
  包、unit/integration 测试和一个极限或外部固定点。
- **高风险**：改变平均场驻点、守恒荷、正则化、生产 baseline 或引入新反馈；除
  上述内容外，必须有 regression/validation 和独立人工审阅。

默认状态流转是 `draft -> candidate -> production_authorized`；任一公式/来源/
数值 gate 未关闭时保持 `candidate`。检查器只验证包的结构和一致性，不自动授予
`production_authorized`。如果路线被替代，使用 `deprecated` 并保留 provenance，
不要删除历史公式或失败证据。

### charged-RPA/BU 试点

当前试点 [ChargedRPA_BU_ProductionRoute.md](../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)
展示了该门禁的最小实例：固定 BQS quark-only 背景，明确
`pi^\pm -> K12`、`K^\pm -> K45`，登记 Rehberg/Blaschke 来源、`2/4` 归一化选择、
`num_s_quark=1` 的兼容边界、BU 支撑和未决的严格 retarded/极点 gate。它仍保持
`candidate`，不改变 `PNJLCore`、`MesonDensity` 稳定默认语义，也不更新 production
baseline。

## Expected Handoff

For a literature-backed implementation question:

1. Define the code/model question.
2. Search for targeted evidence.
3. Review methods for reproducibility and assumptions.
4. Synthesize options into a recommended implementation path.
5. Define validation and regression coverage.
6. Hand citation cleanup or manuscript positioning back to `D:\Desktop\paper` when needed.
