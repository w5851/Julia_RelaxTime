# 公式路线闭合 SOP

状态：active（治理 SOP；不授予任何具体数值路线的 production 资格）

本 SOP 是“每次向新领域拓展”时的最小公式闭合门禁。它与专题计算 SOP 分工：
本文件规定如何把模型、外部来源、代码和测试连成可审阅的公式路线；专题 SOP
仍负责具体运行、收敛、产物和升格。

## 1. 目的与适用范围

本 SOP 适用于以下任一情况：

1. 新增物理领域、可观测量或生产候选路线；
2. 把外部文献中的作用量、传播子、相移、输运核或统计权重接入项目；
3. 改变近似层次、归一化、正则化、边界条件或输入/输出单位；
4. 将已有 diagnostic 路线提议为 production candidate。

完成后应得到一个可独立审阅的公式路线闭合包，而不是只有一段代码或一个
“与文献相似”的数值结果。当前 charged-RPA/BU 试点的闭合包是
[ChargedRPA_BU_ProductionRoute.md](../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)。

## 2. 非适用范围

本 SOP 不替代：

- 只修复拼写、路径或日志格式且不改变物理语义的改动；
- 只改变扫描点、输出目录或已冻结公式参数的普通运行任务；
- 完整文献检索、BibTeX 去重和论文投稿格式；这些事项仍由
  `D:\Desktop\paper\bib` 与相应 literature/paper 流程负责；
- 数值收敛、回归基线或正式产物本身。它们必须由专题 SOP 单独验收。

## 3. 权威入口

- 机读登记：[formula_route_closure.toml](../../../config/governance/formula_route_closure.toml)
- 公式闭合检查：`scripts/dev/check_formula_route_closure.jl`
- 文档权威注册：[docs_authority_map.toml](../../../config/governance/docs_authority_map.toml)
- 文献到实现协议：[literature_to_implementation_protocol.md](../../analysis/governance/literature_to_implementation_protocol.md)
- 公式文档索引：[relaxtime formula README](../../reference/formula/relaxtime/README.md)
- 专题运行入口索引：[脚本入口清单](../scripts/README.md)

本 SOP 不创建新的稳定数值 CLI，也不把历史 PNJL 兼容文件写成当前入口。公式
路线包应放在 `docs/reference/formula/`，决定放在
`docs/decisions/`，执行记录放在 `docs/dev/active/` 或归档任务目录。

## 4. 物理口径、单位与参数约束

每个路线必须先写出以下最小闭合链：

```text
微观模型/拉氏量
  -> 近似与平均场背景
  -> 二次作用量、核或响应函数
  -> 传播子/极点/相移/统计权重
  -> 可观测量与输出字段
```

同时固定：

- 外部输入单位、项目内部自然单位和输出单位的转换；
- 所有化学势、守恒荷、粒子/反粒子的符号；
- 每个耦合、积分、响应函数、传播子和密度的量纲；
- 参数值域、支撑域和不允许静默接受的非法输入；
- 不变的物理/求解语义，例如“quark-only 背景”不应被写成包含介子反馈的
  热力学平衡。

路线包必须列出维度检查和至少一个对称/极限检查。若外部来源使用不同符号，
不能只改变量名；必须显式给出符号、因子、单位和归一化转换。

## 5. 输入配置及优先级

公式来源和输入的优先级固定为：

1. 当前模型的微观定义和已接受 ADR；
2. 公式路线包中的项目规范；
3. 外部文献的原始公式；
4. 代码中的兼容实现和默认值。

若层级之间冲突，先记录冲突与选择，再改代码。配置文件只能提供已在路线包
中定义的参数；CLI 覆盖必须落盘到 manifest 或等价的有效配置快照。未登记的
外部公式不能通过“临时常数”进入生产入口。

## 6. 环境与版本冻结

闭合包必须记录：

- 分支、基线 commit、最终审阅 commit 和日期；
- Julia 主版本、根项目环境和必要的可选 oracle 环境；
- 外部文献的 DOI/arXiv/版本或本地文件 provenance；
- 代码、测试和公式文档的相对路径。

如果路线来自未合并分支，必须把基线 SHA 写入文档，并注明哪些实现不属于本
路线。PR290 的 charged-RPA/BU 试点以 `origin/main` 的 `bc9b2990` 为初始基线。

## 7. Smoke 预检

提交前至少运行：

```powershell
julia --project=. scripts/dev/check_formula_route_closure.jl
julia --project=. scripts/dev/check_sop_governance.jl
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. -e 'ENV["UNIT_FILES"]="config/test_formula_route_closure.jl"; include("tests/unit/runtests.jl")'
```

这些命令只检查登记、章节、来源、路径和文档一致性；它们不证明积分收敛、
极点正确或实验趋势可复现。

## 8. 收敛性验证

公式闭合与数值收敛必须分开。专题 SOP 应针对实际路线选择：

- 求解器 residual、分支 seed 和 continuation；
- 动量/能量节点、截断、`eta`、解析延拓和正则化；
- 相位分支、阈值、Levinson 边界和 Bose 支撑；
- 观测量对上述参数的变化及误差来源。

在这些 gate 未通过前，路线状态最多为 `candidate`；不得因为文档已经完整就
把低节点或截断结果写成 formal production。

## 9. 正式计算命令

本 SOP 没有通用的正式数值命令。路线闭合通过后，必须转到对应专题 SOP 的
稳定入口；对于 charged-RPA/BU 试点，当前只允许 diagnostic 适配和审阅，不
允许用本 SOP 直接生成 production 数据。

## 10. 输出目录与产物合同

一个闭合包至少包含：

- 公式文档（模型起点、推导、来源矩阵、决策和未决项）；
- 机读 route registry 条目；
- 代码/测试映射；
- 验证命令和结果摘要；
- 若有数值试点，单独的 diagnostic manifest/CSV，不与 formal baseline 混放。

对于解耦的介子数密度路线，registry 还必须显式记录可调用的算法族、最终比较
方案和 Bose 支撑策略。当前 charged-RPA/BU 试点的最小集合是
`stable_particle_limit`、`reduced_strict_bw`、`q_pole_strict_bw` 和
`phase_shift_bu`；四者都必须保持可调用，`comparison_scheme` 可设为
`phase_shift_gbu_reference` 作为比较默认，但这不改变其它入口的存在性。
`bose_domain_policy=normal_phase_gate_x_min_cut_diagnostic` 表示冻结线先做
正常相门禁，`x_min_cut` 只用于文献复现或异常点诊断；门禁失败时必须停止
production 升格并另建凝聚零模路线。

公式文档不复制整篇论文。每一条外部公式应有 DOI/arXiv、方程号或页码、访问
版本和对项目约定的转换说明。原始文献和主 BibTeX 库仍保留在项目外部指定
位置。

## 11. Regression / Validation 验收

根据风险选择测试层，并在路线包中写明选择理由：

- `unit`：纯代数、量纲、对称性、字段和路径契约；
- `integration`：代码模块接线、输入/输出顺序和分支状态；
- `regression`：会改变已有数值语义时的固定点漂移；
- `validation`：外部文献固定点、极限或实验趋势；
- `benchmark`：新计算层的成本和缓存边界。

若只新增治理文档/检查器而不改变数值实现，可不改已有 production baseline，
但仍必须有 route registry 的纯契约测试。

## 12. 失败点、断点续算与重跑

以下任一情况都必须保留失败原因并停止升格：

- 微观起点、符号或单位无法闭合；
- 外部来源无法定位到方程/版本；
- 代码实现与公式包不一致；
- 必要测试缺失，或数值层存在未解释的 fallback；
- `production_authorized` 与文档状态不一致。

修复后应从同一基线重跑对应检查。数值扫描的断点、失败行和诊断产物按专题
SOP 管理；本 SOP 不允许删除失败证据或用另一条近似静默替换。

## 13. Diagnostic 与 Formal Production 的边界

route registry 的状态含义为：

- `draft`：结构不完整，不能指导运行；
- `candidate`：公式链和来源已登记，但至少仍有未决实现/数值 gate；
- `production_authorized`：公式、代码、数值验证和人工审阅均已完成；
- `deprecated`：保留历史 provenance，并指向替代路线。

`production_authorized` 是显式人工决策，不能由检查器自动晋升。当前
charged-RPA/BU 试点保持 `candidate`：charged ladder/ordered-bubble 的 `2/4`
归一化已经在公式层闭合，但严格 charged retarded bubble、phase/Levinson、
Bose 支撑、单电荷 `domega/pi` 测度迁移和各算法数值收敛尚未实现或验收；second-sheet 极点只作为
`q_pole_strict_bw` 的后续 oracle，不阻塞实轴 GBU。

checker 对 `production_authorized` 路线要求 `unresolved_items=[]`，对其它状态
要求至少保留一项未决内容；同时拒绝仓库外 registry、文档或测试路径。对于
`charged_rpa_bu_quark_only`，它还要求四类介子数密度算法和广义 BU 比较默认字段
齐全。它仍只是 provenance/contract 检查器，不替代物理推导、极点求解或数值
收敛验证。

## 14. 后处理与作图

后处理只能消费路线包声明的字段和状态。diagnostic 输出必须保留：

- route id、公式版本、代码 commit、配置快照；
- 单位、通道顺序、简并因子和缺失/失败语义；
- support、phase branch、fallback 和收敛状态。

任何图表若要进入论文或正式结果，还需遵守 figure production SOP；公式闭合
检查不替代图表的 hash、manifest、DPI 或 provenance gate。

## 15. 关联公式、API 和测试

- [charged-RPA/BU 候选公式规范](../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)
- [KMT 平均场到二次核](../../reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md)
- [极化函数](../../reference/formula/relaxtime/polarization/Polarization_极化函数byB0.md)
- [BU 相移公式](../../reference/formula/relaxtime/meson_density/MesonDensity_BU相移公式.md)
- [MesonDensity API](../../api/relaxtime/meson_density/MesonDensity.md)
- [公式路线纯契约测试](../../../tests/unit/config/test_formula_route_closure.jl)
- [公式路线检查器](../../../scripts/dev/check_formula_route_closure.jl)

当前稳定数值入口仍由专题 SOP 指定；本 SOP 不修改 `Models` 公共 API。

## 16. 最后验证记录

- 验证日期：`2026-08-30`
- 初始基线 commit：`bc9b2990bcfe3b8c32d2ec0f00066b52b4cf800b`
- 最终验证 commit：以 PR290 合并前最后一次验证的 branch tip 和 PR validation block 为准；
  本文件不自引用自身 commit hash
- 执行命令：见第 7 节及 route registry
- 结果：治理 SOP 可执行；charged-RPA/BU 路线 `candidate`、diagnostic-only
- 备注：本 SOP 试点不实现 strict-support、凝聚零模、second-sheet 极点求解、
  `Omega_M` 反馈或显式 `mu_I` 路线；这些若要引入，必须创建新的闭合任务。
