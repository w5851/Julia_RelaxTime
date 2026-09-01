# Analysis Artifact Index

`docs/analysis` 保存可追溯的诊断分析、历史比较、研究过程和证据包。这里的文件不是统一的 production result root；每个带 `manifest.json`、`decision.json` 或 `AUDIT.md` 的目录都应按自身 provenance 和 verdict 解读。

本索引先建立逻辑分组；已完成迁移的条目使用当前 canonical path，尚未迁移的条目仍保留原路径。脚本、任务单、manifest 和 figure registry 的路径更新必须与对应的物理迁移在同一批次完成。

## Status Vocabulary

- `diagnostic_only`：可用于观察、比较或作者复核，不是 production/reference 晋升结论。
- `candidate`：有明确的后续验证或作者决策入口，不能直接当作最终物理结论。
- `historical`：保留旧比较、旧图或旧实现上下文，不代表当前实现语义。
- `governance`：记录资产注册、清理、迁移或分析流程，不是物理结果。
- `archive_pointer`：指向外部或 Actions/local artifact 的 provenance 入口，本身不是完整原始数据包。

## Logical Map

| 逻辑域 | 当前物理路径 | 角色 | 关系与入口 |
| --- | --- | --- | --- |
| PNJL / Issue #130 | [`pnjl/`](pnjl/README.md) | C1/C2、CEP、Maxwell、phase-reference 证据线 | 按 C1、C2、CEP/Maxwell、phase-reference 分组；各版本目录保持独立 |
| PNJL / C1 surface views | `pnjl/c1_surface_views/` | 同一 C1 source run 的两个坐标投影 | 共享输入但图轴语义不同；保留为两个独立 view |
| PNJL / Mott | [`pnjl/mott/`](pnjl/mott/) | 独立的 Mott/复极点和文献解释链 | 归入 PNJL 域，但不属于 Issue #130 phase-reference 线 |
| Relaxation-time / transport | [`relaxtime/`](relaxtime/README.md) | phase-guided transport、T200 tau spike、传播子机制诊断 | v2 明确继承 v1；T200 已收纳为独立 diagnostic package |
| Historical comparison | `relaxtime/historical/`、`historical/legacy/` | 历史图、文献比较和复用审计 | 保留上下文；不能自动升级为 strict、regression truth 或 external validation gate |
| Figure governance | `governance/figure_asset_registry_v1/` | 资产清理、迁移和 provenance 快照 | 这是治理材料，不与科学分析包合并 |
| Analysis protocol | `governance/literature_to_implementation_protocol.md` | 文献到实现的工作协议 | 流程文档，不是具体 case 的结果包 |

## PNJL Timeline

### C1 and C2 surfaces

- C1 的两个顶层目录是同一 source run 的不同 projection view；`mu_xi_T` 和 `xi_T_mu` 的图轴、CEP 展示和 manifest schema 不同。
- C2 `pnjl/c2_surface_views/` 下的 `c2_phase_surfaces_diagnostic_v1` 到 `v5_no_triangulation` 是连续的语义/显示演进，不是五份可任意去重的副本：v1/v2 修正化学势口径，v3 增加 crossover 物理筛选，v4 处理视觉闭合，v5 保留 native-support/no-triangulation 语义。
- Issue #130 C2 surface 与 phase-reference figure layer 的统一只读系列索引为 `pnjl/phase_surface_series/`；该 namespace 只保存逐字节相等的 snapshot，不移动或改写原始 evidence package。
- `v4_visual_closed_display16` 是 v4 的展示变体。它与 v4 共享主图和多张表，但拥有独立的 `README.md`、`AUDIT.md`、`decision.json` 和 manifest，因此暂保留为独立 provenance 节点。
- `pnjl/c2_audits/` 收纳 C0/C1/C2 convergence audit 与当前 C1/C2 blocking audit；两者保留独立的 run、verdict 和证据表。

### CEP, Maxwell and phase reference

这些目录形成一条证据链，但每个阶段回答的问题不同：

1. `pnjl/cep_maxwell/narrow_pilot/cep_narrow_pilot_v1/`、`pnjl/cep_maxwell/narrow_pilot/cep_narrow_pilot_v2/`：局部 CEP 三态合同和 pilot 演进。
2. `pnjl/cep_maxwell/stagec/`：四个独立的 Stage-C feasibility/certificate/tolerance replay 包。
3. `pnjl/cep_maxwell/cascade_shadow/`：cascade production shadow 和较早的 endpoint-local full shadow；`pnjl/cep_maxwell/endpoint_local/` 保留后续 endpoint-local contract/shadow 版本。
4. `pnjl_maxwell_*`：Maxwell endpoint、limit contract 和 tolerance feasibility。
5. `pnjl/phase_reference/`：对上述证据的冻结、limited-evidence 和 manual-overlay 汇总；它们是 decision/audit 层，不替代输入证据包。

### Supporting and follow-up artifacts

`pnjl/c2_followups/` 收纳 `c2_limited_feasibility_v1`、`c2_cep_xi05_high_side_extension_v1`、`c2_manual_bisection_v1` 和 `c2_targeted_manual_review_v1`；它们分别是输入合同、补点计划和人工复核包，不应与完整 phase-surface case 混成一个结果目录。`raw_curve_archive_v1` 仍是独立的外部归档指针。

`pnjl/algorithmic_feasibility/` 收纳独立的解析算法可行性审计；它不替代 PNJL 数值验证，也不进入 phase-reference promotion。

### Mott and complex-pole evidence

`pnjl/mott/` 收纳各向异性参数 `xi`、Mott 温度、介子质量/宽度、复极点机制和文献定位的独立分析线。该逻辑组与 Issue #130 的 phase-reference 决策层保持分离；其 `figures/` 和说明文档用于诊断与论文讨论准备，不是新的 production result root。

## Relaxation-Time Timeline

1. `relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis/` 是 v1 tau-first 突变和 denominator-chain 分析主包。
2. `relaxtime/phase_guided_transport/phase_guided_transport_v2_pole_sensitive_rendering/` 是基于 v2 on-shell-kernel production 的迁移和派生显示审计；其表格直接引用 v1 机制摘要，同时新增 v2 定点诊断和一阶分支保护。
3. `relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1/` 是从作者接受的 `prod_v2` raw 派生的论文显示候选层；它保留 raw/clean 双值、19 个相邻点插值和当前 direct-coexistence 两侧 marker，`manuscript_eligible=false`。为便于浏览，18 张已审核 PNG 的字节保持镜像位于 [`data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/`](../../data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/)，对应的图层 manifest 记录在该目录内；分析包仍是派生表和完整 provenance 的 canonical source。
4. `relaxtime/transport/t200_tauu_spikes/` 现在收纳 T200 tau-u 双窗口机制包，包括两份 Markdown 说明和 `tauu_pos_*.png` 机制图；其默认绘图输出已同步指向该目录。

## Organization Rules

- 先用本索引表达 `line`, `role`, `status`, `depends_on` 和 `supersedes`，再考虑目录迁移。
- 版本化 evidence package 保持自包含；同 hash 只说明部分派生表或图相同，不自动说明整个 case 可以删除。
- 汇总包和输入包分开：`pnjl/phase_reference/` 下的 phase-reference 包可以指向下游证据，但不能吞并或改写其 manifest/hash。
- 历史图和治理快照只做逻辑归类，不批量重命名或重绘。
- 新增分析包优先使用 `<domain>/<case>_analysis/README.md + manifest.json + figures/ + tables/` 结构，并在本索引登记。

## Safe Follow-up Order

1. 维护 `pnjl/c1_surface_views/` 两个 projection view 的路径、manifest 和轴语义边界。
2. 继续维护 C2、CEP/Maxwell、phase-reference、Mott 和 phase-guided transport 各自的逻辑线边界，不合并不同 verdict 或 provenance 的证据包。
3. 每批物理迁移都必须同步当前索引、任务单和必要的 live 入口，并通过 docs consistency；历史快照和 metadata repair 保持独立。
