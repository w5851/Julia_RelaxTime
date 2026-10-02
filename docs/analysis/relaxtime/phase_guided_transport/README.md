# Phase-guided transport evidence line

本目录收纳 phase-guided relaxation-time transport 的连续诊断证据包。它们都位于 `docs/analysis` 的 diagnostic 边界内，不是新的 production result root，也不改变 `data/outputs` 中的正式 CSV、figure 或 registry。

## Stages

| 阶段 | 当前路径 | 角色 |
| --- | --- | --- |
| v1 tau-first analysis | [`phase_guided_transport_p128_xi001_analysis/`](phase_guided_transport_p128_xi001_analysis/) | 从 tau/channel-rate 突变出发，记录下游输运响应、denominator-chain 候选和 claim ledger |
| v2 pole-sensitive rendering | [`phase_guided_transport_v2_pole_sensitive_rendering/`](phase_guided_transport_v2_pole_sensitive_rendering/) | 在 v2 on-shell-kernel production 上迁移 v1 机制窗口，加入 v2 定点诊断、pole-sensitive mask、论文候选显示和一阶分支保护 |
| publication-clean v6-v9 history | [`publication_review_history_v6_v9_archive_v1.json`](publication_review_history_v6_v9_archive_v1.json) | 已核验可恢复的本地历史归档；完整中间图包不进入当前 Git 主线，生成器依赖保留；不宣称历史公共合同仍有效 |
| publication-clean v10 PNG review candidate | [`phase_guided_transport_publication_clean_v10_png_review/`](phase_guided_transport_publication_clean_v10_png_review/) | 基于同一 v5 点表；参数/端点 key 分放 (a)/(c)，刻度精度与实际位置一致，统一收紧复合图字体/留白，新增曲线和端点遮挡 gate；保留 PNG-only 字形例外，未正式化 |
| publication-clean v11 accepted stage | [`publication_clean_v11_stage_acceptance_v1.json`](publication_clean_v11_stage_acceptance_v1.json) | 作者接受的阶段性显示结果；绑定 74 张 PNG、74 张矢量 PDF 和尺寸报告，不授予投稿或新数值生产资格 |
| publication-clean v11 PNG snapshot | [`phase_guided_transport_publication_clean_v11_png_review/`](phase_guided_transport_publication_clean_v11_png_review/) | 不可变 PNG 审查快照及 v5 点表继承链 |
| publication-clean v11 PDF companions | [`phase_guided_transport_publication_clean_v11_pdf_review/`](phase_guided_transport_publication_clean_v11_pdf_review/) | 相同冻结 renderer 的真正矢量补充，保留复合图字高 preflight 未通过项 |
| publication-clean v11 paper-size review | [`phase_guided_transport_publication_clean_v11_paper_size_review/`](phase_guided_transport_publication_clean_v11_paper_size_review/) | Letter 插入和 A4 打印缩放实测；论文临时 PDF 不入计算仓库 |
| publication-clean current | [`phase_guided_transport_publication_clean_v5/`](phase_guided_transport_publication_clean_v5/) | 正式文稿采用的当前 publication figure layer；只更新图内单位和分支端点措辞，不改变 raw/production 数值 |
| publication-clean accepted parent | [`phase_guided_transport_publication_clean_v4/`](phase_guided_transport_publication_clean_v4/) | v5 的已接受显示父级；保留 v4 的局部显示调整与 raw/display provenance |
| publication-clean intermediate | [`phase_guided_transport_publication_clean_v3/`](phase_guided_transport_publication_clean_v3/) | v4 的不可变父级和中间版本，仅保留 provenance，不作为当前 publication layer |

v2 的 `tables/window_classification.csv` 直接引用 v1 的 `tables/mechanism_window_summary.csv`；v2 是连续审计阶段，不覆盖、不合并、也不重写 v1。

## Reading order

1. 先读 v1 `README.md` 和 `manifest.json`，确认 tau-first 机制窗口及其诊断边界。
2. 再读 v2 `README.md` 和 `manifest.json`，确认 v1 -> v2 transfer gate、pole-sensitive display 规则和一阶保护。
3. 需要核对具体数字时，沿各包的 `tables/claim_ledger.csv`、机制表和 figure manifest 回溯。

## Migration boundary

- 本次整理只改变物理 namespace、分组入口和 live 脚本默认路径。
- 两个包内的 CSV、JSON、PNG、生成时 manifest、图 manifest、旧路径快照、execution/provenance 记录均按字节保持。
- v1 root `manifest.json` 在迁移前已经存在 2/21 个 `outputs` hash mismatch；这不是本次目录移动造成的变化，也不在本批修复。
- metadata 修复已登记为独立 follow-up：[`2026-08-19_docs-analysis-metadata-repair-task.md`](../../../dev/active/2026-08-19_docs-analysis-metadata-repair-task.md)。

## Publication promotion boundary

v5 已由作者明确用于正式文稿，并成为当前 publication-clean **图层**。该接受只覆盖
`author_accepted_formal_layout` 的显示产物，不晋升 phase reference，不修改
solver、transport kernel 或 production registry，也不把显示派生值解释为新的
production-grade 数值结果或 convergence certificate。v5 的
`manuscript_eligible=true` 只适用于 v5 display layer；raw numerical status 仍为
`diagnostic_only`，`raw_manuscript_eligible=false`，且 local high-rate gate 未运行。
v4 及其资格记录保持不变，作为 v5 的父级 provenance。

以下 v6-v9 描述是历史审查沿革；产物已归档，不是当前运行合同或推荐入口。

`publication_clean_v6` 是基于 v5 的历史审查候选，不替换 current 指针，也不进入
正式文稿。它重新生成 72 张图（mode A/B 各 36 张），启用四侧 inward 主/小刻度，
在线性轴使用 `AutoMinorLocator(5)`、在 log 轴使用 `LogLocator` 小刻度，并将源单图
图例设为 20 pt；一阶端点图例改为 `1st-order transition (chirally restored endpoint)` /
`1st-order transition (chirally broken endpoint)`，源图在括号前换行。v6 不调用
solver 或 high-rate convergence gate，`manuscript_eligible=false`；v5 的点表、数值表、
平滑记录、raw provenance 和 `publication_clean_current.json` 均保持不变。

`publication_clean_v7_png_review` 使用同一 v5 冻结点表，遵循
`figure_production_v2`/`candidate_aps_v2`：固定物理尺寸的 600 dpi PNG、圆括号单位、
统一最终字号、四侧 inward 主/小刻度、公共图外图例和颜色+线型区分；额外生成
Figure 1/2 的原生 mode-A 复合 PNG。温度映射与完整端点解释记录在
`caption_handoff.md`。它是 `delivery_stage=png_review`、
`vector_delivery_pending=true`、`manuscript_eligible=false` 的审查层，保持独立 y 轴
范围和已审核相变断线，不新增平滑或任何数值，不修改 v5/current 或 raw provenance。

`publication_clean_v8_png_review` 仍消费同一 v5 冻结点表，不修改当前指针。
它是经 case-level 审查的布局例外：单图使用 `best_in_axes_reviewed`，两张
复合图按曲线可见性选择 designated host panel，使用 `best_in_axes_reviewed`；端点 marker 与实际
数据 series 颜色匹配。线性 panel 使用 5 个 x 主刻度和约 5 个 y 主刻度，曲线
保留 6% y 边距；所有正值且 `data_max/data_min >= 20` 的 panel 使用局部 log-y，
具体 panel 和比值记录在 `caption_handoff.md` 与各图 manifest。v8 仍为
`manuscript_eligible=false`、PNG-only 审查层，等待人工审核后再补矢量交付。

`publication_clean_v9_png_review` 继续消费同一 v5 冻结点表，使用当前 v2 SOP 的
线性小刻度密度和紧凑复合图布局；Figure 1 的 `muB900.0` 完整列强制使用 log-y，
图例明确说明开放圆/方 marker 的一阶相变端点含义。v9 仍为
`manuscript_eligible=false`、PNG-only 审查层，不修改当前指针或 raw provenance。

`publication_clean_v10_png_review` 依照组内博士论文示例采用 panel 内的两类
key：参数曲线在 (a)，带 `1st-order` 标题的 restored/broken marker 在 (c)，
各出现一次。子图编号全图统一位于框外左上方，避免压到图例或曲线。
复合图轴标签、标题和数字刻度为 11 pt，图例为 10.5 pt；普通主文本字高
约 2.5 mm，但数学上下标约 1.69 mm。因此它使用显式的 PNG-review-only
字形例外，不宣称满足完整 APS 投稿合同。曲线及端点 marker 遮挡检查为
独立 gate，且三位小数的 `sigma/T` 刻度使用精确 integer-milli grid。
v10 保持 `manuscript_eligible=false`、`current_publication_layer=false`，不改 v5/current；
72 张单图和两张复合 PNG 的 manifest/hash 与完整温度映射用于作者审核。

`publication_clean_v11_png_review` 保留 v10 的物理尺寸和布局，以最新 SOP
统一 Figure 1（四行弛豫时间图）的 12 个 panel 为 log-y，使用普通数字
log 标签；匹配的 mode-A tau 单图也采用 log-y。Figure 2（三行输运系数图）
仍为线性，各行小数位统一为 eta/s 两位、zeta/s 两位、sigma/T 三位。
`First-order` 标题和图例条目左对齐、同字号，caption 明确参数 key 全图通用、
端点 key 仅对应 muB=900 MeV / alpha_T=1.0 的一阶端点曲线。独立 y 范围
保留，普通数字标签不改变 log 位置。v11 是独立 PNG-only 审查 sibling，
完整记录 v10 保留文件 hash，不改 v5/current、raw provenance 或任何数值；
字形审查例外和 `manuscript_eligible=false` 边界不变。

## v11 阶段性验收与 Git 保留

2026-10-02 作者明确接受 v11 为阶段性结果。验收采用追加的
`publication_clean_v11_stage_acceptance_v1.json`，不重写上述审查 manifest。
验收只覆盖显示版式与同 renderer 的 PDF 补充；上下标字高、实际插入宽度
和 numerical qualification 的未通过/未运行边界继续保留。
`manuscript_eligible=false`，`publication_clean_current.json` 仍为 v5。

当前 Git 保留完整 v10 父级和 v11 PNG/PDF/版面审计，以及 v6-v11 生成器
导入依赖。v10 原 SOP/skill hash 与当前合同存在两处已知漂移；历史快照
检查要求 v11 锁定的 v10 文件字节完全不变，并精确记录该漂移，不把它当作
当前合同通过，也不改旧 hash。v11 的当前输入合同必须全部通过。

v6-v9 中间图包、旧 artifact tests 和 v6 镜像脚本已做 ZIP 逐文件 SHA-256
归档及恢复核验，归档索引保留在本目录。当前归档只在本机可用，
`remote_archive_available=false`；不能声称远端 checkout 能自动恢复历史图包。
恢复必须写入新的空目录，不能覆盖现有工作区。

```powershell
python scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py --check
python scripts/analysis/relaxtime/archive_phase_guided_publication_review_history.py restore --archive <zip> --destination <new-empty-directory>
```

版面报告中的论文路径是历史测量来源，不是 clean checkout 的运行依赖。
便携测试核对仓库内图形 hash 和尺寸公式；有原论文临时文件时，可显式设置
`JRT_LOCAL_PAPER_PREVIEW_CHECK=1` 再复核外部输入及 A4 矢量预览。
