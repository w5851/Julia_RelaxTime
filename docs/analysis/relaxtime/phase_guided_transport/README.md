# Phase-guided transport evidence line

本目录收纳 phase-guided relaxation-time transport 的连续诊断证据包。它们都位于 `docs/analysis` 的 diagnostic 边界内，不是新的 production result root，也不改变 `data/outputs` 中的正式 CSV、figure 或 registry。

## Stages

| 阶段 | 当前路径 | 角色 |
| --- | --- | --- |
| v1 tau-first analysis | [`phase_guided_transport_p128_xi001_analysis/`](phase_guided_transport_p128_xi001_analysis/) | 从 tau/channel-rate 突变出发，记录下游输运响应、denominator-chain 候选和 claim ledger |
| v2 pole-sensitive rendering | [`phase_guided_transport_v2_pole_sensitive_rendering/`](phase_guided_transport_v2_pole_sensitive_rendering/) | 在 v2 on-shell-kernel production 上迁移 v1 机制窗口，加入 v2 定点诊断、pole-sensitive mask、论文候选显示和一阶分支保护 |
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
