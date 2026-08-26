# Issue #130 RS Sharded Production v2 Post-Repair Audit

本目录只记录 `2026-08-26` RS 15+15 shard numerical run 的 solver-free provenance
复核，不包含原始全量曲线，也不改写外部 Actions artifact。作者已完成
candidate/legacy 全部审核；数值结论仍保持 `diagnostic-only`，而已审核的 candidate
图像另行导入正式 figure 路径，二者不混为同一晋升状态。

## 输入与口径

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow head：`22874505877491754eed27519ad8a7b871c82571`
- case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2`
- 选择的唯一 shard：30 个；重复 dispatch 10 个，已取消且不纳入 aggregate
- 外部源目录：`D:\Desktop\Julia_RelaxTime_issue130_artifacts\rs_sharded_production_v2_20260826`
- aggregate：`aggregate_replay_20260826_v4`
- post-repair audit：`post_repair_audit_20260826_v1`
- provenance repair：PR #269，合并 SHA `472e3eb0047b3a9380b27f1935880b0e473ca9b5`

复核命令：

```text
python scripts/analysis/relaxtime/audit_issue130_rs_sharded_provenance.py --source-root D:\Desktop\Julia_RelaxTime_issue130_artifacts\rs_sharded_production_v2_20260826 --output-root D:\Desktop\Julia_RelaxTime_issue130_artifacts\rs_sharded_production_v2_20260826\post_repair_audit_20260826_v1
```

## 结果

审计 verdict 为 `post_repair_audit_pass_diagnostic_only`。

- 30/30 shard 的 scan 行数、finite/converged、重复键、failed points、hard gate 和 timing 均通过。
- mode-A/mode-B 的 direct-coexistence 合同通过：`xi=-0.003` 与 `xi=+0.003` 存在，严格共存的 `xi=0` 不作为唯一输运点。
- aggregate manifest 自校验通过；没有缺失 artifact 或意外 hash mismatch。
- 60 个历史 sidecar hash mismatch 全部属于已知的 `effective_config.json`（30）与
  `failed_points.csv`（30）写入顺序问题。PR #269 修复未来运行的写入顺序，历史文件保持不可变。
- mode-B `alpha_T=NaN` 是固定-T 路由的预期字段语义，不作为非有限数值失败。
- 质量警告仍记录为诊断：mode-A `tau_u_ubar_ratio_high` 230 条，mode-B 240 条。
- candidate/legacy 的差异已由作者审核：mode-A 有 15 个共同键 phase-field mismatch，mode-B 有
  47 个；这些差异不能自动解释为 solver regression，也不能自动晋升为正式 production。
- 审核表进一步拆分了 15/47 个物理 classification mismatch 与大量 reference route metadata
  difference；后者包括 `phase_boundary_xi_used`、`phase_prev/curr`、seed/backend 路由，不应直接
  当作物理分类回归。数值差异极值另列在 `candidate_legacy_numeric_extrema.csv`。

## 作者审核的 muB=900 切片解释

作者接受以下局部输运解释，适用范围是 candidate v2 的 mode-A、`muB=900 MeV`、
`alpha_T=1.0` 审核切片，而不是对所有 `xi`、`muB` 的严格热力学定律：

1. **观察**：在 `xi=0` 附近，`xi=-0.05 ... -0.003` 的点仍沿同一已收敛的
   first-order/quark 侧分支变化；`xi=+0.003` 是为 direct-coexistence 合同保留的另一侧
   输运点，因而在 `xi=0` 附近出现独立的分支跳变。
2. **局部归因**：在分支内部，增大 `xi` 对当前 phase-scaled 路由的响应可粗略看作
   “有效降温”方向的移动。各 `tau` 分量在接近零点前的快速斜率变化，与各向同性基线中
   降低 `T` 时 `tau` 快速增加的熟悉形状同构；它不是由 `xi=0` 的一阶跳变再次产生的
   第二个相变面。
3. **边界**：这里的“有效降温”是局部响应类比，不定义新的 `T_eff(T,xi)` 热力学坐标，
   也不替代 equilibrium solver、Maxwell 或 phase-reference 的坐标。若要将该机制推广到
   其他切片，仍需相同的分支、背景量和输运通道证据。

对应证据为外部 aggregate 的 `mode_a_scan.csv`（其 SHA 记录在本目录
`manifest.json`），以及已审核的 `muB900.0` 逐输运量图。该段落是
`author_accepted_diagnostic_interpretation`，不是新的 production-grade 收敛证明。

## 边界与下一步

本审计不调用 equilibrium solver、不写入 `data/reference` 或 production result、不删除 legacy
fallback，也不完成 RS candidate promotion。作者审核完成后，下一步仍需另建 versioned RS
promotion/import PR；该 PR 必须保留逐键 legacy fallback 和显式 rollback。当前导入的图像仅是
作者接受的正式布局副本，不改变这条数值边界。

详细表格和不可变 hash 位于外部 audit 目录；本目录的 `manifest.json` 保存其指针与摘要哈希。

## Candidate/legacy 图像审阅包

为便于作者审核，另生成了一个只读的同构绘图包：

- 外部目录：`D:\Desktop\Julia_RelaxTime_issue130_artifacts\rs_sharded_production_v2_20260826\candidate_legacy_figures_v2`
- manifest SHA-256：`9a986e1887292309a46a963dfbc76f08cc7fb67d3fa5664f282ee7256853f01b`
- 目录结构：`mode_a_fixed_muB_phase_scaled`、`mode_b_fixed_T_sparse_muB`，各自包含
  `candidate_runtime` 与 `legacy_runtime`；其下保留对应的 `prod_v2`/`prod_v1` case 名、
  `plot_panel=*` 分面和旧图同名 PNG 文件。
- 四棵树各 36 张 PNG，共 144 张；横轴、分面、曲线分组、12 个输运量、PNG 与 600 dpi
  均复用历史入口 `scripts/relaxtime/run_phase_guided_transport_plots.jl`。
- 生成脚本：`scripts/analysis/relaxtime/render_issue130_rs_candidate_legacy_figures.py`。

该包只证明“candidate 与 legacy 可以在旧图视觉合同下逐张对照”。作者审核通过后，candidate
两棵图树已按旧结果同构导入仓库正式 figure 路径：

- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`

每棵树包含 36 张 PNG 和扩展后的 `plot_manifest.json`。图像状态为
`author_accepted_formal_layout`，数值状态仍为 `diagnostic_only`；历史 legacy v1 正式图不被
覆盖，legacy 对照副本继续保留在外部 audit 包中。导入由
`scripts/analysis/relaxtime/formalize_issue130_rs_candidate_figures.py` 完成，逐张复核外部
manifest/hash，`solver_called=false` 且 `production_write=false`。
