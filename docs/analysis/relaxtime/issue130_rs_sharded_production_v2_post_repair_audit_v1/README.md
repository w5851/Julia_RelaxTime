# Issue #130 RS Sharded Production v2 Post-Repair Audit

本目录只记录 `2026-08-26` RS 15+15 shard numerical run 的 solver-free provenance
复核，不包含原始全量曲线，也不改写外部 Actions artifact。

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
- candidate/legacy 的差异仍需作者审核：mode-A 有 15 个共同键 phase-field mismatch，mode-B 有
  47 个；这些差异不能自动解释为 solver regression，也不能自动晋升为正式 production。
- 审核表进一步拆分了 15/47 个物理 classification mismatch 与大量 reference route metadata
  difference；后者包括 `phase_boundary_xi_used`、`phase_prev/curr`、seed/backend 路由，不应直接
  当作物理分类回归。数值差异极值另列在 `candidate_legacy_numeric_extrema.csv`。

## 边界与下一步

本审计不调用 equilibrium solver、不写入 `data/reference` 或 production result、不删除 legacy
fallback，也不完成 RS candidate promotion。下一步是作者审核 candidate/legacy 差异；审核接受后
另建 versioned RS promotion/import PR，保留逐键 legacy fallback 和显式 rollback。

详细表格和不可变 hash 位于外部 audit 目录；本目录的 `manifest.json` 保存其指针与摘要哈希。
