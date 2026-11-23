# 退役归档总索引

本目录下的源码快照 zip 已按 [数据归档策略](../../../../../../config/governance/data_archive_policy.toml) 登记。
策略文件是唯一事实源；本文件只做检索入口与还原指引。

门禁：`scripts/dev/check_data_output_path_guard.jl`（CI: `data-output-path-guard`）——
新增/改名 zip 必须登记，删除必须留墓碑，已登记条目的 sha256 必须与工作区一致。

## 已退役

| 归档 | 条目 | 快照时间 | 还原锚点 |
| --- | --- | --- | --- |
| `full_rs_runner_pilot_v1/calculation-source.zip` | 340 | 2026-10-26 | commit `2e5c411c716d5cbc524c3fc6cfb335fe05b39444`（tag `transport-full-rs-20261026-v1`）；blob `75a52f6884f2dc70b2959896cfb373eece6f6821` |

还原命令（在仓库根执行，任选其一）：

```bash
# 逐字节还原原始 zip
git cat-file blob 75a52f6884f2dc70b2959896cfb373eece6f6821 \
  > data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/calculation-source.zip

# 或按冻结 commit 重建等价内容
git archive --format=zip 2e5c411c716d5cbc524c3fc6cfb335fe05b39444 \
  > data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/calculation-source.zip
```

## 待退役（Phase 2）

同批次的 7 个 `source.zip` 与 `full_rs_runner_review_checks_v1/source.zip` 已在策略文件的
`[[pending_retirement]]` 段登记，等待 Phase 2 一并转入墓碑。

它们都是 2026-10-25 / 2026-10-26 的**脏工作区快照**：相对各自声明 commit 的
`dffe4070eb5c`，每个包有 45 个同名文件内容不同、67–81 个文件尚未入库。
按 CRLF 归一化核对，9 个源码包合计有 13 处内容从未进入任何提交，例如
`ThresholdUpgradeCommon.jl`（7.0 KB 版）、`relaxtime-phase-guided-transport-production.yml`、
`tests/unit/runtests.jl` 的当时版本。

## 保留

| 归档 | 条目 | 状态 |
| --- | --- | --- |
| `full_rs_eta_augmented_20261029_v1/eta_augmented_results_v1.zip` | 67 | `historical_superseded`——v15 正式图层（76 份 PDF）的数值源，已被 `ward_moller_sphere_20261121_v1` 取代；保留为历史数值记录，不作现行数据源 |

## 消费方说明

`full_rs_runner_pilot_v1/calculation-source.zip` 删除后，以下两处仍会引用该文件名，但它们只服务
2026-10 批次的重跑，且该批次在 `config/governance/task_tracks.toml` 中已归档：

- `scripts/analysis/relaxtime/ThresholdUpgradeCommon.jl`：从 `package_manifest.json` 所在目录推导
  `source.zip`，缺失时报 `ArgumentError("Missing frozen source archive")`；
- `scripts/analysis/relaxtime/run_transport_batch_queue.py`：构建 dispatch 包时复制
  `pilot/calculation-source.zip`。

需要重跑时先按上面的还原命令取回文件，或改用对应的 Ward/Møller 现行流程。
