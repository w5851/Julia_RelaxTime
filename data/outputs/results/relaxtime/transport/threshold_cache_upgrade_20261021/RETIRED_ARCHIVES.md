# 退役归档总索引

本目录下的源码快照 zip 已按 [数据归档策略](../../../../../../config/governance/data_archive_policy.toml) 登记。
策略文件是唯一事实源；本文件只做检索入口与还原指引。

门禁：`scripts/dev/check_data_output_path_guard.jl`（CI: `data-output-path-guard`）——
新增/改名 zip 必须登记，删除必须留墓碑，已登记条目的 sha256 必须与工作区一致。

## 已退役

| 归档 | 条目 | 快照时间 | 还原锚点 |
| --- | --- | --- | --- |
| `full_rs_runner_pilot_v1/calculation-source.zip` | 340 | 2026-10-26 | commit `2e5c411c716d5cbc524c3fc6cfb335fe05b39444`（tag `transport-full-rs-20261026-v1`）；blob `75a52f6884f2dc70b2959896cfb373eece6f6821` |
| `full_rs_a_reduction_probe_v1/source.zip` | 1598 | 2026-10-25 | blob `caa0218765dbe11b3e6d2f9fca393a22a0d689b9`（声明 commit 仅 head_base） |
| `full_rs_a_reduction_probe_v2/source.zip` | 1598 | 2026-10-25 | blob `1537c2c11f0f641c2d037f5201647b094e2ec1df` |
| `full_rs_baseline_acceptance_setup_v1/source.zip` | 1631 | 2026-10-25 | blob `24e0475f90fd4a1b2078d48a0b4dc3dabae6d29c` |
| `full_rs_layer_validation_setup_v1/source.zip` | 1601 | 2026-10-25 | blob `a65786ce84e36935465f5f707253f23683c45535` |
| `full_rs_layer_validation_setup_v2/source.zip` | 1602 | 2026-10-25 | blob `f742ee09bd0e80acaf88c140ba33dfcd182dfbbc` |
| `full_rs_mode_b_qualification_setup_v1/source.zip` | 1636 | 2026-10-25 | blob `fa26925b1b76fa9ef7f5b5b47d248f53bdf89b33` |
| `full_rs_transport_validation_setup_v1/source.zip` | 1600 | 2026-10-25 | blob `f49ea7216d7e5a817b8e367035b413df730ba392` |
| `full_rs_runner_review_checks_v1/source.zip` | 2 | 2026-10-26 | blob `243798c10a956ba8ac793faa08568371d6b305ab`；两个文件也可由 commit `2e5c411c716d` 精确重建 |

### 还原命令

在仓库根执行，任选其一：

```bash
# 逐字节还原原始 zip（把 BLOB 与 DEST 换成对应条目的值）
git cat-file blob BLOB > DEST

# runner pilot 的等价重建（该包与 commit 逐字节一致）
git archive --format=zip 2e5c411c716d5cbc524c3fc6cfb335fe05b39444 \
  > data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/calculation-source.zip
```

## 快照性质

8 个 `source.zip` 是 2026-10-25 / 2026-10-26 的**脏工作区快照**：相对各自声明 commit
`dffe4070eb5c`，每个包有 45 个同名文件内容不同、67–81 个文件尚未入库，因此
`dffe4070eb5c` 只是 `head_base`，切回该 commit 无法重建打包态。

按 CRLF 归一化逐 blob 核对，这批快照里有 **13 处内容从未进入任何提交**
（清单见策略文件 `[unique_snapshot_blobs]`），例如 `ThresholdUpgradeCommon.jl` 的 7.0 KB 版本、
`relaxtime-phase-guided-transport-production.yml`、`tests/unit/runtests.jl` 的当时版本；
另 12 处候选差异经归一化确认只是换行符差异，可由历史提交还原。

## 保留

| 归档 | 条目 | 状态 |
| --- | --- | --- |
| `full_rs_eta_augmented_20261029_v1/eta_augmented_results_v1.zip` | 67 | `historical_superseded`——v15 正式图层（76 份 PDF）的数值源，已被 `ward_moller_sphere_20261121_v1` 取代；保留为历史数值记录，不作现行数据源 |

## 消费方说明

`calculation-source.zip` 与各 `source.zip` 删除后，以下两处仍会引用这些文件名，但它们只服务
2026-10 批次的重跑，且该批次在 `config/governance/task_tracks.toml` 中已归档：

- `scripts/analysis/relaxtime/ThresholdUpgradeCommon.jl`：从 `package_manifest.json` 所在目录推导
  `source.zip`，缺失时报 `ArgumentError("Missing frozen source archive")`；
- `scripts/analysis/relaxtime/run_transport_batch_queue.py`：构建 dispatch 包时复制
  `pilot/calculation-source.zip`。

需要重跑时先按上面的还原命令取回文件，或改用对应的 Ward/Møller 现行流程。
