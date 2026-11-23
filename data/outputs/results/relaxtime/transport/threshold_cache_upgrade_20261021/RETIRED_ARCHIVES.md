# 退役归档总索引

本目录下的源码快照 zip 已按 [数据归档策略](../../../../../../config/governance/data_archive_policy.toml) 登记。
策略文件是唯一事实源；本文件只做检索入口与还原指引。

门禁：`scripts/dev/check_data_output_path_guard.jl`（CI: `data-output-path-guard`）——
新增/改名 zip 必须登记，删除必须留墓碑，已登记条目的 sha256 必须与工作区一致。

## 已退役

下表「删除前 commit」是墓碑里的 `retired_before_commit`：该提交在主线可达链上，
删包内容即在其树中，可用 `git rev-parse <commit>:<path>` 复核并取回。

| 归档 | 条目 | 快照时间 | git_blob | 删除前 commit |
| --- | --- | --- | --- | --- |
| `full_rs_runner_pilot_v1/calculation-source.zip` | 340 | 2026-10-26 | `75a52f6884f2dc70b2959896cfb373eece6f6821` | `0483e5db` |
| `full_rs_a_reduction_probe_v1/source.zip` | 1598 | 2026-10-25 | `caa0218765dbe11b3e6d2f9fca393a22a0d689b9` | `b2021ffa` |
| `full_rs_a_reduction_probe_v2/source.zip` | 1598 | 2026-10-25 | `1537c2c11f0f641c2d037f5201647b094e2ec1df` | `b2021ffa` |
| `full_rs_baseline_acceptance_setup_v1/source.zip` | 1631 | 2026-10-25 | `24e0475f90fd4a1b2078d48a0b4dc3dabae6d29c` | `b2021ffa` |
| `full_rs_layer_validation_setup_v1/source.zip` | 1601 | 2026-10-25 | `a65786ce84e36935465f5f707253f23683c45535` | `b2021ffa` |
| `full_rs_layer_validation_setup_v2/source.zip` | 1602 | 2026-10-25 | `f742ee09bd0e80acaf88c140ba33dfcd182dfbbc` | `b2021ffa` |
| `full_rs_mode_b_qualification_setup_v1/source.zip` | 1636 | 2026-10-25 | `fa26925b1b76fa9ef7f5b5b47d248f53bdf89b33` | `b2021ffa` |
| `full_rs_transport_validation_setup_v1/source.zip` | 1600 | 2026-10-25 | `f49ea7216d7e5a817b8e367035b413df730ba392` | `b2021ffa` |
| `full_rs_runner_review_checks_v1/source.zip` | 2 | 2026-10-26 | `243798c10a956ba8ac793faa08568371d6b305ab` | `b2021ffa` |

`retired_before_commit` 的完整值见策略文件（`0483e5db58f68f65914d499398261682ed12e3de`、
`b2021ffa112dad868e073c96163734ef2937a991`）。

### Batch 1：2026-10 批次返回证据包（2026-11-23 退役）

| 归档 | 条目 | 快照时间 | git_blob |
| --- | --- | --- | --- |
| `engineering_delivery_20261104_v1/validation/remote_test_input.zip` | 170 | 2026-10-17 | `5b8ea47ce83b…` |
| `full_rs_baseline_acceptance_setup_v1/returned_evidence.zip` | 422 | 2026-10-25 | `8a186bae28c2…` |
| `full_rs_baseline_review_v1/old_baseline_evidence.zip` | 377 | 2026-10-23 | `76256147d7b5…` |
| `full_rs_joint_representatives_p48_v1/returned_evidence.zip` | 885 | 2026-10-25 | `e077c36fc718…` |
| `full_rs_layer_validation_setup_v1/failed_evidence.zip` | 332 | 2026-10-25 | `13b9bfae601f…` |
| `full_rs_layer_validation_setup_v2/returned_evidence.zip` | 1131 | 2026-10-25 | `ff19e17e6dee…` |
| `full_rs_mode_b_background_setup_v1/returned_evidence.zip` | 410 | 2026-10-25 | `b70c5f8513ec…` |
| `full_rs_mode_b_qualification_setup_v1/returned_evidence.zip` | 2824 | 2026-10-23 | `28a5c749383d…` |
| `full_rs_phase_validation_v1/returned_evidence.zip` | 410 | 2026-10-25 | `d2a79399bd19…` |
| `full_rs_production_replay_v2/returned_failure.zip` | 372 | 2026-10-25 | `39e75d6d6164…` |
| `full_rs_production_replay_v3/returned_evidence.zip` | 942 | 2026-10-25 | `8d5b66d2b22e…` |
| `full_rs_transport_validation_setup_v1/returned_evidence.zip` | 374 | 2026-10-25 | `63f8623afd08…` |
| `remote_preflight_setup_v1/transport_preflight_source_v1.zip` | 879 | 2026-10-23 | `29830920f552…` |

这 13 个包没有 `declared_commit`（`uncommitted_source = false`），是运行返回证据、失败残留与预检源码包。
删除前已用 `git grep -F` 对**完整路径**与**父目录**在 `.github/ scripts/ src/ tests/ config/` 全量检索，均 0 命中；
同目录的 `package_manifest.json` / `returned_archive_manifest.json` / `returned_failure_manifest.json` 等清单保留。
完整 `git_blob` 与 `sha256` 见策略文件。


### 还原命令

**请在 Git Bash / POSIX shell 下执行**。Windows PowerShell 5.1 的 `>` 重定向会按 UTF-16LE
转码，破坏 zip 二进制；若只能用 PowerShell，请改用
`[IO.File]::WriteAllBytes($dest, (git cat-file blob $blob | Out-String))` 之外的等价二进制方式，
或先 `git cat-file blob <blob> | Out-File -Encoding Byte`。

```bash
# 逐字节还原（把 BLOB 与 DEST 换成表里对应条目的值）
git cat-file blob BLOB > DEST

# 复核墓碑本身：删除前 commit 的树里该路径应指向同一个 blob
git rev-parse 0483e5db58f68f65914d499398261682ed12e3de:data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/calculation-source.zip

# runner pilot 的等价重建（该包与 commit 逐字节一致）
git archive --format=zip 2e5c411c716d5cbc524c3fc6cfb335fe05b39444 \
  > data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/calculation-source.zip
```

浅克隆（`--depth`、`--filter=blob:none`）下这些 blob 不在本地，需先
`git fetch --unshallow` 或抓取对应的 `retired_before_commit`。

**禁止对 `data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/**` 执行历史重写**
（`git filter-repo` / BFG）。墓碑的还原能力完全依赖 `.git` 中这些 blob 可达；
一旦重写历史，13 处孤本与 7 个脏快照将不可恢复。

## 快照性质

**7 个** `dffe4070eb5c` 批次的 `source.zip` 是 2026-10-25 的**脏工作区快照**：相对该 commit，
每个包有 45 个同名文件内容不同、67–81 个文件尚未入库，因此 `dffe4070eb5c` 只是 `head_base`，
切回该 commit 无法重建打包态。

第 8 个 `full_rs_runner_review_checks_v1/source.zip` 不同：它是 2026-10-26 的**干净包**，
`uncommitted_source = false`、`declared_commit = 2e5c411c716d`，包内两个文件的 blob 与该 commit
逐字节一致，可直接用 `git archive` 重建。

按 CRLF 归一化逐 blob 核对，7 个脏快照合计有 **13 处内容从未进入任何提交**
（清单见策略文件 `[unique_snapshot_blobs]`），例如 `ThresholdUpgradeCommon.jl` 的 7.0 KB 版本、
`relaxtime-phase-guided-transport-production.yml`、`tests/unit/runtests.jl` 的当时版本。
这 13 处是 7 个包内文件的并集，取回时先按所属包的 `git_blob` 解包，再取对应路径。
另 12 处候选差异经归一化确认只是换行符差异，可由历史提交还原。

## 保留

| 归档 | 条目 | 状态 |
| --- | --- | --- |
| `full_rs_eta_augmented_20261029_v1/eta_augmented_results_v1.zip` | 67 | `historical_superseded`——v15 正式图层（76 份 PDF）的数值源，已被 `ward_moller_sphere_20261121_v1` 取代；保留为历史数值记录，不作现行数据源 |

## 消费方说明

以下三处仍引用已删文件名，但都只服务 2026-10 批次的重跑，且该批次在
`config/governance/task_tracks.toml` 中已归档；CI、`tests/`、`src/` 与活跃脚本均不调用它们
（`git grep` 检查均为 0 命中）：

- `scripts/analysis/relaxtime/ThresholdUpgradeCommon.jl`：从 `package_manifest.json` 所在目录推导
  `source.zip`，缺失时报 `ArgumentError("Missing frozen source archive")`；
- `scripts/analysis/relaxtime/run_transport_batch_queue.py`：构建 dispatch 包时复制
  `pilot/calculation-source.zip`；
- 本目录内 `full_rs_dispatch_preparation_v1/prepare_calculation_commit.py`：把
  `full_rs_runner_review_checks_v1` 的 `source.zip` 列为该批次的历史必需证据文件。

需要重跑时先按上面的还原命令取回文件，或改用对应的 Ward/Møller 现行流程。
