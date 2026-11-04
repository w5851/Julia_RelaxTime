# Full-RS 完整网格执行与验收

状态：2026-10-26 完整网格数值验收通过。135 个新分片和七个复用先导均已完成，
原始产物、合并表和最终审计已取回并核对哈希。所有 1819 点及 30923 道诊断通过，
没有缺点、重复点或失败点；[最终验收](final_review_v1/full_grid_acceptance.json)由 quote-remote 生成。
实际运行及归档见 [收集状态](collection_status.json)，完整点集与逐片检查见
[网格复核](returned/transfer_090/reviews/transfer_090/reviewed_grid.json)。

- 计算提交：`2e5c411c716d5cbc524c3fc6cfb335fe05b39444`。
- 执行标签：`transport-full-rs-20261026-v1`，workflow ID `291078197`。
- 数值配置：`thermal_weighted_v1`、二维完整 RS A、full-A 耦合、`cm` 速度。
- mode-A：910 点，复用先导 6 点，剩余 904 点分为 72 片。
- mode-B：909 点，复用先导 5 点，剩余 904 点分为 63 片。
- 合计：1819 点、每点九项最终量与 17 道诊断；完整诊断数为 30923。

## 已验收的合并表

| 模式 | 点数 | 结果 | 通道诊断 |
| --- | ---: | --- | --- |
| A | 910 | [scan.csv](returned/transfer_090/reviews/transfer_090/mode_a/phase_guided_transport_scan.csv) | [sampling.csv](returned/transfer_090/reviews/transfer_090/mode_a/phase_guided_transport_scan_sampling.csv) |
| B | 909 | [scan.csv](returned/transfer_090/reviews/transfer_090/mode_b/phase_guided_transport_scan.csv) | [sampling.csv](returned/transfer_090/reviews/transfer_090/mode_b/phase_guided_transport_scan_sampling.csv) |

[逐点来源](returned/transfer_090/reviews/transfer_090/point_provenance.csv)关联原始 run、归档和 CSV 哈希。
合并表保留原始数值字符串，原始分片继续独立留存。

[完整背景审计](final_review_v1/background/background_audit.json)通过 1819 点：
质量和熵的最大相对差分别为 `2.17057e-8`、`1.91574e-8`，可区分稳定分支不一致为零。
相邻变化筛出的 10 项全部位于已验证的一阶边界：mode-A 的 μB=900、αT=1、ξ=±0.003，
以及 mode-B 的 T=120、μB=900、ξ=−0.13/−0.12；非边界筛查项为零。
470 条 `tau_u_ubar_ratio_high` 有限比值标记均已核对定义并原样保留。

## 实际 runner 用量

| 运行集合 | jobs | 点数 | 累计 job 小时 | 最长 job [s] |
| --- | ---: | ---: | ---: | ---: |
| A 新分片 | 72 | 904 | 34.0192 | 2558 |
| B 新分片 | 63 | 904 | 32.6736 | 2523 |
| 复用先导 | 7 | 11 | 1.2864 | 809 |
| 全部 | 142 | 1819 | 67.9792 | 2558 |

135 个新分片合计 66.6928 runner 小时，最长 42.63 分钟，均低于 180 分钟规划目标和
360 分钟 job 上限。扫描进程峰值 RSS 为 1,945,735,168 字节（约 1.81 GiB），
占对应主机可见内存约 11.61%。[计时汇总](final_review_v1/runner_accounting.json)及
[逐 run 记录](final_review_v1/runner_accounting.csv)在 quote-remote 生成。
JIT 已包含在阶段耗时中；RSS 只覆盖计时的 Julia 扫描进程，累计 job 时间不等于账户账单，
也不作为逐点稳态计时。首个新 job 开始至最后一个结束为 4 小时 58 分 10 秒，包含派发恢复间隙。

## 来源与恢复记录

[queue_plan.json](queue_plan.json) 冻结请求和先导输入的哈希。
[dispatch_manifest.json](dispatch_manifest.json) 是持续更新的实际运行索引；
`attempts/` 保存请求提交记录，`runs/<id>/` 独立保留原始 API 元数据、ZIP 和解压字节。
`pilot_runs/` 保留已验收先导。各片不得直接解压到同一个结果目录。

本机控制器只执行 Actions 调度、打包、下载及文件身份检查，最多保持 20 个未结束任务。
`relay_reviews.py` 将产物传到 quote-remote 的 `D:/CodexCompute/trfull-20261026-v1`，
逐次复核报告及回执返回 `returned/` 和 `transfers/`。完整网格通过后，
`finish_remote_review.py` 在 W 完成背景审计和计时汇总，返回包见
[最终回执](final_grid_v1.receipt.json)。派发、收集及最终审计进程均已结束。

控制器的十项行为测试、全网格复核器的六项测试及背景审计的三项测试均在 quote-remote 通过。
初版控制器在 POST 前发现 API 布尔值使用字符串；修订校验后沿用冻结请求。
初版全网格复核将 `tau_u_ubar_ratio_high` 启发式标记误作数值失败；修订版核对其真实定义，
保留高比值标记并交给背景和局部变化审计。控制器 v3 补齐最后一片归档延迟可见时的等待；
v4 在 Actions 只读查询发生网络错误时保留状态并继续轮询，POST 仍不自动重发。
传输器增加同一冻结包的重连重试。原版本和失败报告保留，数值任务未因此重算。

两次 POST 网络中断均在完整运行列表确认缺失后，用原冻结请求单次补提交：
[mode-A 恢复](dispatch_recovery/full_modeA_muB450_alpha1.1_002/resolution.json)和
[mode-B 恢复](dispatch_recovery/full_modeB_T120_muB900_004/resolution.json)。
最终每个请求恰好对应一个 run。传输中断及一次状态文件暂时不可读的日志保留，
恢复继续消费已冻结的产物；[接续记录](supervision_handoff.json)记录执行边界。

背景审计复用 [background_inputs/](background_inputs/) 中带原归档哈希的 910/909 点背景证据，
沿用质量/熵的 1% 门槛并核对可区分的稳定分支。相邻变化的 50% 筛查仅选取诊断对象，
不代表积分误差，不修改收敛容差，也不跨越未返回的网格缺口。

完整网格数值验收已通过；本批表已生成
[full-RS 原始 PNG 审阅包](../../phase_guided/full_rs_raw_20261027_png_review/README.md)，
与 v13 同构的 75 对彩色/灰度图和逐图质量检查均通过，全部原始点及有限比值提示保留。
没有平滑、插值、替换或删点，仅保留两处已审计的一阶断线。局部小凸起仍原样供作者查看；
PNG 作者接受、PDF、正式产物晋升与论文使用资格按各自合同记录。
当前冻结计算的代表点收敛、
v3 回归、真实默认入口和 ξ=0 外部对齐证据由
[生产接入任务](../../../../../../../docs/dev/active/2026-10-23_输运采样生产接入与重跑合同.md)维护。
计算源码、配置和物理口径未改变，本批复用这些资格证据；A 解析降维保持独立待办。
