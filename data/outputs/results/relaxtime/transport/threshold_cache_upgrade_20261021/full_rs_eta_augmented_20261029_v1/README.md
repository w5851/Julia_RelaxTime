# η 阈值局部加密正式数据（2026-10-29）

状态：`formal_numerical_data_accepted`。本版本按作者同意的方案纳入48个新ξ点及9个原坐标重算值。
计算核仍为 `2e5c411c716d5cbc524c3fc6cfb335fe05b39444`，全部数值回放及派生审查在 quote-remote 完成。
原v14网格及η诊断归档保持冻结；这是独立数值版本。

| 数据 | 入口 |
| --- | --- |
| mode A，942点 | [主表](returned/grid/mode_a/phase_guided_transport_grid.csv) · [17道采样诊断](returned/grid/mode_a/phase_guided_transport_grid_sampling.csv) |
| mode B，925点 | [主表](returned/grid/mode_b/phase_guided_transport_grid.csv) · [17道采样诊断](returned/grid/mode_b/phase_guided_transport_grid_sampling.csv) |
| 1867点逐行来源 | [point_provenance.csv](returned/grid/point_provenance.csv) |
| 9原坐标变化 | [recomputed_original_changes.csv](returned/grid/recomputed_original_changes.csv) |
| 正式数值资格 | [promotion.json](formal_promotion_v1/promotion.json) |
| 合并与独立复核 | [grid review](returned/grid/review.json) · [assembly verification](returned/assembly_verification.json) |
| 接缝割线 | [boundary_secants.csv](returned/grid/boundary_secants.csv) |
| 57点EOS与τ回放 | [summary](returned/thermo/summary.json) · [checks](returned/thermo/replay_checks.csv) |
| 返回完整性 | [local_return_verification.json](local_return_verification.json) · [package manifest](returned/package_manifest.json) |

1867点由1810个原值、9个原坐标重算值及48个新增点组成，共31739条通道记录。
三条加密曲线为 mode A μB=0/450 MeV、αT=1.0，以及 mode B T=200 MeV、μB=0。
57点同坐标优先选用已完成的更严格结果：39点s积分表rtol=10⁻⁵，18点rtol=10⁻⁶；
其他生产积分设置继承v14。实际背景残差要求10⁻¹⁰、max_iter=80、稳定多种子选解，
与profile文件内默认残差门槛10⁻⁸的区别已记录。所有点均为计算结果，没有平滑、插值或值混合。

独立复核确认1810个原点及30770条对应通道记录逐字段不变；57个选用点的12个绘图量、
质量、密度、τinv和背景变量与来源记录逐字段一致。由保存的凝聚场和Polyakov变量补齐EOS附属列，
质量、熵、热力学势和17道率重建的六种τ及比值全部一致，未调用新gap求解或截面积分。
9个重算点的12个绘图量最大变化为0.009043925%，原生产最终量目标为1%。

470个有限τ比值提示以及两处已审计一阶断线保持不变。6个新旧网格接口检查了τū、η/s、ζ/s、σ/T
的24条割线和相邻区间，没有新的反号或突跳；该检查不提供连续导数误差的统一上界。
57点未保存的求根迭代数留空；9原坐标保留旧相态追踪元数据，48新点标记unknown。
没有重做相态分类，也没有把η反向Mott等同于热力学一阶相变。

本轮复用同一冻结源码的8项生产/回归资格，补入已有ξ细化及s容差比较。
没有重新认证三处全部外层/t/B₀积分或连续导数极值，也未完成η分振幅消融。
回放前后340份源码hash一致；下载归档66文件的大小和hash全部匹配。
返回ZIP的SHA-256为 `213bc222dc2d7e821454c4357bec6cfb7d548ff33bd1202bd28e167d848b893f`。

`grid_attempt_v1`保留第一次合并的未验收输出：旧背景审计按观测量重复记录断线，
初版组装器误把记录行数当成断线数。修正为按物理坐标去重并增加测试后，最终`grid`通过核验；
此修正没有更改数值结果。

新图已生成于同构的[publication_clean_v15_png_review](../../phase_guided/publication_clean_v15_png_review/README.md)。数值正式资格与新PNG作者接受分开记录，
该PNG阶段不自动授予PDF交付或替换当前论文图层。当前分析入口为
[Mott局部结构](../../../../../../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/mott_local_bumps_20261028.md)。
