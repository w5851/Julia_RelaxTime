# `Models.run_external_path_meson_density_scan`

本页说明沿外部显式离散路径点列运行介子数密度 workflow 的统一入口。

## 定位

`Models.run_external_path_meson_density_scan` 负责：

- 接收外部提供的 `(T,\mu_B)` 离散点列
- 保持点列顺序上的 continuation 语义
- 通过 flavor chemical profile 固定 `\mu_u`、`\mu_d`、`\mu_s`
- 通过 meson chemical profile 固定 `\mu_\pi`、`\mu_K`、`d_\pi`、`d_K`
- 复用统一 `MesonDensityWorkflow` 物理核，输出路径扫描结果

它的角色不是新的 solver，而是“**外部路径驱动**”的正式 workflow 入口。

当前首要用途是：

- `Friesen 2019` phase-line 一类文献离散路径复现
- digitized / extracted literature path 上的 `K/\pi` 工作流复跑
- 需要先固定外部路径，再比较物理量而不是内部路径生成器的场景

## 何时使用它

当路径不是：

- 内部规则 T-μ 网格
- freeze-out 参数化
- 内部 crossover locator 自动构建出的连续线

而是来自：

- 外部 CSV
- 文献提取点
- 人工整理路径点

时，应优先使用这条入口，而不是把这些点再手工塞进脚本层循环。

## 输入契约

入口要求 `points` 至少给出：

- `muB_MeV` 或 `muB_GeV`
- `T_MeV` 或 `T_GeV`

可选元数据包括：

- `source_fig` 或 `path_source`
- `case_id` 或 `path_case_id`
- `line_style` 或 `path_line_style`
- `point_index`
- `path_order_key`

若未显式给出 `path_order_key`，默认用 `mu_B` 顺序排序并推进 continuation。

## 当前支持的物理口径

与 freeze-out / crossover meson-density 入口保持一致：

- `:stable`
- `:strict_bw_stage1`
- `:strict_bw_stage2`
- `:phase_shift_current`
- `:phase_shift_gbu`

因此它不是额外维护的一套物理实现，而只是复用了统一密度 workflow 的另一种路径壳。

## continuity 契约

这条入口当前同样复用：

- `MesonMassWorkflow` 返回的 `continuation_state`

这意味着“沿外部路径逐点推进”的连续性机制已经被下沉到正式 workflow 层，而不是继续散落在脚本里。

## 输出合同

CSV 头部固定包含路径元数据：

- `path_source`
- `path_case_id`
- `path_line_style`
- `path_point_index`
- `path_order_key`

以及与其他 meson-density 扫描一致的核心结果：

- `muq_MeV`
- `muB_MeV`
- `T_MeV`
- `T_over_muB`
- `mu_u_MeV`, `mu_d_MeV`, `mu_s_MeV`
- `mu_pi_MeV`, `mu_K_MeV`
- `m_pi_MeV`, `m_K_MeV`
- `gamma_pi_MeV`, `gamma_K_MeV`
- `n_pi`, `n_K`, `kpi_ratio`

## 当前边界

这条入口当前只解决：

- 外部离散路径如何进入统一 workflow
- continuity 如何在外部路径上保持
- 结果如何带着路径元数据稳定落盘

它不负责：

- 自动从文献图像提取路径
- 自动构建连续插值路径
- 自动完成 target comparison / plotting

这些步骤应在脚本或 analysis 层完成，然后把“最终点列”交给本入口。

## 相关主题

- [Overview.md](Overview.md)
- [FreezeoutMesonDensityScan.md](FreezeoutMesonDensityScan.md)
- [../../../dev/archived/2026-05-04_Friesen2019曲线验证口径说明.md](../../../dev/archived/2026-05-04_Friesen2019%E6%9B%B2%E7%BA%BF%E9%AA%8C%E8%AF%81%E5%8F%A3%E5%BE%84%E8%AF%B4%E6%98%8E.md)
