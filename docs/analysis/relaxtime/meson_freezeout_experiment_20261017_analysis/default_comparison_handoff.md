# 化学冻结线默认比较的工作树交接

当前状态（2026-10-23）：实现、默认比较、文献筛选及两族探索已随 [PR #327](https://github.com/w5851/Julia_RelaxTime/pull/327) 合并。
研究与图像保持 `review`，等待下周组会后的作者决定。
当前 ledger item：`meson-freezeout-parameter-envelope`。默认参数阶段保留为
`meson-freezeout-default-comparison`。作者于2026-10-18同意进入带状化，并允许跨文献
汇合参数范围及组合系数。研究继续使用q0外推，不把crossover直接改称化学冻结线。

实现、最终图、数值差异和验证统一见
[仓库内研究链报告](../meson_freezeout_defaults_20261017_analysis/README.md)。
原始计算工作树为 `D:/w/meson-freezeout-defaults`，计算基线
`698feaf8fada670fde81c93bab6093b50093935f` 保留为历史身份；2026-10-23 同步时的提交基线为 `dffe4070eb5c0788e9c8a1e227cf40978b51ce3e`。

历史[默认参数补筛](../meson_freezeout_defaults_20261017_analysis/literature_candidates_20261017.md)
已加入 IHRG 2019、FRG 2024 的两个场景和 functional QCD 2026 mean；标准/加密
共72点完成。共同8个 Au+Au 点上 FRG mu-independent 的 RMSE 最低，但原文未
施加奇异数中性，7.7 GeV也未用于其冻结线拟合。旧14点评分和原图保留。

当前的[跨文献参数包络研究](../meson_freezeout_defaults_20261017_analysis/parameter_envelope_20261018.md)
与 [Andronic 族报告](../meson_freezeout_defaults_20261017_analysis/andronic_family_20261018.md)
统一保存原文约束、采样设计、数值结果和图像入口。下文初始设计为46条完整冻结线的说明保留为历史阶段，
含原文参数组、明确的低温混合场景及包络内部样本。主评分固定为7.7–200 GeV的
九个Au+Au点，不与旧14点指标混比；最优样本和包络均不作统计置信区间解释。

原根工作树保持既有改动；输运 FIG4 和 paper 输运研究继续由其独立台账与原对话管理。
首轮文献与快度口径见 [研究入口](README.md)。
