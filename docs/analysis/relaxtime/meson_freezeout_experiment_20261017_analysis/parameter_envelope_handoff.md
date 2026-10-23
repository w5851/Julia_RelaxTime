# 化学冻结线参数包络交接

Ledger item：`meson-freezeout-parameter-envelope`。作者于2026-10-18同意带状化，
允许汇合不同文献的参数范围并进行跨源组合。主线仍是介子数密度，使用q0外推。

状态（2026-10-23）：`review`。两族各59组参数探索、每族三组standard/fine核验及彩色/灰度PNG已完成，
实现和冻结证据已随 [PR #327](https://github.com/w5851/Julia_RelaxTime/pull/327) 合并；图像继续等待下周组会后的作者决定。

设计、来源、结果和验证统一从[研究链主报告](../meson_freezeout_defaults_20261017_analysis/README.md)进入；
[多项式族](../meson_freezeout_defaults_20261017_analysis/parameter_envelope_20261018.md)与
[Andronic 族](../meson_freezeout_defaults_20261017_analysis/andronic_family_20261018.md)分别保留数值与图像入口。
计算位于 `D:/w/meson-freezeout-defaults`；原始基线 `698feaf8fada670fde81c93bab6093b50093935f` 保留为历史身份，
2026-10-23 同步时的提交基线为 `dffe4070eb5c0788e9c8a1e227cf40978b51ce3e`。

历史默认参数比较见[原交接](default_comparison_handoff.md)。不同实验能区的评分
不能直接混比；探索包络与置信区间分开。输运FIG4及paper输运研究保留其独立台账。

报告中的当前建议为同时展示最佳样本与探索带。峰位及幅度有改善，但高能比值和
horn尖锐程度仍有偏差；跨源混合后的muB映射可超出全部原文中心曲线，独立
物理约束和后续推进待组会决定。所有数值、失败、实验输入勘误和图像链接均以研究链报告为准。
