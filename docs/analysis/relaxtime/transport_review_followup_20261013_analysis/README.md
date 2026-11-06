# 输运重计算前的局部研究

本目录保留 2026-10-14–16 两项仍有独立用途的研究。旧尾项、Julia–Fortran
对齐和整批生产待办已由 PR #331–335 完成；当前数值、v15 图件和未决问题见
[统一状态](../threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)。

## k=0 完整 RS 与线性 B₀

使用相同的旧生产背景、isotropic A 有效耦合、普通在壳分母和实轴宽度约定，
对 RS 占据作角平均，并以主值减法和解析留数计算相对各向同性 B₀ 的修正。
这是占据的完整 RS 比较，范围为 k=0。

[48 个 π 点对照](../../../../data/outputs/results/relaxtime/transport/review_followup_20261013/run_20261014/tables/full_rs_pi_comparison.csv)
的质量最大相对变化为 0.5062%，宽度最大绝对变化 2.2316 MeV；
μB=0、αT=1.1、ξ=0.47 处 Γ 从 7.8434 到 8.5105 MeV（8.5049%）。
这表明近 Mott 的窄宽度可比质量更敏感。旧背景和耦合限制了向 v15 定量外推的范围。

[原始极点表和运行参数](../../../../data/outputs/results/relaxtime/transport/review_followup_20261013/run_20261014/meson_v2/)
与背景输入一并保留。源码包保留对应 FullRSKZero、运行器和测试。
有限 k 的非共线完整 RS、复能量解析延拓仍需专门研究；同节点消融也只量化该项
模型依赖，不能给未计算高阶项提供普遍误差上界。

## Mott 支点与采样偏差

[推导和受控比较](mott_threshold_mechanism.md)保留 P/S 阈值阶数、混合扇区投影、
束缚极点/留数以及带热权积分的可积性。旧 +0.37 窗口的局部残差
−6.5292% → −0.06466% 帮助识别采样误差，后续由完整输运重计算接续。

## 材料与收尾

[关键材料](../../../../data/outputs/results/relaxtime/transport/review_followup_20261013/README.md)
保存最终表格、参数和一个去重源码包。未把过期执行日志、逐点检查点、旧图、
重复源码树或已完成路线重新接入主干。桌面旧六个提交由标签 `archive/desktop-docs-20261106` 保留在 Git 历史，
旧 FIG4/APS 图包可按 `5946aaf73a4d26f0f392bc7d43b49419f0ee6f26` 恢复。
