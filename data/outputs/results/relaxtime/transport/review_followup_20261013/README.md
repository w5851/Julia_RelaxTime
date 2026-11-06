# 早期输运研究的关键材料

2026-10-14–16，代码基线 `d6282e456cbda0d406809b239a4f2e62b35ca6cd`，
Julia 1.12.5。科学说明见
[研究入口](../../../../../../docs/analysis/relaxtime/transport_review_followup_20261013_analysis/README.md)。

- `review_points.csv`：原背景与点位；`meson_v2/manifest.json` 固定其 hash。
- `run_20261014/meson_v2/`：原极点 CSV、参数和状态。
- `run_20261014/tables/full_rs_pi_comparison.csv`：48 点完整 RS/线性 B₀ 对照。
- `mott_threshold_mechanism_20261016/`：最终支点、束缚侧、缓存和同网格速率数据；
  `analysis_v1/` 保存拟合、汇总和原有 13 项输入 hash。
- `research_sources.zip`：相对上述 Git 基线的源码差异及 bound_v2/rates_v1/meson_v2 覆盖层，
  另保留桌面最终独有分析源码与测试。包内 README 说明恢复方式；源码按原字节去重。

这些材料使用旧有效耦合和积分设置。当前正式输入使用 v15 所指向的完整 A^ξ 数值，
后续重算及验收由[当前交付](../../../../../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)维护。
