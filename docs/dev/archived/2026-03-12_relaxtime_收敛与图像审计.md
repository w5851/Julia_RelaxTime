---
title: RelaxTime 收敛与图像审计记录
archived: true
original: docs/dev/active/2026-03-12_relaxtime_收敛与图像审计.md
archived_date: 2026-03-24
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# RelaxTime 收敛与图像审计记录

## 背景

本记录用于跟踪本轮手动 workflow 生成结果中的三类问题：

- 收敛性结论是否已经足够支撑当前 literature 对比图；
- 图像坐标轴是否存在绘图层面的视觉缺口；
- `plan_b` 在 `T=190 MeV` 与 `T=200 MeV` 的波动究竟来自绘图还是来自数值结果本身。

## 本轮结论

### 1. 当前收敛性结论

- 本轮 literature 对比采用的高精度参数已经明显优于早期设置，主要提升来自 `tau_angle_nodes=8`、`tau_phi_nodes=8`、`tau_p_nodes=28`、`sigma_grid_n=128`。
- `plan_a` 主体曲线已经满足当前对比用途，早期大幅毛刺未再出现。
- 剩余异常主要集中在 `plan_b` 的局部 `xi` 扫描点，而不是整条温度扫描都失稳。

### 2. 图像层问题与数值层问题需要分开处理

- `data/outputs/figures/relaxtime/plan_a/sigma_over_T_over_eta_over_s_vs_T_MeV.png` 的“坐标轴似乎缺失”主要是绘图脚本导致的视觉问题。
- 根因是 `scripts/plot_scan_csv.py` 曾将左/下 spine 裁到首末主刻度，再叠加 `bbox_inches="tight"` 导出，容易让坐标轴在边缘看起来像断开。
- 该问题不对应数据缺失，属于展示层缺陷。

### 3. `plan_b combined eta_over_s` 的可读性问题

- `data/outputs/results/relaxtime/plan_b/plan_b_merged.csv` 中，`T=150 MeV` 的 `eta_over_s` 量级约为 `12.63` 到 `20.31`。
- 同一图中其余温度仅约为 `0.06` 到 `0.36`。
- 因此线性坐标下，`T=150 MeV` 会把 `T=190/200/250 MeV` 的曲线压扁。
- 该问题是跨温度量级差导致的展示问题，不是单条曲线失真。

### 4. `T190/T200` 波动属于真实数值波动

- 对 `transport_vs_xi_T190_muB0.csv` 和 `transport_vs_xi_T200_muB0.csv` 的直接检查表明，波动已存在于 CSV 原始数据中，不是绘图伪影。
- `T=200 MeV` 明显比 `T=190 MeV` 更不稳定：
  - `sigma_over_T` 的局部转折次数约为 `21`，而 `T=190 MeV` 约为 `6`；
  - `zeta_over_s` 的局部转折次数约为 `13`，而 `T=190 MeV` 约为 `5`；
  - `eta_over_s` 的局部转折次数约为 `13`，而 `T=190 MeV` 约为 `5`。
- 进一步比较状态量与输运量的相邻点差分可见：
  - `Phi` 在 `T=200 MeV` 的最大相邻变化仅约 `0.0022`；
  - `m_u` 的最大相邻变化约 `0.034`；
  - `tau_u` 的最大相邻变化却达到 `0.43`，集中在 `xi≈0.30-0.38`；
  - `T=200 MeV` 在 `xi≈-0.28` 附近也存在明显的 `tau` 跳变。
- 新增的 `scripts/analysis/relaxtime_xi_fluctuation_study.jl` 在坏区间上给出更强证据：
  - `T=200 MeV` 的 `xi=-0.32,-0.30,-0.28,-0.26,0.30,0.32,0.34,0.36`；
  - `T=190 MeV` 的 `xi=-0.46,-0.44,-0.42`；
  - 这些点上 `models` backend 都未收敛，实际计算全部回退到了 `legacy` backend。
- 在这些坏点上，`current` 与 `reference` 两套 `tau` 数值配置的差异通常只有轻味道扇区约 `0.1%` 到 `4%`，而相邻 `xi` 点之间的 `tau_u`/`tau_ubar` 跳变量却可达到 `0.07`、`0.16`、`0.36` 这一级别。
- 因而当前更像是“平衡态求解 backend 失稳 + 分支连续性破坏”主导了 `tau(xi)` 尖刺，散射率积分分辨率只是次要修正项。

## 已执行修复

- 已在 `scripts/plot_scan_csv.py` 中取消 spine 截断，并适度放宽导出边距。
- 已将 `plan_b combined eta_over_s` 的绘图改为对数 y 轴，以保留所有温度曲线的可读性。
- 已重绘以下受影响图像：
  - `plan_a/sigma_over_T_over_eta_over_s_vs_T_MeV`
  - `plan_b/combined/eta_over_s_vs_xi`
  - `plan_b/T190/sigma_over_T_vs_xi`
  - `plan_b/T200/sigma_over_T_vs_xi`
  - `plan_b/T200/zeta_over_s_vs_xi`

## 待继续跟进

- [ ] 优先检查 `plan_b` 的平衡态解连续性：沿固定 `T` 的 `xi` 扫描复用前一点解、引入相变感知 seed，避免坏区间直接退回到不连续的 `legacy` 解支。
- [ ] 在 `tau` 计算路径继续输出更细诊断量，至少包括主导过程贡献、平均散射率与插值前后 `sigma(s)` 行为，用于区分“解支切换”与“通道积分误差”。
- [ ] 对 `T=200 MeV` 的 `xi≈-0.28` 与 `xi≈0.30-0.38` 做局部加密扫描，但将其视为二级验证，而不是首要修复手段。
- [ ] 评估是否需要为 `plan_b combined` 补充一张 small-multiples 版本，避免后续再次依赖对数坐标解释局部趋势。

## 验收口径

- 若沿 `xi` 复用连续 seed 后，`models` backend 不再大面积失败且 `tau(xi)` 尖刺显著减弱，则说明主因是平衡态分支连续性问题。
- 若即便解连续后波动仍保留，再继续收紧 `sigma_grid_n`、`phi_nodes` 与自适应 `sigma(s)` 细化阈值，判断剩余部分是否来自散射率积分误差。