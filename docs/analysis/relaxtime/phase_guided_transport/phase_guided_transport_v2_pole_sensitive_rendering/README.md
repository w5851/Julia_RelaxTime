# Phase-guided transport v2 极点敏感派生显示审计

## 范围与结论

本分析包消费 `first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`，不修改 production CSV、canonical figure 或 `production_registry.json`。它把旧 xi001 denominator-chain 机制表迁移到 v2，并生成**仅供作者审阅**的非破坏性派生显示候选。

当前结论：8 个旧窗口通过小分母机制迁移门槛；作者复核发现的 2 个 `mu_B=0` 残留窗口通过当前 v2 production 的新定点机制诊断；2 个一阶/上游分支窗口被硬保护。虚线桥接是视觉指南，不是新计算值，也不能静默替换正式图。

## 输入与迁移门槛

| mode | scan rows | diagnostic rows | max tau v1/v2 drift | max rate v1/v2 drift |
| --- | ---: | ---: | ---: | ---: |
| mode_a | 909 | 38178 | 0.00016650389 | 0.00053765629 |
| mode_b | 909 | 38178 | 8.2279794e-06 | 0.00037216686 |

迁移门槛为 `0.001`。v2 的能量语义修改不改变弛豫时间定义；实际 tau 和 channel-rate 漂移均显著低于门槛，因此允许继承已有 denominator-chain 定点证据。该门槛不替代尚未完成的局部 high-rate convergence gate。

新增的两个 `mu_B=0` 窗口不依赖 v1 迁移：`xi=0.37, alpha_T=1.0` 的轻味同步下探由 `mixed_detM` 小分母链支持，`xi=-0.47, alpha_T=1.2` 的奇异味下探由 `simple_1m4KPi` 小分母链支持；两者的上游质量、Polyakov loop 与熵背景均平滑。完整证据见 `supplemental_muB0_noise_mechanism/`。

## Pole-sensitive display mask 候选

| window | panel | series | xi | denominator branch | eligible observables |
| --- | --- | --- | ---: | --- | --- |
| mode_b_T200p0_muB450p0_xip0p31_1 | T200.0 | muB450.0 | 0.31 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;tau_sbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB450p0_alpha1p0_xip0p26_2 | muB450.0 | alpha1.0 | 0.26 | mixed_detM | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB0p0_alpha1p1_xip0p48_1 | muB0.0 | alpha1.1 | 0.48 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB900p0_alpha1p1_xip0p18_1 | muB900.0 | alpha1.1 | 0.18 | simple_1m4KPi | tau_ubar;tau_dbar |
| mode_a_muB450p0_alpha1p1_xip0p35_1 | muB450.0 | alpha1.1 | 0.35 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB900p0_alpha1p2_xip0p49_1 | muB900.0 | alpha1.2 | 0.49 | simple_1m4KPi | tau_ubar;tau_dbar |
| mode_a_muB450p0_alpha1p0_xim0p20_1 | muB450.0 | alpha1.0 | -0.20 | simple_1m4KPi | tau_ubar;tau_dbar |
| mode_b_T200p0_muB0p0_xim0p21_1 | T200.0 | muB0.0 | -0.21 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;sigma;sigma_over_T |
| mode_a_muB0p0_alpha1p0_xip0p37_supplement | muB0.0 | alpha1.0 | 0.37 | mixed_detM | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB0p0_alpha1p2_xim0p47_supplement | muB0.0 | alpha1.2 | -0.47 | simple_1m4KPi | tau_s;tau_sbar;eta_over_s;zeta_over_s;sigma_over_T |

规则：

1. tau 只有在旧分析的 `affected_tau_fields` 中才进入 mask；
2. transport 默认只有在目标点相对相邻两点均值的绝对 log 偏离不小于 `0.03` 时才进入 mask；作者在论文图中明确指出的两个新增残留窗口，经小分母机制确认后，对三个论文展示比值显式纳入 mask；
3. raw 点始终以橙色叉号保留；solid guide 在 eligible 点断开，并以虚线连接两侧真实邻点；
4. `linear_neighbor_guide_value` 只用于画图，不写回 production，也不声称是物理值。

## 一阶窗口硬保护

| window | panel | series | xi | phase structure |
| --- | --- | --- | ---: | --- |
| mode_b_T120p0_muB900p0_xim0p09_1 | T120.0 | muB900.0 | -0.09 | first_order_possible |
| mode_a_muB900p0_alpha1p0_xip0p00_1 | muB900.0 | alpha1.0 | 0.00 | first_order_possible |

这两个窗口保留 raw 曲线和跳变，不应用桥接或 mask。保护依据是背景质量、Polyakov loop 和熵密度的同步快速变化，而不是单纯依赖 phase 字符串标签。

## `zeta/s` 相变前回落的分支一致性审计

对 `mode_a, mu_B=900, alpha_T=1.0` 的 `xi=-0.02,-0.01,0.00`，使用与 production 主平衡态一致的 `p_num=12,t_num=6` 重新调用 `Models.bulk_viscosity_coefficients`，并从强子/夸克种子分别求候选根、比较热力学势。结果见 `tables/bulk_derivative_branch_audit.csv`：

- `xi=-0.02` 与 `xi=0.00` 时，bulk 导数路径的轻味质量与主 production 平衡态处于同一分支；
- `xi=-0.01` 时，主 production continuation 保留 `m_u=0.73435 fm^-1` 的低质量夸克候选，而 bulk 得到 `m_u=1.37470 fm^-1` 的高质量强子候选；两根的 `Omega_h-Omega_q=-1.2548e-3 fm^-4`，因此强子候选才是该离散口径下的稳定平衡态；
- `xi=-0.02 -> -0.01` 虽然 `tau_u` 上升且熵下降，`zeta` 本身却由 `1.85047` 降至 `1.49245`，因此回落来自 bulk `B^2` 核/热力学导数，而不是 tau 或除以熵造成。

当前 workflow 在取得主 equilibrium 后，另行调用不接收该 equilibrium/seed/branch 的 `bulk_viscosity_coefficients`。该回落因此被重新判定为**主 continuation 保留亚稳分支、bulk 选择稳定分支后造成的混合分支结果**，不是已确认的物理非单调趋势，也不属于传播子小分母显示噪点。本 PR 不平滑这一点；后续代码必须先按热力学势选择主稳定态，再让 bulk 复用同一 `base_state`。

旧分析 `phase_guided_transport_p128_xi001_analysis` 只记录了 `xi=0` 主平衡态的一阶跳变，以及 tau 上升和熵下降对 `zeta/s` 跳升的放大；没有在同一点比较两个候选根的热力学势，因此没有覆盖本次发现的亚稳 continuation。

## `alpha_T=1.0` 相变锚点审计

结果见 `tables/phase_anchor_coexistence_audit.csv`。当前 `mu_B=900 MeV` 的 `T=125.06992 MeV` 来自旧 `data/reference/pnjl/boundary.csv` 在 `T=110,130 MeV` 两点之间的线性插值；在相同 `p_num=12,t_num=6` 口径下直接求 `xi=0` 两分支等势，得到 bracket 中点 `T=125.76661 MeV`，比旧锚点高约 `0.69669 MeV`。

在该共存温度 bracket 的上下端，`xi=-0.003` 均稳定为夸克候选，`xi=+0.003` 均稳定为强子候选，首轮双侧夹逼认证通过。该结果支持后续采用“共存点不输出唯一输运量、以认证后的两侧近邻表示单边结果”的设计，但正式 production 前仍需完成双分支连续追踪与热力学积分节点收敛审计。旧 `boundary.csv` 可保留为相图与初始 bracket 证据，不再作为导数敏感 production 的唯一精确锚点。

## 派生图

- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_b_T200p0_muB450p0_xip0p31_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB450p0_alpha1p0_xip0p26_2.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB0p0_alpha1p1_xip0p48_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB900p0_alpha1p1_xip0p18_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB450p0_alpha1p1_xip0p35_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB900p0_alpha1p2_xip0p49_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB450p0_alpha1p0_xim0p20_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_b_T200p0_muB0p0_xim0p21_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB0p0_alpha1p0_xip0p37_supplement.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB0p0_alpha1p2_xim0p47_supplement.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/first_order_protected_mode_b_T120p0_muB900p0_xim0p09_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/first_order_protected_mode_a_muB900p0_alpha1p0_xip0p00_1.png`

每个 pole-sensitive 图包含 primary tau 与 `eta_over_s/zeta_over_s/sigma_over_T` 四个局部面板。每个一阶图包含相同的下游比值与 `tau_u`，用于核对跳变未被删除。

## 论文候选图

- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB0.0/eta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB0.0/zeta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB0.0/sigma_over_T_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB450.0/eta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB450.0/zeta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB450.0/sigma_over_T_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB900.0/eta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB900.0/zeta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_a/plot_panel=muB900.0/sigma_over_T_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T120.0/eta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T120.0/zeta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T120.0/sigma_over_T_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T160.0/eta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T160.0/zeta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T160.0/sigma_over_T_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T200.0/eta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T200.0/zeta_over_s_vs_xi.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/paper_figures/mode_b/plot_panel=T200.0/sigma_over_T_vs_xi.png`

论文候选图沿用正式图的固定 panel、多曲线布局，只绘制 `eta_over_s`、`zeta_over_s` 和 `sigma_over_T`：

1. 共生成 `18` 张 600 DPI 图；
2. `19` 个已确认的小分母下游噪点在派生绘图值中由左右相邻真实样本线性插值替换，图面只显示正常连续实线，不显示叉号、空心点、虚线或修正标签；
3. `6` 个 observable-level 一阶位置以星号标在对应曲线上，不使用竖直虚线；
4. 星号处仍使用 raw production 值，一阶/上游分支跳变没有被平滑；
5. 替换值与星号位置分别记录在 `tables/paper_display_replacements.csv` 和 `tables/paper_first_order_markers.csv`。图面不呈现内部修正痕迹，但仓库内保留完整可追溯记录。
6. `mode_a/plot_panel=muB900.0/zeta_over_s_vs_xi.png` 的 `alpha_T=1.0, xi=-0.01` 已知存在 bulk 导数分支不一致；该图保留用于作者审计，但在代码修复和 production 重跑前不具备论文输入资格。

## 证据边界与作者判断

- `supported`：v2 输入完整、无 failed/NaN/负 rate，v1→v2 tau/rate 迁移门槛通过；一阶窗口保护规则明确。
- `supported_with_scope_limit`：10 个窗口已有 denominator-chain 与生产 rate 复现证据，但没有逐窗口新的 high-rate convergence gate。
- `author_directed_candidate`：论文候选图按作者约定隐藏数值修正痕迹，并用星号标示一阶相变位置；仍需最终视觉审核后决定是否采用。
- `implementation_issue_supported`：`xi=-0.01` 的 bulk 导数分支不一致已由质量和导数定点复算支持；它不应被包装为物理趋势或普通显示平滑。
- `author_check`：是否需要计算层有限宽度/极点正则化和新 production slug。
- 本包不把小分母结构直接定性为随机数值噪声，也不把虚线桥接升级为物理预测。

## 复现

```powershell
julia --project=. scripts/analysis/relaxtime/audit_phase_guided_bulk_branch_consistency.jl
python scripts/analysis/relaxtime/build_phase_guided_pole_sensitive_rendering.py
```

关键表：

- `tables/input_inventory.csv`
- `tables/window_classification.csv`
- `tables/rendered_point_audit.csv`
- `tables/pole_sensitive_mask.csv`
- `tables/first_order_protection.csv`
- `tables/paper_display_replacements.csv`
- `tables/paper_first_order_markers.csv`
- `tables/bulk_derivative_branch_audit.csv`
- `tables/phase_anchor_coexistence_audit.csv`
- `tables/claim_ledger.csv`
- `figures/plot_manifest.json`
- `paper_figures/plot_manifest.json`
