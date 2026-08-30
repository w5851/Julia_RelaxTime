# PNJL Thermal Integral Legacy Workflow Notes

本文档记录隔离区 `Julia_test` 仓库中的 PNJL 热项积分实验 workflow，供以后判断是否值得在 `Julia_RelaxTime` 主线中重做。本文只做能力沉淀与迁移边界说明，不表示这些旧脚本已经成为当前主项目的有效入口。

## 来源与状态

- 隔离区来源（已删除）：`D:\Desktop\_cleanup_quarantine\2026-05-27\Julia_test`
- 远端来源：`https://github.com/w5851/Julia_test.git`
- 当前状态：`main...origin/main` 无 ahead，但工作区 dirty。
- 本次审计使用的完整归档（已删除）：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\full_archives\Julia_test_worktree_2026-08-16.tar.gz`
- 完整归档 SHA-256：`8CA458C13D9BD58CBB09166C02D35A08BE68DE0B4F37A05FA5991B3B84EE7EFE`
- 本次复核日期：`2026-08-16`
- 审计解包目录（已删除）：`D:\Desktop\_cleanup_quarantine\2026-05-27\_archive_review_2026-08-16\Julia_test`
- tracked dirty 文件：
  - `Manifest.toml`
  - `Project.toml`
  - `integral_test/Pnjl_pure.jl`
  - `integral_test/test_focus_transform.jl`
- 关键 untracked 文件：
  - `integral_test/Pnjl_minimal.jl`
  - `integral_test/Pnjl_pure_nofreeze.jl`
  - `integral_test/热项对数函数性质分析*.md`
  - `scripts/verify_analytic_ppeak_pmax.jl`
  - `scripts/bayesopt_estimate_pmax.jl`
  - `scripts/bench_node_freeze_vs_fixed.jl`
  - `scripts/regression_pnjl_freeze_vs_nofreeze.jl`
  - `tests/test_*`
  - `outputs/data/*.csv`
  - `outputs/figures/*.png`
- 安全导出（已删除）：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\Julia_test`

清理前的轻量 `_safety_exports` 只保存 git bundle、tracked diff、status 和 untracked 文件列表；本次完整归档曾包含 `integral_test/`、`scripts/`、`tests/` 和 `outputs/` 下的 untracked 正文。完整归档和轻量安全导出均已在审计完成后删除。

## Workflow 意图

这套旧 workflow 的核心问题是：PNJL 热项动量积分在低温、高化学势或轻质量区域会出现明显的费米面结构和尾部截断敏感性；如果动量节点、积分上限或节点位置随参数变化而直接参与 ForwardDiff，可能带来导数不连续或 AD 代价异常。

它尝试回答四个问题：

1. 热项被积函数的主要贡献区域在哪里，如何估计 `p_peak` 与有效尾部上界 `p_sup`。
2. 是否可以用解析近似或离散边界策略得到稳定的 `p_max`。
3. 如何冻结节点位置和权重，使 ForwardDiff 只对物理被积函数求导，而不对节点选择或节点映射求导。
4. 冻结策略对数值值、导数和运行时间的影响是否可接受。

## 旧代码组件

### 最小 PNJL Omega 与 AD 友好积分

`integral_test/Pnjl_minimal.jl` 是最完整的可读草稿。它把 PNJL 热项积分拆成几个层次：

- `MINIMAL_LEGENDRE_CACHE`：缓存标准 `[-1, 1]` Gauss-Legendre 节点。
- `MINIMAL_BOUNDARY_CACHE`：缓存离散边界 `[250, 500, 1000, 2000, 4000, 6000] MeV` 下的映射节点。
- `estimate_pmax(...)`：用参数化公式估计连续 `pmax`。
- `analytic_p_max(...)`：用 `ForwardDiff.value` 冻结 `mass/T/mu/Phi1/Phi2`，再用低温快捷路径或离散边界节点选择积分上限。
- `calculate_thermal_term_cutoff_fixednodes(...)`：用固定节点数和选定边界计算热项。
- `Omega_minimal(...)` 与 `dOmega_*`：把热项接入最小 PNJL Omega，并用 ForwardDiff 计算导数。

最值得保留的设计不是具体常数，而是边界：

- 截断选择属于控制逻辑，应在 AD 图之外完成。
- 被积函数本身仍保留 Dual 传播，避免把整个热项退化成纯 Float64。
- 离散 `pmax` 边界让缓存命中、结果可审计，也减少连续截断导致的导数噪声。

### 双区费米面聚焦积分

`integral_test/Pnjl_pure.jl` 和 `integral_test/test_focus_transform.jl` 中存在另一条积分策略：

- 当 `mu > mass` 时计算费米动量 `pF = sqrt(mu^2 - mass^2)`。
- 约 40% 节点映射到 `[0, pF)`，并用幂律把节点聚焦到费米面附近。
- 剩余节点用 Mobius 型变换覆盖 `[pF, infinity)`。
- `beta = max(delta_mev / hc, 5T, 0.5)` 控制费米面外侧伸展尺度。

tracked dirty diff 的核心是对这些动态节点和权重追加 `ForwardDiff.value.(...)`：

- `_linear_transform_eval(...)` 冻结 `p_transformed` 与 `weights`。
- `thermal_focus_transform(...)` 冻结 `p_low/weights_low` 和 `p_high/weights_high`。
- `calculate_thermal_term_cutoff(...)` 冻结 `p_mapped/w_scaled`。
- `test_focus_transform.jl` 补回 `if abspath(PROGRAM_FILE) == @__FILE__` 脚本保护。

这说明旧实验已经识别到一个关键风险：费米面位置和变换尺度依赖参数，如果直接参与 AD，导数将包含“节点移动”贡献。是否应包含该贡献是一个数值方法选择，而不是普通性能优化。

### `p_peak` / `p_sup` 估计与贝叶斯调参

`scripts/verify_analytic_ppeak_pmax.jl` 提供数值扫描：

- 对 `p^2 * [log_AA(E-mu) + log_AAbar(E+mu)]` 扫描峰值。
- 以峰值的相对阈值寻找尾部有效上界 `p_sup`。
- 对比解析或半解析 `p_peak` / `p_max` 估计。
- 输出 `verify_ppeak_pmax.csv` 和散点图。

`scripts/bayesopt_estimate_pmax.jl` 在此基础上用 `GaussianProcesses` / `Distributions` 做简易 Bayesian optimization，调参对象是：

- `log10(eps)`
- `alpha`
- `beta`
- `pmin_mev`

目标函数同时惩罚相对误差和低估 `p_sup` 的情况。这个设计可以保留为“离线调参 harness”思路，但不适合直接进入主线：生产代码不应在运行时依赖 GP 调参，也不应在脚本开头自动 `Pkg.add`。

### 冻结策略审计脚本

`scripts/bench_node_freeze_vs_fixed.jl` 对比三种口径：

- 外部预计算固定节点。
- 函数内部根据参数计算节点，但对节点和权重做 `ForwardDiff.value` 冻结。
- 不冻结节点和权重。

`scripts/regression_pnjl_freeze_vs_nofreeze.jl` 通过机械替换生成 `Pnjl_pure_nofreeze.jl`，比较 freeze / no-freeze 的值、导数和耗时，并输出摘要 CSV。

`tests/` 下的 untracked 测试覆盖了：

- `analytic_p_max` 是否返回预定义节点或低温快捷路径。
- 节点缓存幂等性。
- 热项积分边界条件。
- Omega 组成项。
- AD 导数与有限差分的一致性，典型容差是 5%。

## 结果与实验文件

旧工作区保存了大量输出：

- `thermal_integrand_summary.csv`
- 多组 `thermal_integrand_T*_muB*_flavor*.csv`
- `benchmark_performance.csv`
- `compare_minimal_vs_pure.csv`
- `extreme_params_omega_compare.csv`
- `integrand_bounds.csv`
- `p_max_heatmap.png`
- `p_sup_heatmap.png`
- `verify_pmax_scatter.png`
- `verify_ppeak_scatter.png`

这些结果适合做历史实验参考，但不应直接进入当前主项目 baseline。原因是它们缺少当前 `Julia_RelaxTime` 的 commit provenance、模型 profile、参数指纹、脚本合约和回归准入说明。

## 2026-08-16 复核结果

本次从完整归档重新读取了旧实现、dirty diff、CSV 摘要，并用当前主项目实现做了 focused verification。

旧归档结果的可复核摘要：

- `compare_minimal_vs_pure.csv` 有 50 个常规样本；`rel_cut_qg` 最大约 `4.33e-4`，但 `dmuB_rel_min_qg` 最大约 `3.71e-2`。
- `extreme_params_omega_compare.csv` 有 110 个极端样本；`rel_min_qg` 最大约 `0.263`，异常样本的 `dPhi1` 相对误差最高约 `54.8` 倍，集中在 `T=1e-8 MeV` 一带。
- `benchmark_performance.csv` 有 100 个样本；记录的 `speedup_omega` 范围约为 `0.697--4.284`，`speedup_grad` 范围约为 `0.327--1.925`，因此不能据此声称冻结策略对梯度有稳定性能收益。

旧方案的边界问题：

- `analytic_p_max` 在 `T <= 0` 或 `mu <= 0` 时直接返回 `0.0`；这不是一般有限温度热尾部的通用处理。
- 低温快捷路径只取 `1.3 * pF`，且离散边界策略、阈值和经验参数没有当前主线所需的 profile、误差预算与 solver 级 gate。
- Bayesian 脚本固定粗扫描和随机样本，并在脚本开头尝试 `Pkg.add`；其输出只能是离线探索，不能成为生产配置或 baseline。
- `ForwardDiff.value` 冻结节点选择确实是有价值的数值方法提示，但旧实现只证明了旧最小模型上的有限样本行为，不能直接推导当前各向异性、磁场和介子 workflow 的结论。

当前主项目复核通过：

- `src/models/pnjl_physics/PNJLIntegrals.jl` 的 `NODE_CACHE` 键已经包含 `(p_num, t_num, p_max_inv_fm)`；不同热项上限不会复用同一组节点。
- `PNJLModel` 已有 `:tensor_gauss` 与 `:rs_reduced_adaptive` 两条策略；后者使用仓库内自适应 Gauss-Legendre、费米动量断点和 compactified infinite tail，并让节点选择使用 primal 值而不污染 AD。
- 当前 unit 测试 `PNJLIntegrals`：`31/31` 通过；当前 PNJL integral ForwardDiff/integration smoke：`2/2` 与 `3/3` 通过；phase thermal quadrature validation：`66/66` 通过。
- 额外检查 `(p_num, t_num, p_max_inv_fm)` 的不同键、节点对象区分和热项 AD 导数：通过。

结论：旧归档没有需要立即移植的生产代码、参数或数值结果。可迁移内容已经沉淀为本文件中的方法边界；当前主项目已有更完整的实现和验证，不再启动独立的 `pmax`/节点冻结移植任务。

## 当前主项目对应能力

`Julia_RelaxTime` 当前已经吸收了其中的通用方法要求，但没有复制旧 workflow 的脚本和经验参数：

- `src/models/pnjl_physics/PNJLIntegrals.jl` / `src/models/pnjl_physics/PNJLCore.jl`：提供当前 models 侧 PNJL 热项积分；`src/models/pnjl/core/` 仅保留兼容性路径，不能作为当前生产入口。
- `src/Constants_PNJL.jl`：已有 `thermal_p_max_inv_fm` 配置项，支持把热项上限从常量层下沉到 profile/config。
- `src/models/pnjl_physics/PNJLModel.jl`：已提供 `:tensor_gauss` 和 `:rs_reduced_adaptive`，包含零温固定态、费米断点和无限尾部处理。
- `tests/unit/pnjl/test_pnjl_integrals.jl`、`tests/integration/models/test_pnjl_integrals_forwarddiff_smoke.jl` 与 `tests/validation/pnjl/test_phase_thermal_quadrature_validation.jl`：覆盖缓存、AD、误差和极端温度/化学势场景。
- `docs/dev/archived/2026-01-20_MesonMass_A_Integral_Cutoff20_Node16.md`：记录过 A 热项从 `pmax=10` 到 `pmax=20` 的收敛治理。
- `docs/dev/archived/2026-05-04_Friesen2019曲线验证口径说明.md`：记录过 Friesen 2019 公式层的 `∫_0^∞ dp` 热项口径，以及当前用有限上限近似时的稳定性限制。

因此，旧 workflow 不应按文件迁入。若以后需要研究动态热项截断，应围绕当前 `Models` / `src/models/pnjl_physics/PNJLIntegrals.jl` 的稳定入口另立 analysis 任务，而不是导入 `Julia_test` 的 include 驱动实验脚本。

## 不建议直接迁移的原因

- `Julia_test` 的 README 仍是早期“神经网络学习简单积分”说明，已经明显落后于当前 dirty 内容。
- 关键实现大量位于 untracked 文件中，未形成仓库内稳定提交历史。
- `scripts/bayesopt_estimate_pmax.jl` 会在脚本中检查并尝试安装包，不符合主项目脚本治理。
- `integral_test/Pnjl_minimal.jl` 里有大量设计说明和测试需求编号，但这些需求并没有进入 `Julia_RelaxTime` 的 `docs/dev/active/` 或测试分层。
- 离散 `pmax` 与低温快捷路径的物理误差、导数误差和 solver 稳定性尚未在当前主项目 profile 下验证。
- 双区费米面聚焦方法与当前 PNJL 各向异性热项、rotation 热项、磁场热项的接口边界尚未定义。
- `outputs/` 中的大型 CSV/PNG 是实验产物，不应直接作为主项目 regression target。

## 未来重做建议

如果以后要把这套能力主线化，建议按以下顺序新建任务：

1. 先写 `docs/dev/active/` 任务单，标题可设为“PNJL 热项积分自适应截断与 AD 稳定性治理”。
2. 在主线中明确目标口径：
   - 固定大上限 `thermal_p_max_inv_fm`。
   - 离散 `pmax` cache 边界。
   - 费米面双区聚焦。
   - 或仅作为 analysis/probe，不改变生产默认。
3. 先实现独立、可测试的小模块，避免直接改 solver：
   - `estimate_thermal_pmax(...)`
   - `thermal_cutoff_policy(...)`
   - `thermal_nodes_for_policy(...)`
   - `calculate_log_sum_with_cutoff_policy(...)`
4. 单位和输入契约必须显式：
   - `mass/T/mu` 使用 `fm^-1`。
   - 外部展示的边界使用 MeV。
   - `mu` 是 quark chemical potential 还是 `mu_B / 3`。
   - `Phi1/Phi2` 的允许域与越界保护。
5. 验证应分层：
   - unit：`pmax` 边界、缓存幂等性、边界条件。
   - integration：PNJL 单点 Omega 和 gap residual 不崩溃。
   - regression：固定点值与导数不漂移。
   - validation/analysis：对 `∫_0^∞ dp` 参考或高上限参考做误差审计。
6. 冻结策略必须明确记录：
   - `ForwardDiff.value` 冻结的是节点选择 / 节点位置 / 权重。
   - 不冻结的是热项被积函数对 `mass/T/mu/Phi` 的依赖。
   - 允许的导数误差必须用有限差分或更高精度参考证明，而不是只凭运行成功。
7. 如果保留 Bayesian tuning，只应作为离线 analysis script：
   - 固定随机种子。
   - 输出参数、样本范围、objective 版本和失败统计。
   - 不作为生产路径依赖。

## 处理建议

对当前隔离区 `Julia_test`：

- 不迁移旧代码、旧配置或旧输出；本次没有新增主项目生产文件。
- 保留本文件作为主项目侧的长期方法审计记录；完整归档仅作为本次清理前的短期可追溯证据。
- 完整归档清理前已验证可读，共 621 个 tar 条目；已核对 SHA-256 与审计解包内容。
- 本项已完成彻底清理：删除 `Julia_test` 原目录、完整归档、文件清单、review 解包目录和轻量 `_safety_exports\Julia_test`，不清空 Windows 回收站，也不影响其他隔离项目或主项目现有 dirty 文件。
