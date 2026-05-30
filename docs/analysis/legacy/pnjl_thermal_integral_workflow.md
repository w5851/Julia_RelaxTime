# PNJL Thermal Integral Legacy Workflow Notes

本文档记录隔离区 `Julia_test` 仓库中的 PNJL 热项积分实验 workflow，供以后判断是否值得在 `Julia_RelaxTime` 主线中重做。本文只做能力沉淀与迁移边界说明，不表示这些旧脚本已经成为当前主项目的有效入口。

## 来源与状态

- 隔离区来源：`D:\Desktop\_cleanup_quarantine\2026-05-27\Julia_test`
- 远端来源：`https://github.com/w5851/Julia_test.git`
- 当前状态：`main...origin/main` 无 ahead，但工作区 dirty。
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
- 安全导出：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\Julia_test`

注意：当前 `_safety_exports` 中有 git bundle、tracked diff、status 和 untracked 文件列表，但不包含 untracked 文件正文。因此在决定删除完整目录前，必须先确认是否还需要保留 `integral_test/`、`scripts/`、`tests/` 和 `outputs/` 下的 untracked 内容。

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

## 当前主项目对应能力

`Julia_RelaxTime` 当前已经有若干相关能力，但没有完整吸收这套旧 workflow：

- `src/models/pnjl/core/PNJLIntegrals.jl`：提供 `cached_nodes` 与 `calculate_log_sum`，用于 models 侧 PNJL 热项积分；当前仍是固定节点和固定热项区间口径。
- `src/Constants_PNJL.jl`：已有 `thermal_p_max_inv_fm` 配置项，支持把热项上限从常量层下沉到 profile/config。
- `docs/dev/archived/2026-01-20_MesonMass_A_Integral_Cutoff20_Node16.md`：记录过 A 热项从 `pmax=10` 到 `pmax=20` 的收敛治理。
- `docs/dev/archived/2026-05-04_Friesen2019曲线验证口径说明.md`：记录过 Friesen 2019 公式层的 `∫_0^∞ dp` 热项口径，以及当前用有限上限近似时的稳定性限制。

因此，旧 workflow 不应按文件迁入。若以后要恢复，应围绕当前 `Models` / `src/models/pnjl/core/PNJLIntegrals.jl` 的稳定入口重做，而不是导入 `Julia_test` 的 include 驱动实验脚本。

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

- 不建议现在删除完整目录。
- 建议标记为“中高保留价值实验参考”，直到确认 untracked 文件已经以文档、归档包或其他方式保留。
- 当前最值得后续吸收的是设计思路：AD 节点冻结、离散 `pmax` cache、费米面聚焦、`pmax` 参数离线调参 harness。
- 当前不值得直接吸收的是脚本组织、自动安装依赖、旧 README、实验输出和 include 驱动入口。
- 若以后磁盘清理需要删除，应先额外打包或复制 untracked 内容；仅依赖 `_safety_exports\Julia_test` 不足以恢复这些脚本和结果。
