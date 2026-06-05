# Models 扫描入口总览

本页回答三个问题：

- 我应该用 `Models.run_tmu_scan` 还是 `Models.run_trho_scan`？
- 什么时候应该改用 `Models.run_freezeout_fixedmu_scan`？
- 什么时候应该直接进入 `Models.run_crossover_meson_density_scan`？
- 什么时候应该直接进入 `Models.run_freezeout_meson_density_scan`？
- 什么时候应该直接进入 `Models.run_external_path_meson_density_scan`？
- 什么时候应该直接进入 `Models.run_freezeout_meson_mass_scan` 或 `Models.run_isentropic_meson_mass_scan`？
- 什么时候需要 `Models.build_default_rho_grid`？
- 扫描结果会进入哪些下游链路？

## 首选入口

对大多数调用方，扫描主题的首选入口只有两个：

- `Models.run_tmu_scan`：在 `(T, μ, ξ)` 网格上扫描，适合热力学趋势、边界线预扫、与固定化学势链路对接
- `Models.run_trho_scan`：在 `(T, ρ, ξ)` 网格上扫描，适合相结构、S 形检测、Maxwell 构造、CEP 与 phase pipeline 前置数据生成

当输入天然是 freeze-out 路径参数，而不是独立的 `(T,\mu)` 网格时，应改用：

- `Models.run_freezeout_fixedmu_scan`：在 `\sqrt{s_{NN}}` 路径上扫描，内部先映射到 `(T,\mu_B)`，再复用 `FixedMu` 主链

当你不只需要 equilibrium 点，而是要直接产出 `n_\pi(T)`、`n_K(T)`、`K/\pi(T)` 或其 BW/BU 版本时，应进一步改用：

- `Models.run_crossover_meson_density_scan`
- `Models.run_freezeout_meson_density_scan`

当路径已经由外部数据明确给出，例如文献提取后的 phase-line 离散点，而你不希望再让本仓库内部路径生成器介入时，应改用：

- `Models.run_external_path_meson_density_scan`

当目标不是介子数密度，而是沿正式路径直接输出介子质量/宽度与阈值时，应改用：

- `Models.run_freezeout_meson_mass_scan`
- `Models.run_isentropic_meson_mass_scan`

`Models.build_default_rho_grid` 是公开导出，但它的定位是辅助入口：当你需要默认的多分辨率 `ρ` 网格，或要在相图流程前自定义低密度加密策略时再直接调用。

## 怎么选入口

### 选择 `Models.run_tmu_scan`

适用于：

- 输入自然以平均化学势 `μ` 组织
- 希望沿 `μ` 或 `T` 方向做 continuity / phase-aware 跟踪
- 需要先得到一组规则的 T-μ 网格结果，再交给后续分析脚本

关键特征：

- 支持 `use_phase_aware=true`
- 首点可配合 `bootstrap_multiseed=true` 稳定选分支
- 输出 CSV 使用固定的 T-μ 契约

### 选择 `Models.run_trho_scan`

适用于：

- 输入自然以归一化密度 `ρ/ρ₀` 组织
- 目标是相图、Maxwell、spinodal、crossover 或 CEP 前置数据
- 需要显式控制 `reverse_rho`、`seed_policy`、`constraint_mode`

关键特征：

- 默认 `reverse_rho=true`，降低低密度端连续性跟踪失败风险
- 默认 `rho_values=DEFAULT_RHO_VALUES`
- 导出的 `build_default_rho_grid` 与此入口天然配套

### 选择 `Models.run_freezeout_fixedmu_scan`

适用于：

- 输入自然以 `\sqrt{s_{NN}}` 或等价一维路径参数组织
- 需要沿 chemical freeze-out baseline 做 continuity 求解
- 需要把 freeze-out 路径结果接到 charged / meson-density workflow

关键特征：

- 不是独立 solver，而是 `path parameter -> (T,\mu_B) -> FixedMu solve`
- 支持 `profile_name` 选择 freeze-out 参数化
- 支持 `traversal` 固定 continuation 方向
- 续扫键以 `sqrt_s_NN_GeV, muB_MeV, xi` 为主

### 选择 `Models.run_crossover_meson_density_scan`

适用于：

- 路径天然来自内部 crossover locator
- 目标是 charged / neutral `K/\pi` 的 crossover-line reproduction
- 需要让 stable / strict BW / current BU / generalized BU 共用同一条内部 path shell

关键特征：

- 上游路径由 `Models.build_crossover_line` 统一生成
- continuation 契约与其他 meson-density workflow 保持一致
- 输出 CSV 会保留 `crossover_method`、`crossover_variable` 与 `T/\mu_B` 对齐字段

### 选择 `Models.run_freezeout_meson_density_scan`

适用于：

- 输入天然是 `\sqrt{s_{NN}}`
- 目标已经是介子数密度或 `K/\pi`
- 需要同时固定 freeze-out baseline 与 meson chemical profile
- 需要在有需要时显式指定 flavor-level `\mu_s` phenomenology profile
- 需要在 stable / strict BW / current BU / generalized BU 间复用同一条路径扫描壳

关键特征：

- 仍复用统一 meson workflow，而不是脚本层平行拼接
- continuation 契约来自 `MesonMassWorkflow.continuation_state`
- 当前已经能表达 charge-resolved `\mu_\pi` profile
- 当前也已能表达最小 flavor-level `\mu_s` profile

### 选择 `Models.run_external_path_meson_density_scan`

适用于：

- 路径来自文献提取、外部 CSV、人工整理点列
- 你要先固定外部 `(T,\mu_B)` 点，再比较物理量
- 你希望 continuation 仍由正式 workflow 管理，而不是在脚本层自己缓存种子

关键特征：

- 输入是离散路径点列，而不是内部 path generator
- 输出会保留 `path_source / path_case_id / path_line_style` 元数据
- 可直接复用 stable / strict BW / current BU / generalized BU 同一套物理核

### 选择 combined meson-density `--path trho_asymmetric`

适用于：

- 需要把现有 `FixedAsymmetricRho` equilibrium source 接到介子数密度后处理
- 需要在同位旋不对称条件下检查 `mu_u != mu_d`、signed `mu_pi` / `mu_K`
- 当前只做 smoke / diagnostic 小网格，不生产正式高精度产物

关键特征：

- CLI 入口为 `scripts/relaxtime/run_combined_meson_density_scan.jl --path trho_asymmetric`
- 路径参数为 `--rho-values` 或 `--rhomin/--rhomax/--rhostep`
- 约束参数为 `--asym-ud-ratio-target` 与 `--asym-s-target`
- 内部先调用 `Models.solve(model, FixedAsymmetricRho(...), T_fm)`，再调用 `Models.solve_meson_point_from_equilibrium`
- 输出会额外记录 `constraint_mode`、`rho_target`、`rho_norm`、`rho_u_fm3`、`rho_d_fm3`、`rho_s_fm3`、`rho_u_over_rho_d`、`constraint_residual_norm`、`mu_u_MeV`、`mu_d_MeV`、`mu_s_MeV`、`muB_MeV`、`muQ_MeV`、`muS_MeV`

### 选择 `Models.run_freezeout_meson_mass_scan` / `Models.run_isentropic_meson_mass_scan`

适用于：

- 目标是沿正式路径输出 meson mass / width / threshold / gap
- 不希望继续在规则网格脚本上手工拼装 continuation
- 路径天然来自 chemical freeze-out 或 fixed-`σ` 等熵线

关键特征：

- 两条入口都复用 `MesonMassWorkflow.continuation_state`
- freeze-out 路径直接消费 `(T,\mu_B)` 点列
- isentropic 路径先通过 `FixedSigma` 求出路径点，再回到统一 meson workflow
- 输出 CSV 在现有 meson-mass 字段前补齐路径元数据

## 典型工作流

### 1. 快速流程验证

- T-μ：较粗 `T_values` 与 `mu_values`
- T-ρ：较粗 `rho_values`，保持 `reverse_rho=true`
- 降低 `p_num` 与 `t_num`，但仅用于 smoke 或接口联调

### 2. 常规研究扫描

- T-μ：`use_phase_aware=true`
- T-ρ：使用 `Models.build_default_rho_grid()` 或其加密变体
- 输出后统一按本主题约定的 CSV 字段合同读取：
	- T-μ：`T_MeV, mu_MeV, xi, pressure_fm4, rho, entropy_fm3, energy_fm4, phi_u, phi_d, phi_s, Phi1, Phi2, M_u_MeV, M_d_MeV, M_s_MeV, iterations, residual_norm, converged, message`
	- T-ρ：`T_MeV, rho, xi, mu_u_MeV, mu_d_MeV, mu_s_MeV, mu_avg_MeV, mu_B_MeV, mu_Q_MeV, mu_S_MeV, pressure_fm4, entropy_fm3, energy_fm4, rho_u_fm3, rho_d_fm3, rho_s_fm3, phi_u, phi_d, phi_s, Phi1, Phi2, M_u_MeV, M_d_MeV, M_s_MeV, iterations, residual_norm, converged, message`

### 3. 相图前置数据生成

- 优先使用 `Models.run_trho_scan`
- 保留 `reverse_rho=true`
- 在低密度端或 CEP 附近按需局部加密 `ρ` 网格
- 再接 phase 主题文档中的 `run_phase_pipeline`、Maxwell 与 crossover 相关链路

## 结果去向

扫描结果不会停留在“生成 CSV”这一步。它们通常继续流向：

- [../phase/README.md](../phase/README.md)：phase 主题主入口
- [../phase/Algorithms.md](../phase/Algorithms.md)：相图产线中的扫描依赖与加密策略

读取约定补充：

- 下游应按列名而不是列位置解析 CSV
- `converged=false` 的失败点仍保留在输出中，`message` 用于区分未收敛、种子失败或近似接受等原因
- 若后续追加新列，应追加在末尾，避免破坏现有读取脚本

## 进阶入口

如果你需要直接控制网格生成逻辑，而不是接受默认 `DEFAULT_RHO_VALUES`，再使用 `Models.build_default_rho_grid`。这类调用更偏维护者或研究脚本，而不是多数用户的首选入口。
