# Models 扫描入口总览

本页回答三个问题：

- 我应该用 `Models.run_tmu_scan` 还是 `Models.run_trho_scan`？
- 什么时候需要 `Models.build_default_rho_grid`？
- 扫描结果会进入哪些下游链路？

## 首选入口

对大多数调用方，扫描主题的首选入口只有两个：

- `Models.run_tmu_scan`：在 `(T, μ, ξ)` 网格上扫描，适合热力学趋势、边界线预扫、与固定化学势链路对接
- `Models.run_trho_scan`：在 `(T, ρ, ξ)` 网格上扫描，适合相结构、S 形检测、Maxwell 构造、CEP 与 phase pipeline 前置数据生成

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