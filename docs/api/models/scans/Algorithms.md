# 扫描主题的职责核心

本页是 scan 主题的第二层。文件名使用 `Algorithms.md`，但语义上属于统一的“职责核心层”：解释扫描顺序、采样策略、输出契约与下游流程之间的边界。

## 1. 为什么扫描主题需要独立治理

`run_tmu_scan`、`run_trho_scan` 与 `build_default_rho_grid` 并不是三个彼此独立的小函数。它们共同定义了：

- 网格如何被组织
- continuity / phase-aware 跟踪沿哪条路径推进
- 默认采样分辨率如何设置
- 下游 phase、分析与回归脚本能否稳定读取结果

因此，本主题的核心不是“再写一遍函数签名”，而是把职责边界说明清楚。

## 2. T-μ 扫描的核心职责

`Models.run_tmu_scan` 负责在 `(T, μ, ξ)` 网格上组织求解，并在需要时启用 phase-aware 跟踪。

职责重点：

- 校验 `T_values`、`mu_values`、`xi_values`
- 管理断点续扫和 CSV 追加写入
- 在 `use_phase_aware=true` 时建立 `PhaseAwareContinuitySeed`
- 首点可配合 `bootstrap_multiseed=true` 做更稳定的分支自举
- 统一输出 T-μ CSV 契约

这意味着它不是单点求解器，而是“网格组织 + 分支跟踪 + 输出合同”三件事的聚合入口。

## 3. T-ρ 扫描的核心职责

`Models.run_trho_scan` 负责在 `(T, ρ, ξ)` 网格上组织求解，并把低密度端不稳定、种子连续性、约束模式等问题集中在统一入口处理。

职责重点：

- 默认使用 `DEFAULT_RHO_VALUES`
- 默认 `reverse_rho=true`，避免从 `ρ=0` 正向推进导致 continuity 容易失稳
- 接管 `seed_policy`、`constraint_mode`、`hybrid_weighted_fallback` 等高级控制
- 输出 T-ρ CSV 契约，供相图与 Maxwell 类下游读取

和 T-μ 相比，T-ρ 更接近相图产线前置层，因此它的“网格策略”比单点求解配置更重要。

## 4. `build_default_rho_grid` 的职责边界

`Models.build_default_rho_grid` 的职责不是执行扫描，而是提供一个多分辨率默认密度网格。

默认策略包含四段：

- 全区间粗网格
- `ρ<=1.0` 的中等网格
- `ρ<=0.3` 的细网格
- `ρ<=0.15` 的超细网格

它存在的意义是：

- 把经验性的推荐分辨率从调用脚本中抽出来
- 为 `run_trho_scan` 提供默认可用的研究级网格
- 为相图链路提供可重复、可解释的默认采样起点

## 5. 扫描顺序为什么属于主题核心

扫描顺序会直接影响 continuity 或 phase-aware 跟踪的物理结果，因此它必须在主题核心页解释，而不能只埋在底层实现里。

关键规则：

- T-μ 场景中，continuity 沿外层循环定义的路径推进
- T-ρ 场景中，`reverse_rho=true` 是默认推荐，因为低密度端更容易引发跟踪失败
- 若目标是相图主线，不应同时使用过粗网格和过低积分分辨率

这些不是“实现细节”，而是 API 的使用合同。

## 6. 输出契约与下游边界

扫描 API 的返回值只是快速统计；真正跨模块传递的稳定产物是 CSV。

因此，输出契约本身就是 scan 主题的一部分，而不是附属旧页：

- T-μ 输出公共字段：`T_MeV`、`xi`、`pressure_fm4`、`entropy_fm3`、`energy_fm4`、`phi_u/phi_d/phi_s`、`Phi1/Phi2`、`M_u_MeV/M_d_MeV/M_s_MeV`、`iterations`、`residual_norm`、`converged`、`message`
- T-μ 特有字段：`mu_MeV`、`rho`
- T-ρ 特有字段：`rho`、`mu_u_MeV`、`mu_d_MeV`、`mu_s_MeV`、`mu_avg_MeV`、`mu_B_MeV`、`mu_Q_MeV`、`mu_S_MeV`、`rho_u_fm3`、`rho_d_fm3`、`rho_s_fm3`
- [../phase/Algorithms.md](../phase/Algorithms.md)：phase 产线如何消费 T-ρ 结果

兼容性规则：

- 下游读取方应按列名而不是列位置解析
- 失败点也必须写入 CSV，而不是静默丢弃
- 如需扩展 schema，优先追加列而不是重命名现有列

## 7. 推荐模板与禁用组合

### smoke 模板

- T-μ：`T = 100:20:140 MeV`，`μ = 0:60:240 MeV`，`p_num=12`，`t_num=4`，`use_phase_aware=true`
- T-ρ：`T = 80:20:160 MeV`，`ρ = 0.0:0.2:2.0`，`reverse_rho=true`

### standard 模板

- T-μ：`T = 50:10:200 MeV`，`μ = 0:10:400 MeV`，`p_num=24`，`t_num=8`，`use_phase_aware=true`
- T-ρ：温度同上，`rho_values = DEFAULT_RHO_VALUES`，`reverse_rho=true`

### phase-diagram 模板

- 优先使用 T-ρ 主扫描，再进入 Maxwell / CEP / crossover
- 基线可从 `T = 30:10:350 MeV` 与 `ρ = 0.0:0.05:4.0` 起步
- 若问题集中在某一温度区段，应优先局部降低 `Δρ`，而不是盲目扩大全局网格

### 禁用或不推荐组合

- T-ρ 含 `ρ=0` 且 `reverse_rho=false`
- 一阶相变敏感区关闭 `use_phase_aware`
- 相图主线使用 `Δρ > 0.1` 的过粗均匀网格
- 同时降低积分分辨率与网格分辨率后直接下物理论断

边界划分如下：

- `Models.run_*scan`：生成扫描结果
- `Models.build_default_rho_grid`：提供默认多分辨率采样
- phase 主题：消费扫描结果，继续执行 Maxwell/crossover/CEP 等流程
- `docs/api/pnjl/*` 旧页：仅保留迁移说明与历史兼容定位