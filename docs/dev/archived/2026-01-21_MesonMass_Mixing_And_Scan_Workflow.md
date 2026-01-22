---
title: 完善介子质量计算功能（混合介子 + 扫描工作流）
archived: true
original: docs/dev/active/2026_1_21完善介子质量计算功能.md
archived_date: 2026-01-21
---

以下为原始内容（保留，以便审阅与历史参考）：

***

当前只实现了对 π/K（以及 σ_π、σ_K）介子的介子质量计算与 Mott 相变相关的计算；需要补充“混合介子”的计算能力。

本计划阶段只基于现有工程实现进行完善（先不写 Fortran 对齐目标）。

## 当前状态（2026-01-21）

核心目标已完成并通过单元测试：

- MesonMass 支持混合介子：`:eta`, `:eta_prime`, `:sigma`, `:sigma_prime`
- MottTransition 增加混合介子阈值/差值接口：`mott_threshold_masses`、`mott_gaps`（保留旧接口兼容）
- 混合矩阵元素计算已收敛为单一来源：`EffectiveCouplings.mixing_matrix_elements`，并由 MesonMass/MesonPropagator 复用，避免公式漂移
- 文档已同步：`docs/api/*` 与公式文档中的 M08 符号约定

已补齐 Step 4（扫描脚本/工作流串联输出）：

- 工作流：src/pnjl/workflows/MesonMassWorkflow.jl
- 扫描脚本：scripts/relaxtime/run_gap_meson_mass_scan.jl

## 现状盘点（基于当前代码）

- 介子质量求解入口在 [src/relaxtime/MesonMass.jl](src/relaxtime/MesonMass.jl)：
	- 目前支持 `:pi`, `:K`, `:sigma_pi`, `:sigma_K`。
	- 使用极点方程 `f = 1 - 4K * Π(p0)`，并用 NLsolve 解 `Re(f)=0, Im(f)=0` 得到 `(M, Γ)`。
- 极化函数在 [src/relaxtime/PolarizationAniso.jl](src/relaxtime/PolarizationAniso.jl)：
	- 已支持 `channel=:P/:S`，以及宽度输入（通过 `polarization_with_width`）。
	- 目前一次调用对应某一对夸克味道（由 `(m1, m2, μ1, μ2, A1, A2)` 指定）。
- 混合通道所需的有效耦合系数已经具备：
	- [src/relaxtime/EffectiveCouplings.jl](src/relaxtime/EffectiveCouplings.jl) 已提供 `K0_±, K8_±, K08_±` 与 `det_K_±`。
- 混合介子的公式在 [docs/reference/formula/relaxtime/MesonMass.md](docs/reference/formula/relaxtime/MesonMass.md) 已给出（η/η′ 与 σ/σ′）。

## 目标范围（建议的“混合介子”最小可用集）

1. 赝标量混合：η / η′（P 通道，使用 `K0_plus/K8_plus/K08_plus`）。
2. 标量混合：σ / σ′（S 通道，使用 `K0_minus/K8_minus/K08_minus`）。

命名建议（避免 Unicode/引号问题）：
- `:eta`, `:eta_prime`
- `:sigma`, `:sigma_prime`

## 开发步骤（基于现有模块最小侵入式扩展）

### Step 1：在 MesonMass 中增加“混合介子”的方程构造

目标：新增一条与现有 `meson_mass_equation(::Symbol, ...)` 同风格的路径，但内部基于 2×2 混合结构构造一个标量“逆传播子”并返回复数残差。

建议实现方式：

1) 复用现有极化函数，先分别计算（同一 `p0`、同一 `k_norm`）：
- `Π_uu`：`polarization_with_width(channel, k0, gamma, k_norm, m_u, m_u, μ_u, μ_u, ...)`
- `Π_ss`：`polarization_with_width(channel, k0, gamma, k_norm, m_s, m_s, μ_s, μ_s, ...)`

2) 用 EffectiveCouplings 里已有的 `K0/K8/K08/det_K` 按 [docs/reference/formula/relaxtime/MesonMass.md](docs/reference/formula/relaxtime/MesonMass.md) 组装 `M00/M08/M88`。

3) 由 `(M00, M08, M88)` 构造两个“本征逆传播子”（η 与 η′，或 σ 与 σ′）：
- `Minv_light = M00 + M88 - sqrt((M00 - M88)^2 + 4*M08^2)`
- `Minv_heavy = M00 + M88 + sqrt((M00 - M88)^2 + 4*M08^2)`

4) 在 `meson_mass_equation` 中：
- 对 `:eta` 返回 `Minv_light(p0)`
- 对 `:eta_prime` 返回 `Minv_heavy(p0)`
- 对 `:sigma`/`:sigma_prime` 同理，但使用 S 通道的耦合与极化

### Step 2：默认初值与求解稳定性

目标：在 `default_meson_mass_guess` 中补齐混合介子初值，保证 NLsolve 不至于跑偏。

建议：
- `:eta` 初值靠近 `2*m_u` 或一个常量（比如按真空物理质量 547 MeV 换算到 fm⁻¹）；
- `:eta_prime` 初值靠近 `2*m_s` 或 958 MeV；
- `:sigma`、`:sigma_prime` 初值可用 `2*m_u` 与 `2*m_s`。

并在计划实现时补一个“重试策略”（可选）：若第一次 NLsolve 不收敛，则尝试用更小 `initial_gamma` 或调整初值到阈值附近。

### Step 3：MottTransition 对混合介子的定义

当前 [src/relaxtime/MottTransition.jl](src/relaxtime/MottTransition.jl) 用“组分夸克阈值”定义 Mott 点。

混合介子存在两个潜在阈值（uu 与 ss）。需要在实现前明确一个可执行定义：

- 选项 A（保守、便于用）：用更低阈值作为分解起点
	- `threshold(:eta/eta_prime) = min(2m_u, 2m_s)`
	- 标量同理
- 选项 B（更细、信息更完整）：同时输出两个 gap
	- `gap_uu = M - 2m_u`, `gap_ss = M - 2m_s`
	- Mott 判据可按任一 gap 过零来判定，并在输出中标注是哪一个通道先开口

建议先做 B（更不容易“拍脑袋”），但需要扩展 API/输出格式。

### Step 4：脚本/工作流层的串联输出（后续可选，但建议在计划里写清）

现有对比脚本（只算 π/K）已从 `scripts/relaxtime` 清理；代码集中归档在 [docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md](docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md)。

建议新增一个纯 Julia 的扫描脚本（或一个 workflow 模块）来统一输出：
- π/K/η/η′/σ/σ′ 的 (M, Γ)
- 以及 Mott gap（按 Step 3 的定义）

可参考已有工作流结构 [src/pnjl/workflows/TransportWorkflow.jl](src/pnjl/workflows/TransportWorkflow.jl)，但是否放入 `pnjl/workflows` 还是 `relaxtime` 下，需要结合“是否依赖 PNJL 求解平衡解”来定：
- 若脚本/模块需要先跑 gap 方程得到 `m_u,m_s,Φ,Φbar`，更适合放在 `src/pnjl/workflows`。
- 若仅消费外部给定的 `quark_params/thermo_params`，放在 `src/relaxtime` 更清爽。

## 预计会遇到的问题（提前标注风险点）

1. 归一化/符号一致性：
	 - 当前 π/K 使用 `1 - 4K Π = 0`；混合通道的 `M00/M08/M88` 公式是否与工程里 `Π` 的定义完全同归一化，需要实现时做一次“退化检查”（例如 `G_u==G_s` 时 `K08→0`，混合应消失）。
2. 复数平方根分支：
	 - `sqrt` 在复数上有分支切换，可能导致 η/η′ 在某些参数点发生“标签翻转”（本征值交换）。
	 - 需要一个稳定的判别策略：例如按连续性跟踪、或按与初值更接近的根来选择“η vs η′”。
3. Mott 判据在混合介子上不唯一：
	 - 必须在 API 层明确输出哪个阈值/哪个 gap，避免后续分析歧义。
4. 数值敏感区：
	 - 混合通道依赖 `det_K` 与两套极化的组合，接近 `det_K→0` 或极点靠近阈值时会非常敏感；需要在输出里记录 `det_K_±` 以便诊断。

## 验证建议（计划阶段先写清，后续实现时落实）

1. 退化极限检查：设置 `G_u == G_s`（或构造一个近似 SU(3) 对称点），验证 `K08≈0` 时混合消失。
2. 连续性检查：沿 T 扫描，η/η′ 质量曲线不应出现明显的非物理跳变（除非确有 Mott/共振转换）。
3. 与公式文档逐项对齐：对比 [docs/reference/formula/relaxtime/MesonMass.md](docs/reference/formula/relaxtime/MesonMass.md) 中的 `M00/M08/M88` 构造，确保变量含义与工程 `Π` 的定义一致。

***
