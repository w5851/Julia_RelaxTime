# PNJL 涨落计算 AD 接入路径（含隐式求导边界）

本文档用于落实开发待办中的“梳理涨落计算自动微分接入路径”，目标是把现有能力与边界整理成可执行步骤，避免在涨落量实现时重复走弯路。

## 1) 当前可复用 AD 能力

### 1.1 隐函数求解与 `dx/dθ`

- 入口：`PNJL.create_implicit_solver`（`src/pnjl/solver/ImplicitSolver.jl`）
- 形态：`ImplicitFunction(forward_solve_mu, conditions_mu)`
- 参数：`θ = [T_fm, μ_fm]`
- 状态：`x = [φ_u, φ_d, φ_s, Φ, Φ̄]`
- 现状：可直接通过 `ForwardDiff.jacobian(θ -> solver(θ)[1], θ)` 取得 `dx/dθ`

### 1.2 导数量聚合（已在主线使用）

- `PNJL.ThermoDerivatives.thermo_derivatives`
  - 基于隐式求解 + ForwardDiff 产出：
    - `dP_dT`, `dP_dmu`, `dEpsilon_dT`, `dEpsilon_dmu`, `dn_dT`, `dn_dmu`
    - 组合导数：`dP_depsilon_n`, `dP_dn_epsilon`
- `PNJL.ThermoDerivatives.bulk_viscosity_coefficients`
  - 使用 `dx/dθ` + 链式法则构造 `ds/dθ`, `dn/dθ`, `dM/dθ`

### 1.3 现有验证覆盖

- `tests/unit/pnjl/test_implicit_jacobian.jl`
  - 已验证隐函数定理计算的 `dx/dθ` 与手工公式（`-(∂F/∂x)^{-1}∂F/∂θ`）一致。
- `tests/unit/pnjl/test_thermo_derivatives.jl`
  - 已验证热力学导数接口可用与基础一致性。

---

## 2) 隐式求导边界（必须显式约束）

以下场景必须在涨落计算中做防护或降级：

1. **雅可比近奇异边界**
   - 条件：`∂F/∂x` 条件数极高或不可逆。
   - 影响：`dx/dθ` 放大噪声，涨落量可能爆发式抖动。

2. **多分支/相变切换点**
   - 条件：解分支切换（特别是一阶相变附近）。
   - 影响：局部可导性变差，单分支导数不再代表物理解的全局趋势。

3. **非光滑算子穿越 AD 边界**
   - 条件：`if/clip/min/max/argmin` 等离散分支选择直接参与导数链路。
   - 影响：导数不连续或定义失真。

4. **分母退化**
   - 条件：组合导数分母接近 0（如 `dP_depsilon_n`、`dP_dn_epsilon` 的分母）。
   - 影响：返回 `NaN/Inf`，必须转为“不可用点”语义处理。

---

## 3) 接入策略（建议主线）

### 3.1 分层实现

- **第 0 层：状态导数层**
  - 输入：`(T_fm, μ_fm, xi, p_num, t_num)`
  - 输出：`x`, `dx_dT`, `dx_dμ`，以及可选 `cond_dFdx`。
- **第 1 层：基础热力学导数层**
  - 复用 `thermo_derivatives` 已有返回。
- **第 2 层：涨落组合层（新增）**
  - 仅做组合与后处理，不重复求解隐函数。
  - 输出应标注每个量的 `status`（`ok/singular/non_smooth`）。

### 3.2 边界守卫（最低要求）

- 对 `∂F/∂x` 条件数设置阈值（例如 `κ > 1e10` 记为 `singular`）。
- 对关键分母设置 `abs(denom) < ε` 防护（建议 `ε = 1e-10` 起步，可配置）。
- 保留 `converged/iterations/residual_norm` 并随涨落输出一并落盘。
- 对分支切换点输出标记，避免“静默成功”。

### 3.3 最小回归策略

- 固定点回归：
  - 低 μ 平滑区 1 点；
  - 近相变区 1 点（允许标记为 `singular/non_smooth`，但不能崩溃）。
- 验证目标：
  - 可导区数值稳定（rtol/atol）；
  - 非可导区行为可解释（状态码一致）。

---

## 4) 最小 PoC 任务清单

1. 在 `PNJL.ThermoDerivatives` 增加“涨落导数聚合”入口（不改旧 API，新增函数）。
2. 返回结构统一为：`values + diagnostics + status`。
3. 将该入口接到扫描链路（先 T-μ，后 T-ρ）。
4. 为新增入口补两类单测：
   - 平滑区一致性；
   - 奇异/分支边界容错与状态码。

---

## 5) 与现有文档关系

- 扫描输出字段契约：`docs/api/pnjl/ScanOutputContract.md`
- 热力学导数 API：`docs/api/pnjl/ThermoDerivatives.md`
- 顶层 PNJL 导航：`docs/api/pnjl/PNJL.md`

本文件只定义“如何接入 AD 与边界治理”，不替代具体公式文档。