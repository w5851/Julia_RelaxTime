# 求解失败问题记录 (2026-03-05)

## 概述

本次 Plan A/B 补算中遇到的 ζ (bulk viscosity) 求解失败问题。η/s 和 σ/T **全部有效**（0 NaN），问题仅影响 ζ/s。

## 问题分类

### 1. `bulk_coeffs_isentropic.masses` 含负值

- **位置**: `src/relaxtime/TransportCoefficients.jl:700` (`_validate_transport_inputs`)
- **原因**: `bulk_viscosity_coefficients`（`ThermoDerivatives.jl`）内部使用 **ForwardDiff 自动微分 + ImplicitFunction（ImplicitDifferentiation.jl）隐函数微分** 计算 dM/dT、dM/dμ_B。隐函数微分通过 `ImplicitFunction(forward_solve_impl, conditions_impl)` 包装 `solve_gap`，再用 `ForwardDiff.jacobian(solve_state, θ)` 对 gap 方程解求导。在高温区 AD chain 产生的质量值可能略微为负（非有限差分问题，而是 AD + NLsolve 精度限制所致）。
- **触发条件**: T ≥ 215 MeV（ξ ≠ 0 时更早），μ_B = 0
- **影响**: Plan A 中占绝大多数 ζ 相关 warning 点
- **已修复 (2026-03-05)**:
  1. `_validate_transport_inputs` 中 mass ≥ 0 检查从 `error()` 改为 `@warn`
  2. `bulk_viscosity_isentropic` 中 `masses = abs.(bulk_coeffs_isentropic.masses)` 钳制负值
  3. `TransportWorkflow.solve_transport_from_equilibrium` 中 `bulk_viscosity_coefficients(...)` 调用加 try/catch，失败则 ζ=NaN
- **注意**: 相变临界温度附近（T~215 MeV），abs() 钳制可能导致 ζ 数值尖峰（artifact），后续需人工检查/标记

### 2. `solve_gap did not converge` (AD chain 内部)

- **位置**: `src/models/gap_solver.jl:359`
- **原因**: ImplicitFunction 隐函数微分需要在扰动参数 θ=(T, μ) 处重解 gap 方程。ForwardDiff 对 solve_gap 求 Jacobian 时，其内部 NLsolve 可能不收敛（**不是**因为有限差分偏移 T±δT，而是 AD dual number 运算路径上的数值问题）。
- **触发条件**: T ~ 230-270 MeV，ξ = -0.3，μ_B = 0；`residual_norm` 在 1.6-14.1 范围
- **影响**: Plan A 中部分 ζ=NaN 点（abs 钳制无法解决，因为 solve_gap 本身不收敛）
- **已修复 (部分)**: TransportWorkflow 层 try/catch 保证不 crash，但 ζ 仍为 NaN
- **后续方向**: 改善 AD chain 内 NLsolve 收敛性（增大迭代上限、改善初始猜测）或考虑解析 Jacobian 替代隐函数微分

## Plan A: ζ 有效范围统计

| ξ | 首个 NaN T (MeV) | ζ 有效点数/总 61 | 最高有效 T |
|---|---|---|---|
| -0.3 | 215 | 28 | 355 |
| -0.2 | 215 | 30 | 375 |
| -0.1 | 220 | 27 | 340 |
| 0.0 | 220 | 31 | 360 |
| 0.1 | 220 | 28 | 400 |
| 0.2 | 220 | 25 | 355 |
| 0.3 | 225 | 28 | 390 |

> **注**: T < 215 MeV 范围（包含手征对称性破缺相和伪临界温度附近）的 ζ/s 数据完整。论文核心图表数据不受影响。

## Plan B: ζ 有效范围统计

| 温度 | ζ=NaN 数/总 51 | 说明 |
|---|---|---|
| T=150 | **0** ✅ | 全部有效 |
| T=190 | **0** ✅ | 全部有效 |
| T=200 | **0** ✅ | 全部有效 |
| T=250 | **48** ❌ | 仅 ξ=-0.48, 0.08, 0.34 三点有效 |

## 当前容错措施

- `run_gap_transport_scan.jl`: 双层 try/catch
  1. bulk 失败 → 回退到 `compute_bulk=false`，ζ=NaN，η/σ/τ 正常输出
  2. 整点失败 → skip 并 @warn，不中断后续扫描
- `scan_transport_vs_xi_T150_muB800.jl`: 同上双层容错

## 后续永久性修复方向

1. **短期**: ✅ 已完成 — mass ≥ 0 检查改为 warning + abs 钳制；try/catch 兜底
2. **中期**: 提高 AD chain 内 NLsolve 的 max iterations 和收敛判据；考虑为 ImplicitFunction 提供更好的初始猜测
3. **长期**: 考虑用解析 Jacobian（手动推导 ∂F/∂x, ∂F/∂θ）替代 ImplicitDifferentiation.jl 的自动隐函数微分，以增强收敛稳定性
