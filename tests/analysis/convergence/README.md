# 收敛性测试说明（RelaxTime）

本目录包含针对 RelaxationTime 输运计算的收敛性回归测试脚本，主要用于评估数值积分参数（动量、角度、σ 网格等）对关键物理量（`τ`、`η`、`σ`、`η/s`）的影响，并给出推荐的默认参数。

新增了 A 函数（单线积分）在 `p ∈ [0, 10]`、`p_nodes=32` 下的收敛性验证脚本，用于确认默认积分上限与节点数的可靠性。
新增了 B0 函数收敛性脚本，用更严格积分参数作为参考值，验证默认参数的稳定性。

## 覆盖范围

脚本 [tests/analysis/convergence/run_convergence.jl](tests/analysis/convergence/run_convergence.jl) 对以下参数做系统扫描：

- `tau_p_nodes`：平均散射率积分的动量节点数
- `tau_n_sigma_points`：计算 `σ(s)` 时的 `t` 积分点数
- `tau_angle_nodes` / `tau_phi_nodes`：角度/φ 积分节点数
- `sigma_grid_n`：`σ(s)` 预计算网格点数

测试点包含：
- (T=200 MeV, μB=0, ξ=0)
- (T=200 MeV, μB=800 MeV, ξ=0)
- (T=200 MeV, μB=800 MeV, ξ=-0.4)
- (T=200 MeV, μB=800 MeV, ξ=0.4)

## 运行方式

```bash
julia --project=. tests/analysis/convergence/run_convergence.jl
```

```bash
julia --project=. tests/analysis/convergence/test_A_convergence.jl
```

```bash
julia --project=. tests/analysis/convergence/test_B0_convergence.jl
```

可通过环境变量调节容忍度：
- `CONV_TOL`（默认 0.08）
- `CONV_TOL_ANGLEPHI`（默认 0.12）
- `CONV_TOL_SIGMA_GRID`（默认 0.20）

对 A 函数收敛性脚本可使用：
- `A_CONV_TOL`（默认 5e-3）

对 B0 函数收敛性脚本可使用：
- `B0_CONV_TOL`（默认 5e-4）

## 最新结论（摘要）

- `tau_p_nodes`：
  - `p=20` 相对 `p=24` 的差异通常在 0.5%~1.4% 内；
  - `p=24` 更稳健，建议作为基准或高精度配置。
- `tau_n_sigma_points`：
  - `6→8→10` 的变化非常小（多数 <0.1%），`8` 已较稳健，`10` 更保守。
- `tau_angle_nodes` / `tau_phi_nodes`：
  - 由 `4/8` 提升到 `6/12` 时，在部分点（尤其 `μB=0` 或 `ξ≠0`）差异可达约 4%~11%；
  - 建议 `6/12` 作为更可靠的默认或用于精度验证。
- `sigma_grid_n`：
  - `40` 显著偏低（误差可达 ~8%~20%），`60` 仍可能有 2%~8% 波动；
  - 建议 `80` 作为较稳健的选择，`60` 可用于更快的扫描。

## 推荐默认参数（若追求稳健性）

- `tau_p_nodes = 24`
- `tau_n_sigma_points = 8`（或 10）
- `tau_angle_nodes = 6`
- `tau_phi_nodes = 12`
- `sigma_grid_n = 80`

如需更高性能，可将 `tau_p_nodes=20`、`sigma_grid_n=60` 作为折中方案，但建议在目标区域做一次快速收敛对照。
