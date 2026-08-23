# `Models.run_magnetic_scan`

这是外磁场 PNJL 的生产级 `(T, mu, eB)` 扫描入口。它与普通 `run_tmu_scan`、
`run_trho_scan` 分离，因为非零 `eB` 必须调用 magnetic Omega 的完整五维
`solve_magnetic_gap`，不能落入共享 `ProblemSpec` 约束链。

## 推荐调用

```julia
using .Models

result = Models.run_magnetic_scan(
    T_values=[50.0, 150.0],
    mu_values=[0.0, 60.0],
    eB_values=[0.0, 2.0e4],
    xi_values=[0.0],
    output_path="data/outputs/results/pnjl/scan/magnetic/selected.csv",
    candidates_output_path="data/outputs/results/pnjl/scan/magnetic/candidates.csv",
)
```

等价 CLI 为：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 `
  scripts/models/run_unified_scan.jl scan magnetic `
  --model_kind=PNJLMagnetic --solver_mode=fixed_mu `
  --T_values=50,150 --mu_values=0,60 --eB_values=0,20000 --xi_values=0 `
  --output_path=data/outputs/results/pnjl/scan/magnetic/selected.csv `
  --overwrite=true --resume=false
```

## 输入合同

- `T_values`：MeV，必须全部 `>0`。
- `mu_values`：MeV 的共同夸克化学势；adapter 将其转换为 `mu_u=mu_d=mu_s`。
  这里不是 `mu_B` 字段。
- `eB_values`：MeV^2；内部转换为 `eB_fm2=eB_MeV2/hbarc^2`。
- `xi_values`：目前只能为 `[0.0]`。
- `model_kind`：只能是 `:PNJLMagnetic`。
- `solver_mode`：只能是 `:fixed_mu`。
- `p_num`、`t_num`、`pz_max`、`n_max`、`cutoff_N` 和 solver 控制项会原样传给
  magnetic solver；默认值是保守生产配置，但不自动构成全域收敛证明。

每个点调用 `solve_magnetic_gap`。扫描沿每个 `(T,eB)` 的 `mu` 顺序传递上一点
选中态作为 continuation seed，同时保留 magnetic solver 的默认多 seed 候选。
`classify_stability=true` 只增加诊断标签，不把 Hessian 正定性提升为默认过滤条件。

## 输出合同

`output_path` 写入每个参数点一行的选中态，包含：

- `Omega/pressure/rho/entropy/energy`；
- 三味净密度、五个平均场量和三味质量；
- `selected_candidate_index`、`branch_label`、`attempt_count`、`failed_attempts`、
  `residual_norm`、`n_max`、`converged`。

`candidates_output_path` 写入每个点的全部去重后收敛候选，包含 seed、分支标签、
稳定性诊断标签、方法、迭代数、Omega、状态、残差和 `n_max`。没有候选的点仍会在
selected CSV 中写入 `converged=false` 和错误信息；它不会伪造候选分支。

返回值至少包含：

```text
total, success, failure, skipped,
output, selected_output, candidates_output
```

`resume=true` 按 `(T,mu,eB,xi)` 四元键跳过已写入的 selected 行；候选 CSV 是同一次
运行的伴随 artifact，不应单独解释为完整分支全集覆盖证明。

## 物理与流程边界

- 这是完整五维 `FixedMu` equilibrium 路径：
  `dOmega/d(phi_u,phi_d,phi_s,Phi,PhiBar)=0`。
- `FixedRho`、普通 `T-rho`、phase/Maxwell/CEP magnetic 路径尚未实现；
  `run_tmu_scan`、`run_trho_scan` 和 phase pipeline 遇到 `:PNJLMagnetic` 会显式报错，
  不会静默使用普通 PNJL residual。
- 当前 `omega/pressure` 是固定外部磁场背景下的标量物质量，不含 Maxwell 自能，
  也不提供横向/纵向压力张量或磁化强度。
- 运行结果需要结合 `n_max` 收敛报告和代表性点审计；本入口本身不宣称全域外部
  validation 或全分支完备性。

## 追溯导航

磁场路线的固定阅读顺序是：

1. 总账：[`implemented_capabilities.md`](../../../../reference/implemented_capabilities.md#14-路线-l外磁场-pnjl)
2. API 细账：本页与 [Magnetic 主题总览](../variants/magnetic/README.md)
3. 公式审核表：[`PNJL_magnetic_core.md`](../../../../reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md#开发者审核表公式到实现的对应关系)
4. 源码入口：[`MagneticScan.jl`](../../../../src/models/scans/MagneticScan.jl)、
   [`PNJLMagneticModel.jl`](../../../../src/models/pnjl_physics/PNJLMagneticModel.jl)、
   [`MagneticGapSolver.jl`](../../../../src/models/pnjl_physics/MagneticGapSolver.jl)
5. 测试证据：`tests/unit/models/test_magnetic_scan.jl`、
   `tests/unit/pnjl/test_pnjl_magnetic_model.jl`、
   `tests/regression/pnjl/test_magnetic_fixedpoint_regression.jl`

旧入口 [`run_magnetic_point.jl`](../../../../scripts/pnjl/run_magnetic_point.jl)、
[`run_magnetic_eb_scan.jl`](../../../../scripts/pnjl/run_magnetic_eb_scan.jl) 和
[`run_magnetic_stability_scan.jl`](../../../../scripts/pnjl/run_magnetic_stability_scan.jl)
仍是固定 `x_state` 的内核、`n_max` 或稳定性诊断脚本；它们不是 equilibrium 扫描的
替代品。
