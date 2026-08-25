# `pnjl_mag` / Julia 磁场 profile 收敛与重求解诊断 v1

本目录记录三个位点、两种参数 profile 和有限积分节点档位的定向诊断。它用于区分
外部跨求解器复核口径与 Julia 生产口径，不创建正式 validation target，也不代表全域
磁场扫描已经完成收敛。

## 两种 profile

`source-parity` 是跨求解器 diagnostic profile：严格使用外部仓库约定的
`hc=197.33 MeV fm`、外部 PNJL 参数映射、`a=0.0108805` 的 IMC、外部 `eB` 单位换算、
`muB/3` 语义，以及固定的 MFIR/Hurwitz-zeta + Landau thermal route。它的目标是复核
公式、单位和 solver branch 是否可以对齐。

`production-parity` 使用 Julia 当前默认 `hbarc=197.3269804 MeV fm`、当前运行时模型
配置和相同的 MFIR route。除了 `hbarc`，该 profile 还会随 Julia 默认参数注入、内部
fm 单位换算和生产边界一起变化；因此它不是把 source CSV 的数值简单改写一个常数。
两种 profile 的输出严格写入不同 CSV，禁止混合比较。

## 结果摘要

- fixed-state 诊断：每个 profile 12 行，覆盖 `T={50,240} MeV`、`eB={0.2,0.8} GeV^2`
  的三个锚点和四档节点。
- source-parity 在 matched 节点对外部固定状态保持很小 residual；`T=240 MeV,
  eB=0.8 GeV^2` 的 `n_max=79` 与 `n_max=95` 结果也显示应单独制定 Landau 截断策略。
- production-parity 把外部状态作为初值重新求解后，所有定向点均收敛；其根与外部
  状态的距离约为 `4.5e-5`--`1.6e-4`，说明 production profile 必须独立重求解。
- solver 诊断使用 3 个 seed、trust-region primary、Newton fallback 和候选去重。
  source-parity 在高温高场点保留两个物理候选；production-parity 本轮三个锚点均只
  得到一个候选。该结果不证明全分支完备性。

## 文件

- `source_parity_fixed_state_v1.csv`：source-parity 外部状态固定点评估；
- `production_parity_fixed_state_v1.csv`：production-parity 对同一外部状态的固定点评估，
  仅作 profile 差异诊断，`omega_comparable=false`；
- `source_parity_solver_v1.csv`：source-parity 独立多 seed 重求解和候选输出；
- `production_parity_solver_v1.csv`：production-parity 独立多 seed 重求解和候选输出；
- `manifest.json`、`provenance.json`：schema、输入、脚本、hash 和结果边界。

生成命令：

```powershell
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode both
```

本目录所有结果均为 `diagnostic_only`。在获得更完整的 n_max / pz / zeta 收敛矩阵、
分支覆盖证据和明确的生产 profile 审核前，不得将任何 CSV 作为正式 acceptance baseline。
