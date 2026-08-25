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
- 高温高场锚点枚举了 3 个 seed 的全部 6 种顺序。source-parity 每种顺序均得到相同
  的两个候选，production-parity 每种顺序均得到相同的一个候选；两个输出在第二次
  独立执行后 SHA-256 逐字节不变。该检查只证明已声明 seed 集合内的顺序不敏感和
  确定性，不证明 seed 集合或物理解分支已经完备。

## 截断与节点诊断

在 `p_num=128, zeta_num=256, pz_max=40 fm^-1` 下，`n_max` 从 `31` 扫到 `511`：

- `T=50 MeV, eB=0.2 GeV^2` 在全部档位上稳定；
- `T=240 MeV, eB=0.8 GeV^2` 在 `n_max=63` 后已达到约 `1e-8 fm^-4` 的稳定性；
- `T=240 MeV, eB=0.2 GeV^2` 到 `n_max=383` 后才稳定，`383 -> 511` 的 Omega
  变化约 `4.8e-11 fm^-4`。外部 source state 在 `n_max=79` 恰好满足外部
  source-parity，但在收敛截断上的 residual 约为 `1.68e-4`，所以外部
  `n=0..79` 不能作为截断收敛证明。

真正重求解也确认这一点。production-parity 的低场高温根从 `n_max=79` 到
`n_max=383` 发生约 `1.53e-5` 的五维状态位移和 `-1.773e-4 fm^-4` 的 Omega 变化；
低温点不变，高场点的 Omega 只变化约 `-6.36e-10 fm^-4`。

解耦节点矩阵显示，在这三个锚点上 `p_num=128` 和 `pz_max=40 fm^-1` 已与更高探针
对齐到约 `1e-9 fm^-4`；高场点需要至少单独审查 `zeta_num`，`64/128` 明显不足，
`384` 与 `512` 的 Omega 差约 `1.6e-11 fm^-4`。

## 生产 blocker

当前 production 默认是 `p_num=96, zeta_num=64, pz_max=25 fm^-1, n_max=nothing`。
auto `n_max` 在三个锚点分别解析为 `4,3,3`。在两个 `T=240 MeV` 点，外部状态上的
fixed residual 分别约为 `1.284` 和 `0.150`；默认 solver 得到的状态与外部状态距离
分别约为 `0.253` 和 `0.0133`。这不是可由 profile 小差异解释的误差。

实现还会为每个 seed 独立解析 auto `n_max`，再把这些可能来自不同截断问题的候选
按 Omega 排序。因此，当前 auto-`n_max` 生产路线不能作为正式 equilibrium validation
target，也不能用本目录的 `n_max=383` 直接替换成全域常数。后续必须建立共同的、
温度感知的 thermal Landau 尾项截断合同，并让同一点的全部 seed/attempt 使用同一
截断后，才能重新评估 production acceptance。

作为候选合同的短诊断，使用同一点共享的
`E_tail=max(|mu|+k*T,mass)` 推导 `n_max`，并以 `n_max=511` 做局部参考：低场高温
`k=24` 给出 `n_max=249`、Omega 差约 `1.13e-8 fm^-4`，`k=30` 给出 `n_max=389`、
差约 `3.84e-11 fm^-4`；高场点对应差异约 `9.7e-9` 和 `3.4e-11 fm^-4`。低温点
`k=16` 后已稳定。这个结果支持温度感知、同点共享截断的设计方向，但不是全域
接受标准；`k`、尾项误差、真实 root solver 和更多 `T/mu/eB` 点仍需单独审计。

## 文件

- `source_parity_fixed_state_v1.csv`：source-parity 外部状态固定点评估；
- `production_parity_fixed_state_v1.csv`：production-parity 对同一外部状态的固定点评估，
  仅作 profile 差异诊断，`omega_comparable=false`；
- `source_parity_solver_v1.csv`：source-parity 独立多 seed 重求解和候选输出；
- `production_parity_solver_v1.csv`：production-parity 独立多 seed 重求解和候选输出；
- `source_parity_branch_repeatability_v1.csv`：source-parity 高温高场 6 种 seed 顺序；
- `production_parity_branch_repeatability_v1.csv`：production-parity 高温高场 6 种 seed 顺序；
- `source_parity_nmax_fixed_state_v1.csv`、`production_parity_nmax_fixed_state_v1.csv`：
  `n_max=31..511` 的固定态矩阵；
- `source_parity_quadrature_fixed_state_v1.csv`、`production_parity_quadrature_fixed_state_v1.csv`：
  `p_num`、`pz_max`、`zeta_num` 解耦矩阵；
- `source_parity_solver_cutoff_v1.csv`、`production_parity_solver_cutoff_v1.csv`：
  `n_max=79/383` 的独立重求解；
- `production_default_nmax_fixed_state_v1.csv`、`production_default_nmax_solver_v1.csv`：
  当前默认 auto-`n_max` 的实际解析层数和求解结果；
- `source_parity_thermal_cutoff_fixed_state_v1.csv`、`production_parity_thermal_cutoff_fixed_state_v1.csv`：
  `k={12,16,20,24,30,36}` 的温度感知截断候选诊断；
- `manifest.json`、`provenance.json`：schema、输入、脚本、hash 和结果边界。

生成命令：

```powershell
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode both
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode branch
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode cutoff
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode quadrature
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode solver_cutoff
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode default
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl --mode thermal_cutoff
```

本目录所有结果均为 `diagnostic_only`。在修复 auto-`n_max` production blocker、补充
更广的 `T/mu/eB` 收敛矩阵并复核分支覆盖前，不得将任何 equilibrium CSV 作为正式
acceptance baseline。
