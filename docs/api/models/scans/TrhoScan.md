# `Models.run_trho_scan`

本页从 `Models` 统一入口说明 T-ρ 扫描。对于相结构、Maxwell、CEP 与 phase pipeline 前置数据生成，它通常比 T-μ 扫描更关键。

## 定位

`Models.run_trho_scan` 是 T-ρ 网格扫描的首选业务入口。它适合：

- 生成 `(T, ρ, ξ)` 网格结果
- 为 S 形检测、Maxwell、spinodal、crossover、CEP 提供前置数据
- 在统一入口层控制 `reverse_rho`、`seed_policy` 与约束模式

实现转发位于 `src/models/entrypoints.jl`，底层执行位于 `src/models/scans/TrhoScan.jl`。

## 入口形态

从当前实现看，核心关键字参数包括：

- `T_values`
- `rho_values`
- `xi_values`
- `output_path`
- `overwrite`
- `resume`
- `reverse_rho`
- `seed_policy`
- `continuation`（`:auto`、`:legacy`、`:single_pass`、`:palc`）与 `palc_options`
- `constraint_mode`
- `hybrid_weighted_fallback`
- `solver_backend`
- `p_num`, `t_num`
- `thermo_quadrature_policy`（`:tensor_gauss` 或 PNJL 标量热核的 `:rs_reduced_adaptive`）
- `thermo_quadrature_rtol`, `thermo_quadrature_atol`, `thermo_quadrature_maxevals`

返回值与 T-μ 扫描一致，也是 `(total, success, failure, skipped, output)` 风格的命名元组。

## 为什么它是相图前置层入口

和 `Models.run_tmu_scan` 相比，`Models.run_trho_scan` 更直接地决定相图链路的数据质量，因为：

- `ρ` 网格分辨率会影响 S 形检测与 Maxwell 闭合
- `reverse_rho` 会影响连续性跟踪是否稳定
- `DEFAULT_RHO_VALUES` 本身就是研究级默认采样模板的一部分

这也是为什么 scan 主题要把 `build_default_rho_grid` 作为公开辅助入口显式写出来。

## 默认网格与扫描顺序

`:auto` 对 PNJL、tensor 积分、默认 ground_state 选择器且每条等温线至少 21 个目标点的
`FixedRho` / `FixedAsymmetricRho` 扫描使用 `:palc`；较短路径、RPNJL 及其它配置沿用
`:legacy`。该阈值来自已完成的 21–25 点与中等网格对照；显式 `:legacy` 可复用旧策略。
`:palc` 对每个 `(T, ξ)` 进行多起点双向延拓，保留回弯
上的候选，随后在请求的 ρ/ρ₀ 点作精确联立校正并按既有硬约束和最大压力规则选优。
`:single_pass` 每个不同种子求解一次；FixedAsymmetricRho 复用既有候选展开实现。
新策略支持 PNJL/RPNJL、tensor 积分与默认 ground_state selector。其它组合可以用
`:legacy`；显式不兼容组合报 `ArgumentError`。
自定义 `solver` 或路径未支持的求解器选项由 `:auto` 保留在旧链路；显式新策略在扫描前报错。

`palc_options` 传入 `PALCContinuation`；`path_diagnostics_cb((T, xi), path_result)`
可读取起点、分支与工作量。CSV/callback 仍按 `reverse_rho` 控制的输入顺序输出。
resume 保留已有行；部分完成的等温线会重建 PALC 路径，仅写出缺失点。

若不显式提供 `rho_values`，入口默认使用 `DEFAULT_RHO_VALUES`。若目标是相图主线，一般不建议跳过这一默认多分辨率策略，除非你已经明确知道需要怎样的局部加密。

同时，默认 `reverse_rho=true`，因为从高密度向低密度推进通常更稳，能降低低密度端 continuity 失败风险。

这条规则不是经验小技巧，而是主题级使用合同：

- `ρ=0` 附近更容易触发连续性失稳
- 在 S 形或 Maxwell 敏感区，过粗的 `ρ` 网格会直接破坏判据质量
- 若 Maxwell 失败集中在少数温度切片，优先局部加密 `ρ`，不要先盲目增大全局积分参数

## 输出合同

T-ρ 扫描的稳定跨模块产物同样是 CSV。当前 header 为：

```text
T_MeV,rho,xi,mu_u_MeV,mu_d_MeV,mu_s_MeV,mu_avg_MeV,mu_B_MeV,mu_Q_MeV,mu_S_MeV,pressure_fm4,entropy_fm3,energy_fm4,rho_u_fm3,rho_d_fm3,rho_s_fm3,phi_u,phi_d,phi_s,Phi1,Phi2,M_u_MeV,M_d_MeV,M_s_MeV,iterations,residual_norm,converged,message
```

字段语义要点：

- `mu_*_MeV` 来自 flavor 化学势分解，`mu_avg_MeV` 为三味平均
- `rho_u_fm3`、`rho_d_fm3`、`rho_s_fm3` 是分味粒子数密度
- 失败或近似接受的点仍保留在输出中，并由 `message` 区分原因

同一个点的 joint solve、压力/密度/熵后处理、residual 复算和分味密度写出必须共享同一组热积分参数；入口不会把自适应策略只用于残差而在 CSV 后处理退回固定网格。
当前 `:rs_reduced_adaptive` 在本入口明确限定为 `model_kind=:PNJL` 且 `constraint_mode=:fixed_rho`；其他模型或约束模式会提前报错，不会静默回退到 tensor。

## 推荐与禁用组合

- 常规研究或相图主线：`reverse_rho=true`，`rho_values=Models.build_default_rho_grid()`
- 仅 smoke 或联调：可降低 `p_num`、`t_num` 并改用更粗 `rho_values`
- 不建议：`ρ=0` 与 `reverse_rho=false` 同时出现，或用 `Δρ > 0.1` 的粗网格直接做 Maxwell / CEP 结论

## 相关主题

- [Overview.md](Overview.md)
- [Algorithms.md](Algorithms.md)
- [SamplingGrid.md](SamplingGrid.md)
- [../phase/README.md](../phase/README.md)
