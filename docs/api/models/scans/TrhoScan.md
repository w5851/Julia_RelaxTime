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
- `constraint_mode`
- `hybrid_weighted_fallback`
- `solver_backend`
- `p_num`, `t_num`

返回值与 T-μ 扫描一致，也是 `(total, success, failure, skipped, output)` 风格的命名元组。

## 为什么它是相图前置层入口

和 `Models.run_tmu_scan` 相比，`Models.run_trho_scan` 更直接地决定相图链路的数据质量，因为：

- `ρ` 网格分辨率会影响 S 形检测与 Maxwell 闭合
- `reverse_rho` 会影响连续性跟踪是否稳定
- `DEFAULT_RHO_VALUES` 本身就是研究级默认采样模板的一部分

这也是为什么 scan 主题要把 `build_default_rho_grid` 作为公开辅助入口显式写出来。

## 默认网格与扫描顺序

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

## 推荐与禁用组合

- 常规研究或相图主线：`reverse_rho=true`，`rho_values=Models.build_default_rho_grid()`
- 仅 smoke 或联调：可降低 `p_num`、`t_num` 并改用更粗 `rho_values`
- 不建议：`ρ=0` 与 `reverse_rho=false` 同时出现，或用 `Δρ > 0.1` 的粗网格直接做 Maxwell / CEP 结论

## 相关主题

- [Overview.md](Overview.md)
- [Algorithms.md](Algorithms.md)
- [SamplingGrid.md](SamplingGrid.md)
- [../phase/README.md](../phase/README.md)