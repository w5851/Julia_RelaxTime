# TransportSamplingProfiles

`Models.TransportSamplingProfiles` 把输运计算的背景、A、热平均及末端积分绑定为同一版本配置。
配置文件位于 `config/models/pnjl/transport_sampling/`。两个扫描 CLI 均默认使用 `thermal_weighted_v1`；
旧固定网格、扣除/加回、w0cdf/PCHIP 和无身份 CSV 缓存加载实现已移除。
完整数值验证状态见[统一状态说明](../../../analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)。

```julia
profiles = Models.TransportSamplingProfiles
profile = profiles.load_profile("thermal_weighted_v1")
runtime = profiles.build_runtime(profile)
result = profiles.solve_profile_transport(runtime, equilibrium, T_inv_fm, mu_inv_fm;
    xi=xi, compute_bulk=true, propagator_xi_policy=:match_thermo)
```

`mu_inv_fm` 是夸克化学势，`xi > -1`。`equilibrium` 必须由同一背景模型、温度、化学势及 ξ 求得；
它需要提供 `x_state`、`masses`、`mu_vec` 等 workflow 平衡态字段。
上层扫描使用实际模型下的多种子稳定解选择，旧相图仅用于起始参考。

## 配置与实际积分

| 层次 | `thermal_weighted_v1` |
| --- | --- |
| PNJL 背景 | 热上限 25 fm⁻¹，p96/cos24，残差上限 10⁻⁸ |
| 极化与耦合 A | 半无穷 p32/cos16，scale=5 fm⁻¹；极化与耦合均使用 full-RS A |
| τ 热平均 | 半无穷 p48/cos36/phi36，scale=10 fm⁻¹；截面的 t 积分 24 点 |
| 热平均内部密度 | 半无穷 p32，scale=10 fm⁻¹，cos36 |
| workflow 密度 | 有限上限 25 fm⁻¹，p64/cos8 |
| η、ζ、σ 末端积分 | 半无穷 p32/cos16，scale=5 fm⁻¹ |
| σ 表预算 | 加权误差 `atol + rtol * abs(rate)`，rtol=10⁻³、atol=10⁻¹² |

物理散射支持仍使用模型的 Λ 截断。半无穷动量积分不代表放宽物理散射支持。
`effective_config(profile)` 返回结构化配置、源路径及 SHA-256；
`profile_model(profile; template=..., thermal_p_max_inv_fm=...)` 只替换 PNJL 模型的热积分上限。
自定义 TOML 可以改变数值参数，但不会自动继承已验收配置的收敛资格。

## 缓存、诊断与失败

按 `RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU` 顺序逐道清空极化缓存、执行 GC 并建表，
随后将完整 17 道 `ThermalWeightedCache` 传给 workflow。极化缓存采用有效位量化，
改变跨道缓存生命周期或调用次序需要重新比较数值。

`on_channel(row)` 在每道成功建表后调用，提供 `channel`、`seconds`、`rate`、`estimated_error`、
节点/求值/事件数和停止原因。返回值包含原 workflow 结果以及
`sampling_diagnostics` 和 `sampling_profile`。

任一道预算耗尽、无效输入或 bulk 导数失败都会抛出异常；不会退回旧缓存或无 ζ 的结果。
缓存误差仅控制给定外层节点下的 σ 插值，不证明外层求积、背景或末端输运已收敛。
最终六种 τ 和 η/s、ζ/s、σ/T 的 1% 资格来自独立及联合数值比较。

两个扫描入口均接受 `--sampling-profile thermal_weighted_v1`，拒绝与它冲突的旧网格、扣除和有限末端上限选项。
`<result_stem>_sampling.csv` 保存逐道诊断；`<result_stem>_sampling_contract.json` 绑定配置和数值源码。
续算时身份不一致或缺失会报错。`effective_config.json` 沿用现有包装格式，实际选项位于 `options` 中。

完整验证入口见[采样验证指南](../../../guides/scripts/transport_sampling_qualification.md)。

## 固定背景与物理处方对照

`ward_moller_aligned_v1` 和 `ward_moller_sphere_v1` 保留上述求积设置，
使用 `sampling.relative_velocity="medium_moller"`。可选 `[physics]` 节要求
`a_method` 与 `b0_direction` 两项；未提供该节时沿用 `quadrature` 与 `aligned`。
新 profile 使用 `a_method="ward_matched"`，以变分 RS 能隙方程的热凝聚定义
同步极化与耦合 A。`b0_direction` 取 `aligned` 或 `sphere_average`。
后者在构造线性两线圈时把占据中的极角平方替换为 `1/3`，保持在壳色散、
虚部支集及 `k=0` 处方。它是平均圈函数的近似，不是对完整散射振幅平方取方向平均。

```julia
profile = profiles.load_profile("ward_moller_sphere_v1")
runtime = profiles.build_runtime(profile)
state = profiles.fixed_point_state(profile, quark_params, thermo_params)
result = profiles.solve_frozen_transport(runtime, quark_params, thermo_params,
    species_densities, entropy_inv_fm3; bulk_coeffs=isentropic_coefficients,
    on_channel=(row, cache)->nothing)
```

该入口接收原始背景、六种物种密度、熵与同背景的等熵导数，逐道重建散射率后
计算六种 τ 及 η、ζ、σ；不求解新能隙或相图。`existing_rates` 仅供上层已验证
源码／profile／输入身份的断点续算；输入单位沿用 fm 自然单位。
`bulk_coeffs` 缺省时 ζ 返回 `NaN`；要求全部三个输运量的生产入口必须提供它。
逐道回调可保存插值节点、通道率和时间。新处方的数值资格须由对应对照建立。

冻结 CSV 重跑入口为 `scripts/relaxtime/run_frozen_transport_scan.jl`：

```powershell
julia --project=. scripts/relaxtime/run_frozen_transport_scan.jl --input mode_a_grid.csv --outdir output --sampling-profile ward_moller_sphere_v1 --shard-index 1 --shard-count 32
```

点表包含原稿的质量、Polyakov 环、T、μ、ξ、密度、熵和原 τ/输运量。
程序先以原 τ 重放三个末端积分，再使用新通道率；体黏滞导数锁定在传入分支。
同一输出目录的断点续算要求输入、profile、源码树和点集合身份完全一致。
