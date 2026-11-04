# 平均散射率与热权重缓存

`average_scattering_rate(process, quark_params, thermo_params, K_coeffs; ...)`
默认使用热权重自适应采样和精确角度合并。外层默认 p48/cos36/phi36，截面的 t 积分默认 24 点；
完整生产配置还包括背景、A 和末端积分，见 [TransportSamplingProfiles](TransportSamplingProfiles.md)。

## 积分与角度合并

平均散射率是分布加权的相对速度与截面平均。将整体方位角积分提出后，系数为
`d_q²/(32π⁵ n_i n_j)`，其余积分为
`∫dp_i dp_j p_i²p_j² ∫dc_i dc_j dφ f_i f_j v_rel σ(s)`。
数密度默认在半无穷动量区间积分；各向异性分布对极角余弦为偶函数。

外层动量使用 `p=scale*t/(1-t)`，默认 scale=10；
截面的物理支持仍由模型 Λ 决定，二者不是同一个截断。
严格等质量入射态使用 `v_rel=2sqrt((s-4m²)/s)`，其他情况使用一般表达式。
单位：s 为 fm⁻²，截面和平均散射率为 fm²。

角度合并只合并**完全相同的相对夹角**，把各自的分布和求积权重相加。
它保留各向异性占据数，不要求等质量、相同化学势或对称的自定义网格。
热权重采样在此离散测度上分配截面插值误差；两项优化同时使用。

## 缓存接口

```julia
using .AverageScatteringRate
config = ThermalSamplingConfig(rtol=1e-3)
options = (p_nodes=48, angle_nodes=36, phi_nodes=36,
           n_sigma_points=24, sampling_config=config)
cache = build_thermal_weighted_cache(:uubar_to_uubar,
    quark_params, thermo_params, K_coeffs; options...)
rate = average_scattering_rate(:uubar_to_uubar,
    quark_params, thermo_params, K_coeffs; options..., cs_cache=cache)
sigma = get_sigma(cache, s)
direct_sigma = get_sigma(cache, s; direct=true)
```

未传 `cs_cache` 时自动构建；已构建的 `ThermalWeightedCache` 保存小表、平均速率、误差估计和统计，
释放构建时的大型离散热测度。复用时强制核对物理参数、传播子输入、A、耦合、
全部有效节点及权重、物理截断、t 节点和误差预算。上下文改变必须重建，表数组应当只读。

`get_sigma` 在物理区间外返回 0；在物理区间内但超出已认证表支持时抛出 `DomainError`。
`direct=true` 可局部直算，不改变表、速率或认证状态。它不是另一套完整热积分。

## 误差与失败语义

`ThermalSamplingConfig` 默认 `rtol=1e-3, atol=1e-12`。
阈值窗口 `weighted_window=0.6` 的单位为 fm⁻²，在该窗口插值 `(s-s_th)σ`。
取点结合阈值、传播子实轴特征、热积分贡献及两轮错位探测，不使用 ξ 例外清单。

误差估计只控制**给定外层求积下的截面插值误差**，不认证外层、t 积分、B₀、背景或末端积分。
评估数、节点、轮数或事件预算耗尽时抛出异常，不能得到有效缓存。
入射密度必须有限且正；自定义网格必须成对提供，权重必须为正。

`sigma_cache_policy=:thermal_weighted` 是唯一策略；
`interpolation_mode=:auto|:thermal_weighted` 和 `angular_reduction=:grouped` 仅保留输入校验。
旧 `CrossSectionCache`、w0cdf/PCHIP 建表、阈值扣除/加回、无身份 CSV 表加载和旧网格参数已移除。
历史实验使用各自保存的源码快照重现，不能以新默认运行后覆盖旧结果。

## 传播子与密度

默认 `propagator_xi_policy=:match_thermo` 使传播子和外层分布使用同一 ξ。
诊断选项 `:isotropic` 只把传播子上下文置为 ξ=0，外层分布和密度仍用真实 ξ；
可显式提供 `propagator_quark_params`。未给 A 时沿用独立的 A 构造默认值（角度 4 点），
完整生产配置则显式提供已验证的 A，避免把外层角度节点误当作 A 节点。

`number_density` 默认 p32、极角 36、scale=10。也接受显式的成对网格；
有限动量区间必须独立检查尾项。τ 聚合所用外部密度由调用方负责，见 [RelaxationTime](RelaxationTime.md)。

分层收敛、同点复现、源码快照和性能证据见
[输运采样验收](../../../guides/scripts/transport_sampling_qualification.md)。
