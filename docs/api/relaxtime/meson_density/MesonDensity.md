# MesonDensity 模块 API 文档

## 模块概述

`MesonDensity` 模块提供当前介子数密度主线的最小实现入口，当前仅覆盖稳定粒子极限：

- `π/K` 聚合通道默认简并因子
- 玻色分布函数
- 稳定粒子极限数密度
- `K/π` 比值与温度扫描

后续 BU / BW / 各向异性扩展将在该模块基础上继续演化。

## 单位约定

与项目其他 `relaxtime` 模块一致，内部计算使用自然单位制：

- `mass`, `T`, `μ`, `q`：`fm^-1`
- 返回的数密度：`fm^-3`

## API 参考

### `meson_degeneracy(meson; charge_resolved=false)`

返回当前主线的 `π/K` 简并因子。

- 聚合通道：`d_π = 3`、`d_K = 4`
- 电荷分辨通道：`d = 1`

### `bose_distribution(E, μ, T)`

计算玻色分布函数：

```math
g(E) = \frac{1}{\exp((E-\mu)/T)-1}.
```

约束：

- `T >= 0`
- `E > μ`

### `stable_meson_number_density(mass, T; μ=0.0, degeneracy=1, qmax=nothing, num_q_nodes=256)`

计算稳定粒子极限介子数密度：

```math
n_M = d_M \int_0^\infty \frac{dq\,q^2}{2\pi^2}
\frac{1}{\exp((E_M-\mu_M)/T)-1},
\qquad E_M=\sqrt{q^2+m_M^2}.
```

### `stable_kpi_ratio(m_pi, m_K, T; μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...)`

返回稳定粒子极限的聚合通道 `K/π` 数密度比值。

### `stable_kpi_scan(temperatures; m_pi, m_K, μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...)`

对一组温度执行稳定粒子极限扫描，返回：

- `temperatures`
- `n_pi`
- `n_K`
- `kpi_ratio`
