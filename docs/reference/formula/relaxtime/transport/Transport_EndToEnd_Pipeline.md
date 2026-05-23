# 输运主链路：从平衡态到输运系数

本文固定主项目 `relaxtime` 主链路中，从平衡态求解到 `eta / zeta / sigma` 的端到端语义。

适用范围：

- `Models.solve_gap_and_transport(...)`
- `src/models/workflow_apps/TransportWorkflow.jl`
- 以及它向下游委托的 `A / B0 / Pi / propagator / scattering / tau / transport` 模块

不适用范围：

- Beth-Uhlenbeck 相移链
- 介子数密度 / `K/pi` 比工作流
- 把 strange-light 通道显式拆成 charge-resolved `K^+ / K^-` 的外部文献记号

## 1. 主链路总图

```text
equilibrium
-> quark_params / thermo_params
-> A_f and B0_ff'
-> Pi_ff'^(P/S)
-> propagator D
-> |M|^2
-> dσ/dt
-> σ(s)
-> wbar_ij
-> tau_i
-> eta / zeta / sigma
```

其中：

- `f, f'` 是有序 flavor 对
- `P/S` 分别表示赝标量 / 标量通道
- strange-light 非混合通道在主项目里写作 `K` 和 `sigma_K`

## 2. 逐步输入/输出合同

| 步骤 | 入口/对象 | 主要输入 | 主要输出 | 单位与语义 |
|---|---|---|---|---|
| 1 | `solve_gap_and_transport` / 平衡层 | `T_fm`, `mu_fm`, `xi` | `equilibrium` | `T_fm, mu_fm` 为 `fm^-1`；默认 transport workflow 使用 `μ_u=μ_d=μ_s=mu_fm` |
| 2 | `build_equilibrium_params` | `equilibrium`, `T_fm`, `mu_fm` | `quark_params`, `thermo_params` | `quark_params=(m, μ)`，`thermo_params=(T, Φ, Φbar, ξ)` |
| 3 | `ensure_quark_params_has_A` | `quark_params`, `thermo_params` | 带 `A=(u,d,s)` 的 `quark_params` | `A_f` 按项目文档口径记为 `fm^-2` |
| 4 | `B0_ff'` | `k0`, `k`, `m_f`, `μ_f`, `m_f'`, `μ_f'`, `T`, `Φ`, `Φbar` | `B0_ff'` | `k0, k, m, μ, T` 为 `fm^-1`；`B0` 在主项目文档口径下视为无量纲复数 |
| 5 | `Pi_ff'^(P/S)` | `A_f`, `A_f'`, `B0_ff'` | `Pi_ff'^(P/S)` | `Pi` 为复数，文档口径单位 `fm^-2` |
| 6 | `meson_propagator_*` | `Pi` 与 `K` 系数 | `D_pi, D_K, D_sigma_pi, D_sigma_K, ...` | `D` 为复数，文档口径单位 `fm^2` |
| 7 | `scattering_amplitude_squared` | 过程键 `process`, `s`, `t`, `quark_params`, `thermo_params`, `K_coeffs` | `|M|^2` | `s,t,u` 为 `fm^-2`；`|M|^2` 文档口径单位 `fm^-4` |
| 8 | `differential_cross_section`, `total_cross_section` | `|M|^2` 与运动学量 | `dσ/dt`, `σ(s)` | `σ(s)` 文档口径单位 `fm^2`；`dσ/dt` 由 `DifferentialCrossSection` 模块消费 |
| 9 | `average_scattering_rate` | `σ(s)`，分布函数，数密度 | `wbar_ij` | 通过 `tau_i^-1 = sum_j rho_j wbar_ij` 进入下一层 |
| 10 | `relaxation_times` | `densities`, `rates` | `tau`, `tau_inv` | `tau` 单位 `fm`，`tau_inv` 单位 `fm^-1` |
| 11 | `transport_coefficients` | `quark_params`, `thermo_params`, `tau`, `bulk_coeffs` | `(eta, zeta, sigma)` | `eta,zeta` 为 `fm^-3`，`sigma` 为 `fm^-1` |

## 3. 平衡层与 flavor 化学势边界

### 3.1 transport workflow 的默认口径

当前主项目 transport workflow 在入口 `solve_gap_and_transport(T_fm, mu_fm; ...)` 中，默认把单个 `mu_fm` 扩成：

```text
μ = (μ_u, μ_d, μ_s) = (mu_fm, mu_fm, mu_fm)
```

这意味着：

- 默认主链路不会因为 `μ_u != μ_s` 而产生 charge-resolved strange-light splitting
- 在默认 workflow 里，`K` / `sigma_K` 只是 generic strange-light 通道标签

### 3.2 更低层 API 的能力边界

更低层的 `quark_params.μ` 本身允许传入：

```text
μ = (μ_u, μ_d, μ_s)
```

因此在模块能力上：

- `B0_ff'`
- `Pi_ff'`
- propagator
- scattering

都可以消费 `μ_u != μ_s` 的 flavor-asymmetric 输入。

但需要注意：

- 当前主链路并没有把 `(u,s)` 与 `(s,u)` 显式暴露成两个独立的 meson symbol
- 所以默认 workflow 不是“只算某个带电 kaon”，而是“没有做 charge-resolved meson labeling”

## 4. `quark_params` / `thermo_params` 的对象边界

主链路中最常见的两个对象是：

```text
quark_params = (
    m = (u=..., d=..., s=...),
    μ = (u=..., d=..., s=...),
    A = (u=..., d=..., s=...)   # 进入 Pi / propagator / scattering 前补齐
)

thermo_params = (
    T = ...,
    Φ = ...,
    Φbar = ...,
    ξ = ...
)
```

其中：

- `m, μ, T` 都是 `fm^-1`
- `A` 在当前项目文档口径中是 `fm^-2`
- `Φ, Φbar, ξ` 无量纲

## 5. `A_f` 与 `B0_ff'` 的 flavor-order 语义

### 5.1 `A_f`

`A_f` 只带单 flavor 标记：

```text
A_u, A_d, A_s
```

它依赖：

```text
(m_f, μ_f, T, Φ, Φbar)
```

当前项目数值口径下：

- 真空/发散部分使用 cutoff 处理
- 热分布部分按半无穷积分目标做大截断近似

### 5.2 `B0_ff'`

`B0_ff'` 是有序的两线泡图对象，`(f,f')` 不能随意交换。

它的核心组合量是：

```text
lambda_ff' = k0 + μ_f - μ_f'
```

并且 `(f,f')` 的顺序还同时决定：

- 哪一条内部线使用 `(m_f, μ_f)`
- 哪一条内部线使用 `(m_f', μ_f')`
- 四项分解里每一项的分布函数归属

因此当 `μ_f != μ_f'` 时，通常不能把：

```text
B0_ff'
```

直接当成：

```text
B0_f'f
```

在某些 single-strange 传输实现路径里，当前项目会再套一层 `k0 <-> -k0` 对称化修正；但那是建立在有序 `B0_ff'` 之上的额外数值处置，不改变 `B0_ff'` 本身是 ordered object 这一点。

## 6. `Pi_ff'^(P/S)` 的语义

主项目极化函数写作：

```text
Pi_ff'^(P/S)
```

其结构可概括为：

```text
Pi_ff'^(P/S)
= A_f + A_f' + [mass_combo_(P/S) - lambda_ff'^2 + k^2] * B0_ff'
```

其中：

- 赝标量 `P` 通道使用 `mass_combo_P = (m_f - m_f')^2`
- 标量 `S` 通道使用 `mass_combo_S = (m_f + m_f')^2`

这一步是 strange-light `K` 与 `sigma_K` 开始分叉的地方：

- `K` 走 `P` 通道
- `sigma_K` 走 `S` 通道

## 7. propagator 层的 strange-light 通道

主项目里 strange-light 非混合传播子是：

```text
D_K       = 2 K4_plus  / (1 - 4 K4_plus  * Pi_us^P)
D_sigma_K = 2 K4_minus / (1 - 4 K4_minus * Pi_us^S)
```

这里的默认主线约定是：

- `K`：strange-light 赝标量通道
- `sigma_K`：strange-light 标量通道
- 主线默认 flavor order 取 `(f,f') = (u,s)`

如果 `μ_u = μ_s`，那么 `(u,s)` 与 `(s,u)` 会简并，generic `K` 标签足够使用。

如果 `μ_u != μ_s`，则从理论对象上讲 `(u,s)` 与 `(s,u)` 应被视为不同 ordered channel；只是当前主链路没有把它们拆成两个独立 symbol。

## 8. scattering 到 tau 的对象边界

### 8.1 `|M|^2`

`|M|^2` 不是直接从某个单一 propagator 得到，而是：

- 先按过程键 `process` 查散射道
- 再按 `t/u/s` 道组合对应的 propagator
- 最后做色自旋平均

因此：

```text
process + kinematics + quark_params + thermo_params + K_coeffs
-> |M|^2
```

### 8.2 `σ(s)`

总截面层把：

```text
|M|^2
```

变成：

```text
σ_ij->cd(s)
```

这是后续 `wbar_ij` 的直接输入。

### 8.3 `wbar_ij`

平均散射率对象：

```text
wbar_ij
```

是对分布函数加权后的过程平均。它不直接区分单个传播子，而是已经把：

- 截面
- 运动学
- 分布函数
- 数密度归一

全部折叠进来。

### 8.4 `tau_i`

最后按：

```text
tau_i^-1 = sum_j rho_j * wbar_ij
```

把过程平均折叠到 flavor / antiflavor 级别：

```text
tau = (u, d, s, ubar, dbar, sbar)
```

## 9. `K / sigma_K / K4_plus / K^+` 的符号边界

这几个符号不能混用。

| 符号 | 在主项目主链路中的含义 | 不是 |
|---|---|---|
| `K` | strange-light 赝标量 meson channel | 不是自动等同于 charge-resolved `K^+` 或 `K^-` |
| `sigma_K` | strange-light 标量 meson channel | 不是“带负电 kaon” |
| `K4_plus` | strange-light 赝标量有效耦合分支 | 不是 charged `K^+` |
| `K4_minus` | strange-light 标量有效耦合分支 | 不是 charged `K^-` |
| `K^+ / K^-` | 只有在显式讨论 charged meson labeling 的文献或 workflow 中才使用 | 不能回代成主链路里的 `K4_plus / K4_minus` |

最容易出错的地方是：

```text
K4_plus  != K^+
K4_minus != K^-
sigma_K  != K^-
```

### 9.1 历史实现与当前代码语义

当前仓库里，generic strange-light 非混合通道与 charged label 不是同一批引入的：

- 最早进入介子质量 / Mott 主线的是 generic `K` 与 `sigma_K`
- 后来 charged / BU / meson-density workflow 才补入 `K_plus` 与 `K_minus`

因此在当前代码语义下，需要再区分两层：

1. generic strange-light label
   - `K`
   - `sigma_K`
2. charged / ordered label
   - `K_plus`
   - `K_minus`

其中当前实现最容易让人误读的是：

- 在已支持 charged label 的介子质量 / 介子数密度路径里，`K` 与 `K_plus` 都落在 ordered `(u,s)` 赝标量分支
- `K_minus` 落在 reversed `(s,u)` 赝标量分支

所以若你后来采用某篇文献的 charged-kaon 记号，并把 ordered `(u,s)` 解释成 `K^+`，那么：

```text
旧的 generic K 在当前实现上更接近 “K_plus-like branch”
```

但这一点**不能**直接平移到 `sigma_K`，因为：

- `sigma_K` 当前仍是 generic strange-light 标量标签
- 主项目尚未把它拆成 charge-resolved `sigma_K_plus / sigma_K_minus`

因此更稳妥的说法是：

```text
旧的 generic K   在当前实现上与 charged-plus 语义近似同向
旧的 generic sigma_K  还不能直接叫作 charged-plus scalar kaon
```

## 10. 建议的阅读顺序

若你要沿主链路排查一个 strange-light 问题，建议按这个顺序读：

1. [docs/api/relaxtime/workflow/TransportWorkflow.md](../../../../api/relaxtime/workflow/TransportWorkflow.md)
2. [OneLoopIntegral_B0.md](../integrals/OneLoopIntegral_B0.md)
3. [Polarization_极化函数byB0.md](../polarization/Polarization_极化函数byB0.md)
4. [Propagator_传播子byPolarization.md](../propagator/Propagator_传播子byPolarization.md)
5. [docs/api/relaxtime/scattering/ScatteringAmplitude.md](../../../../api/relaxtime/scattering/ScatteringAmplitude.md)
6. [docs/api/relaxtime/scattering/DifferentialCrossSection.md](../../../../api/relaxtime/scattering/DifferentialCrossSection.md)
7. [docs/api/relaxtime/scattering/TotalCrossSection.md](../../../../api/relaxtime/scattering/TotalCrossSection.md)
8. [docs/reference/formula/relaxtime/transport/RelaxationTime_FromAverageScatteringRate.md](./RelaxationTime_FromAverageScatteringRate.md)
9. [docs/reference/formula/relaxtime/transport/TransportCoefficients_FromRelaxationTime.md](./TransportCoefficients_FromRelaxationTime.md)
