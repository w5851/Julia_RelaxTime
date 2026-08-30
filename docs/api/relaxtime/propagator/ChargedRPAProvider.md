# ChargedRPAProvider

`ChargedRPAProvider` 是 Phase C 的第一层 bubble-provider 契约。它把
`ChargedRPAKernelSpec` 的有序味道对映射到当前
`PolarizationAniso.polarization_aniso`/`polarization_with_width` 调用，并返回
极化值及完整输入元数据。

当前实现是可替换 provider 的诊断适配器：它复用已有的 `A`/`B0` 正则化，但不
声称完成严格 retarded 解析延拓、极点求解、相移或 BU 积分。

## 入口

```julia
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec
using Main.RelaxTime.ChargedRPAProvider: charged_polarization

spec = charged_rpa_spec(:K_plus; retarded_convention=:external_retarded)
sample = charged_polarization(
    spec,
    k0_inv_fm,
    q_inv_fm,
    (u=m_u, d=m_d, s=m_s),
    (u=mu_u, d=mu_d, s=mu_s),
    (T=T_inv_fm, Φ=Phi, Φbar=PhiBar, ξ=xi),
    (u=A_u, d=A_d, s=A_s),
)
Pi_us = sample.value
```

## 有序输入

provider 不会把共轭通道的输入折叠成同一个数组：

| `spec.meson` | `spec.pair` | `m1, μ1, A1` | `m2, μ2, A2` | 默认 `num_s_quark` |
|---|---|---|---|---:|
| `:K_plus` | `(:u,:s)` | `u` | `s` | 1 |
| `:K_minus` | `(:s,:u)` | `s` | `u` | 1 |
| `:pi_plus` | `(:u,:d)` | `u` | `d` | 0 |
| `:pi_minus` | `(:d,:u)` | `d` | `u` | 0 |

返回值同时保留 `pair`、`kernel_pair`、`m1_inv_fm`、`m2_inv_fm`、
`μ1_inv_fm`、`μ2_inv_fm`、`A1_inv_fm2`、`A2_inv_fm2` 以及 `num_s_quark`，
便于审计 `Pi_us` 与 `Pi_su` 的来源。

## 模式与单位

- `mode=:real_axis`（默认）直接调用 `polarization_aniso`，并要求
  `gamma_inv_fm=0`。
- `mode=:finite_width` 调用现有 `polarization_with_width`。它是在实轴 `B0`
  上组合有限宽度项的兼容 proxy，不等于已证明的 retarded continuation。
- `masses`、`chemical_potentials`、`k0_inv_fm`、`q_inv_fm`、`thermo.T` 的单位为
  `fm^-1`；`A_values` 为 `fm^-2`；返回的 `value` (`Pi_a`) 为 `fm^-2`。

所有输入必须有限，`q`、`gamma` 非负且 `T>0`。provider 不自动计算或重积分
`A_f`；调用方需显式传入同一背景下的 `A_u/A_d/A_s`。

## 后续门禁

本模块已锁定当前 A/B0 计算链的有序参数传递，并用直接调用 parity 测试验证数值
一致性。以下事项仍未完成：

1. 严格 retarded `Pi_{us}`/`Pi_{su}` 的解析延拓和半平面约定；
2. 有限宽度、阈值和电荷共轭关系的物理验证；
3. 传播子极点、Mott 分支、相移边界和 Levinson 计数；
4. strict-support/凝聚处理以及 BU 密度接入。

实现：`src/relaxtime/ChargedRPAProvider.jl`。测试：
`tests/unit/relaxtime/test_charged_rpa_provider.jl`。
