# ChargedRPAProvider

`ChargedRPAProvider` 是 Phase-C 的 ordered real-axis bubble 适配层。它把
`ChargedRPAKernelSpec` 的有序味道映射到现有
`PolarizationAniso.polarization_aniso`，并把 strict-route candidate 与 legacy
`num_s_quark=1` oracle 明确分开。

该模块仍是 diagnostic adapter：它复用当前 `A/B0` 计算，不声称已经完成
second-sheet 极点、有限宽度物理验证、相移或 BU 积分。

## 入口

```julia
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec
using Main.RelaxTime.ChargedRPAProvider: charged_polarization

spec = charged_rpa_spec(:K_plus)
ordered = charged_polarization(
    spec,
    k0_inv_fm,
    q_inv_fm,
    (u=m_u, d=m_d, s=m_s),
    (u=mu_u, d=mu_d, s=mu_s),
    (T=T_inv_fm, Φ=Phi, Φbar=PhiBar, ξ=xi),
    (u=A_u, d=A_d, s=A_s),
)
Pi_us = ordered.value
```

## 处方边界

`prescription` 只有两个值：

| 处方 | strange channel 的 `num_s_quark` | 定位 |
|---|---:|---|
| `:ordered_retarded`（默认） | `0` | strict GBU 的 ordered real-axis candidate |
| `:legacy_symmetrized_B0` | `1` | Rehberg/旧 Fortran/Cpp 对齐 oracle |

对 pion 通道两种处方均使用 `num_s_quark=0`。legacy 处方只在每个有序 strange
通道内部做 `k0/-k0` 平均；它不会把 `(u,s)` 与 `(s,u)` 合并。因此 `K+` 和
`K-` 始终分别返回 `Pi_us` 与 `Pi_su`。

显式 legacy 调用：

```julia
legacy = charged_polarization(
    spec,
    k0_inv_fm,
    q_inv_fm,
    masses,
    chemical_potentials,
    thermo,
    A_values;
    prescription=:legacy_symmetrized_B0,
)
```

## 返回字段与单位

返回值保留 `pair`、`kernel_pair`、`prescription`、`num_s_quark`、有序
`m1/m2`、`mu1/mu2`、`A1/A2` 和热力学背景。`analytic_scope=:real_axis` 明确表示
当前没有把 `B0` 延拓到复平面。

- 质量、化学势、`k0`、`q` 和 `T`：`fm^-1`；
- `A_values` 和返回的 `value=Pi_a`：`fm^-2`；
- `Phi`、`PhiBar`、`xi`：无量纲。

所有输入必须有限，`q>=0` 且 `T>0`。provider 不自动计算 `A_f`，调用方必须传入
同一平均场背景和同一正则化口径下的 `A_u/A_d/A_s`。

## 未完成门禁

1. ordered 与 legacy 在 `mu=0`、有限化学势固定点和冻结线上的对照；
2. legacy `Gamma` 的来源绑定固定点；
3. second-sheet pole、阈值和 Mott 分支；
4. 相移、Levinson、Bose support 和 BU 密度接入。

权威公式见
[`ChargedRPA_BU_ProductionRoute.md`](../../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)。
实现：`src/relaxtime/ChargedRPAProvider.jl`。测试：
`tests/unit/relaxtime/test_charged_rpa_provider.jl`。
