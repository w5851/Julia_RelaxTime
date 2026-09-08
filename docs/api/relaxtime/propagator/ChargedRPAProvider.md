# ChargedRPAProvider

`ChargedRPAProvider` 是 Phase-C 的 ordered charged bubble 适配层。默认 strict
处方在 `z=omega+i*eta` 上调用 `OneLoopIntegrals.B0_retarded`；显式的
`:ordered_pv_cut` 使用实轴主值加解析 cut 的 retarded 边界值；历史实轴
`PolarizationAniso` 则作为两个显式 legacy oracle 保留。

该模块仍是 diagnostic adapter：它不声称已经完成
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
    eta_inv_fm=1e-3,
    energy_nodes=128,
)
Pi_us = ordered.value
```

## 处方边界

`prescription` 有四个值：

| 处方 | 后端 | strange channel 的 `num_s_quark` | 定位 |
|---|---|---:|---|
| `:ordered_retarded`（默认） | `B0_retarded(z=omega+i*eta)` | `0` | strict GBU 上半平面探针 |
| `:ordered_pv_cut` | `B0_pv_cut`（PV 实部 + retarded cut） | `0` | strict 实轴 PV/cut diagnostic |
| `:ordered_legacy_B0` | 历史 `PolarizationAniso` | `0` | PR289 ordered adapter 兼容 oracle |
| `:legacy_symmetrized_B0` | 历史 `PolarizationAniso` | `1` | Rehberg/旧 Fortran/Cpp 对齐 oracle |

对 pion 通道三种处方均使用 `num_s_quark=0`。symmetrized legacy 处方只在每个有序 strange
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
`m1/m2`、`mu1/mu2`、`A1/A2` 和热力学背景。legacy 返回
`analytic_scope=:real_axis_legacy`；strict 返回
`analytic_scope=:upper_half_plane_probe`，并记录有效 `eta_inv_fm` 和
`energy_nodes`。

- 质量、化学势、`k0`、`q` 和 `T`：`fm^-1`；
- `A_values` 和返回的 `value=Pi_a`：`fm^-2`；
- `Phi`、`PhiBar`、`xi`：无量纲。

所有输入必须有限，`q>=0` 且 `T>0`。strict 后端当前只支持 `xi=0`；非零各向
异性只能显式选择 legacy 处方。provider 不自动计算 `A_f`，调用方必须传入
同一平均场背景和同一正则化口径下的 `A_u/A_d/A_s`。

`:ordered_pv_cut` 返回 `analytic_scope=:real_axis_pv_cut`、`eta_inv_fm=0` 和
`energy_nodes=0`。它不表示已经通过 `eta->0+`、端点、Levinson/Mott 或节点/截断
门禁；这些比较必须由分析脚本显式记录。

## 未完成门禁

1. ordered 与 legacy 在冻结线上的扩展对照；
2. legacy `Gamma` 的来源绑定固定点；
3. second-sheet pole、阈值和 Mott 分支；
4. 相移、Levinson、Bose support 和 BU 密度接入。

权威公式见
[`ChargedRPA_BU_ProductionRoute.md`](../../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)。
实现：`src/relaxtime/ChargedRPAProvider.jl`。测试：
`tests/unit/relaxtime/test_charged_rpa_provider.jl`。
