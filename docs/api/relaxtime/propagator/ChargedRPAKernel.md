# ChargedRPAKernel

`ChargedRPAKernel` 是 charged-RPA/BU 公式路线的 Phase-B 契约层。它显式记录
`pi^±/K^±` 的有序味道、KMT 通道和项目已经闭合的 charged 单通道归一化，
并在调用方提供 `K_a` 与 `Pi_a` 后执行纯代数传播子计算。

这是 diagnostic/experimental backend，不会自动替换
`MesonPropagator.meson_propagator_simple` 或 `MesonDensity`，也不计算夸克泡、
相移、Beth-Uhlenbeck 积分或 `Omega_M` 反馈。

## 入口

```julia
using Main.RelaxTime.ChargedRPAKernel

spec = charged_rpa_spec(:K_plus; channel=:P)
inverse = charged_rpa_inverse(spec, K_45, Pi_us)
D = charged_rpa_propagator(spec, K_45, Pi_us)
```

从完整 KMT 核读取通道耦合：

```julia
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction

kernel = build_full_kmt_interaction((phi_u, phi_d, phi_s); G=G_fm2, K=K_fm5)
K_45 = charged_rpa_coupling(kernel, spec)
```

## 通道映射

| `meson` | 有序味道 | `kernel_pair` | 物理组成 |
|---|---|---|---|
| `:pi_plus` | `(:u, :d)` | `:K12` | `u\bar d` |
| `:pi_minus` | `(:d, :u)` | `:K12` | `d\bar u` |
| `:K_plus` | `(:u, :s)` | `:K45` | `u\bar s` |
| `:K_minus` | `(:s, :u)` | `:K45` | `s\bar u` |

两个 charged kaon 都使用 `K45`。旧 `K4567` 在 `phi_u != phi_d` 时代数上对应
`K67`，只能作为 legacy/neutral oracle；`K03/K38` 仍属于中性 `(0,3,8)` 块。

## 固定归一化

本模块只实现公式路线闭合包已经确定的形式：

```math
D_a=\frac{2K_a}{1-4K_a\Pi_a}.
```

中性矩阵写法 `2K[I-2K Pi_matrix]^{-1}` 与该式并不冲突；charged ladder 顶角
给出 `Pi_matrix=2Pi_a`。因此不能把同一个 `Pi_a` 原样放入 `1-2KPi_a`，本 API
也不再暴露这一已否决的候选。

`charged_rpa_spec` 不提供 normalization 或 retarded convention 开关，并固定记录：

- `normalization_source=:charged_ladder_goldstone`；
- `retarded_convention=:retarded_e_minus_iwt`，仅作为后续 provider 的约定元数据；
- 纯代数始终使用分子因子 `2` 和分母因子 `4`。

## 单位与错误语义

`K_a` 为 `fm^2`，`Pi_a` 为 `fm^-2`，所以逆传播子无量纲，`D_a` 为 `fm^2`。
输入必须有限。精确零分母会抛出 `DomainError`，不会使用旧传播子的 epsilon
静默平滑。

## 边界

本模块只锁定通道映射、归一化和纯代数。以下事项仍属于后续 gate：

1. ordered real-axis `Pi_us/Pi_su` 的固定点和电荷共轭验证；
2. legacy finite-width 与 second-sheet pole 的分离实现；
3. 相移常数、Levinson/Mott 边界和 generalized-BU 权重；
4. Bose support、凝聚零模和节点收敛；
5. 介子热力学反馈或 production 晋升。

权威公式见
[`ChargedRPA_BU_ProductionRoute.md`](../../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)。
实现：`src/relaxtime/ChargedRPAKernel.jl`。测试：
`tests/unit/relaxtime/test_charged_rpa_kernel.jl`。
