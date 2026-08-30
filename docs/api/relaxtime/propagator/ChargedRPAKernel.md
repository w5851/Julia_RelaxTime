# ChargedRPAKernel

`ChargedRPAKernel` 是严格 charged-RPA/BU 复核的 Phase-B 契约层。它把每个
`π^±`/`K^±` 的味道顺序、KMT 通道、传播子归一化候选和 retarded 约定作为显式
元数据，并在给定 `K_a` 与 `Pi_a` 后执行纯标量代数。

这是 diagnostic/experimental backend，不会自动替换
`MesonPropagator.meson_propagator_simple` 或 `MesonDensity`，也不计算夸克泡、极点、
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

| `meson` | `pair`（传给极化 provider 的有序味道） | `kernel_pair` | 物理组成 |
|---|---|---|---|
| `:pi_plus` | `(:u, :d)` | `:K12` | `u\bar d` |
| `:pi_minus` | `(:d, :u)` | `:K12` | `d\bar u` |
| `:K_plus` | `(:u, :s)` | `:K45` | `u\bar s` |
| `:K_minus` | `(:s, :u)` | `:K45` | `s\bar u` |

因此在非同位旋对称背景下，两个 charged kaon 都使用 `K45`；旧
`K4567` 只在兼容语义下等价于 `K67`，不能作为 `K^±` 的默认完整耦合。
`K03/K38` 仍属于中性 `(0,3,8)` 矩阵，本模块不会把它们注入 charged 标量分母。

## 归一化候选

`charged_rpa_spec` 的 `normalization` 只选择候选公式并留下来源标签：

| `normalization` | `normalization_source` | 分子/分母 |
|---|---|---|
| `:legacy_scalar`（默认） | `:legacy_scalar_diagnostic` | `D_a=2K_a/(1-4K_a Pi_a)` |
| `:matrix_scalar_candidate` | `:neutral_matrix_candidate` | `D_a=2K_a/(1-2K_a Pi_a)` |

第二行不是已确认的生产定义；中性矩阵式的因子 `2` 不能在没有重新定义 charged
泡函数和介子场归一化时机械搬用。后续必须用同位旋对称极限、真空介子质量/
Goldstone 条件和极点/相移一致性决定唯一形式。

## API 与单位

- `charged_rpa_spec(meson; channel=:P, normalization=:legacy_scalar,
  retarded_convention=:external_retarded)`：返回不可变的
  `ChargedRPAKernelSpec`。`retarded_convention` 当前仅作 provenance，不执行解析延拓。
- `charged_rpa_coupling(kernel, spec)`：从 `FullKMTInteraction` 读取 `K12/K45/K67`。
- `charged_rpa_inverse(spec, K_a, Pi_a)`：返回无量纲的
  `1 - c_den K_a Pi_a`，不添加极点 epsilon。
- `charged_rpa_propagator(spec, K_a, Pi_a)`：返回
  `c_num K_a/(1-c_den K_a Pi_a)`；`D_a` 的单位是 `fm^2`。

`K_plus` 与 `K_minus` 的 `pair` 字段分别固定为 `(:u,:s)` 和 `(:s,:u)`。
调用方应把对应的 `Pi_us`、`Pi_su` 分别传入同一函数；契约不会假设两者相等。
在 `phi_u=phi_d` 时，测试还核对了默认 `:legacy_scalar` 公式与旧
`MesonPropagator` 的代数值和近极点位置相容。

项目自然单位约定为：`K_a` 为 `fm^2`，`Pi_a` 为 `fm^-2`，所以 `K_a Pi_a` 无量纲，
传播子 `D_a` 为 `fm^2`。输入必须有限；精确零分母会抛出 `DomainError`，不会用
`MesonPropagator` 的历史 epsilon 静默平滑。

## 边界与验证

本模块只验证通道映射、候选标量代数和输入单位/有限性边界。它不声称已经解决：

1. `Pi_{us}`/`Pi_{su}` 的 retarded 解析延拓和有限宽度；
2. 极点残差、阈值和 Mott 分支；
3. 相移常数、Levinson 边界及 generalized-BU 权重；
4. 正常相 Bose 支撑、凝聚零模和节点收敛；
5. 介子热力学反馈或生产路线晋升。

实现：`src/relaxtime/ChargedRPAKernel.jl`。测试：
`tests/unit/relaxtime/test_charged_rpa_kernel.jl`。
