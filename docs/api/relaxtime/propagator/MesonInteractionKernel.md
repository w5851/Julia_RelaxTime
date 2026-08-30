# MesonInteractionKernel

`MesonInteractionKernel` 是 Phase 1.5 的并行、纯代数 KMT 相互作用核后端。它从给定平均场凝聚
`(phi_u, phi_d, phi_s)` 和原始耦合 `(G, K)` 生成：

- 中性 `(lambda_0, lambda_3, lambda_8)` 基底的 3×3 P/S 耦合矩阵；
- `pi` (`K_12`)、`u-s` kaon (`K_45`) 和 `d-s` kaon (`K_67`) 的味道分辨带电/中性通道。

它不求解 PNJL 平衡，不计算极化函数/相移，也不把介子压力反馈进 `Omega`。新类型与旧
`EffectiveCouplings.calculate_effective_couplings` 和 `MesonPropagator` 2×2 API 并行存在。

`MesonRPA` 是其并行的 Phase 2 代数消费者：`neutral_polarization_matrix` 将
`(Pi_u, Pi_d, Pi_s)` 组装到 `(0, 3, 8)` 基底，`neutral_rpa_propagator` 按
`2K * inv(I - 2K * Pi)` 计算完整中性矩阵。它仍不计算夸克泡、不找极点、不做相移，也不修改
PNJL 平衡或旧传播子入口。

## 入口

```julia
using Main.RelaxTime.MesonInteractionKernel

kernel = build_full_kmt_interaction(
    (-0.30, -0.28, -0.12);
    G=G_fm2,
    K=K_fm5,
)

Kp = neutral_coupling_matrix(kernel, :P)
Ks = neutral_coupling_matrix(kernel, :S)
Kpi = charged_coupling(kernel, :K12, :P)
Kus = charged_coupling(kernel, :K45, :P)
Kds = charged_coupling(kernel, :K67, :P)
K03 = get_Kab(kernel, 0, 3, :P)
```

也可从现有 equilibrium facade 返回值读取 `x_state[1:3]`：

```julia
kernel = build_full_kmt_interaction_from_equilibrium(equilibrium; G=G_fm2, K=K_fm5)
```

该适配器只读取前三个标量凝聚，不触发新的 gap solve。

## 类型与字段

- `KernelConvention`：记录 `lambda_0` 归一化、凝聚定义和候选公式来源。
- `ChargedKMTChannels`：字段 `K12`、`K45`、`K67`，单位 `fm^2`。
- `FullKMTInteraction`：保存 `G`、`K`、`phi`、P/S 中性矩阵、P/S 带电通道、基底和约定元数据。
- `neutral_coupling_matrix(kernel, channel)`：`channel` 可用 `:P`/`:S`，也接受 `:pseudoscalar`/`:scalar`。
- `charged_couplings(kernel, channel)`：返回对应 `ChargedKMTChannels`。
- `charged_coupling(kernel, pair, channel)`：读取 `:K12`、`:K45` 或 `:K67`。
- `get_Kab(kernel, a, b, channel)`：读取 `a,b ∈ (0,3,8)` 的中性矩阵元。

数值均采用项目内部自然单位；`phi=(phi_u,phi_d,phi_s)` 是物理夸克凝聚
`<bar q_f q_f>`，单位为 `fm^-3`（不是无量纲量）；`G` 为 `fm^2`，
`K` 为 `fm^5`，输出耦合为 `fm^2`。因此 `K*phi_f` 的单位也是 `fm^2`。
调用方必须显式传入与 PNJL 配置一致的符号和单位。

## 非对称背景下的通道映射

完整核直接使用平衡解中的 `phi_f`，不需要先把它转换成旧 helper `H_f`。
若为了兼容旧 `EffectiveCouplings` 接口而定义 `H_f=-phi_f`，则有：

| 旧接口/代数通道 | 完整 KMT 通道 | 物理味道组成 |
|---|---|---|
| `K123`（`H_s`） | `K12` | `pi^±` |
| `K4567`（`H_u`） | `K67` | `K^0=d\bar s`、`\bar K^0=s\bar d` |
| 无对应旧独立字段 | `K45` | `K^+=u\bar s`、`K^-=s\bar u` |

在 `phi_u=phi_d` 的同位旋对称极限，`K45=K67`，所以旧的 `K4567`
数值与 charged/neutral 两个完整通道耦合；在 `phi_u!=phi_d` 时，
`K4567` 只能按 `K67`（u spectator）解释，不能继续作为 `K^±` 的默认耦合。

## 与 RPA 分层的关系

这里的 `K_ab` 由固定平均场背景上的 KMT Hartree 收缩得到；RPA 随后在同一
背景上构造二次介子涨落并重求和夸克泡图。极化函数 `Pi` 不仅依赖质量和外部
动量，还依赖 `T`、各味 `mu_f`、Polyakov 背景、正则化及所选外部运动学。
因此本模块只负责核的纯代数构造，不能被解释为已经完成完整 RPA 或热力学反馈。

## 稳定性边界

目前这两个模块属于 diagnostic/Phase 1.5-2 后端。它们不会自动替换旧生产接口，也不会因为 `phi_u != phi_d` 就修改 PNJLCore 的 stationarity。当前 RPA 接口只完成矩阵代数；夸克泡数值接入、极点/相移和带电 K/π 生产 A/B 对照属于后续阶段。

对应公式、符号风险和文献来源见 `docs/reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md`；测试见 `tests/unit/relaxtime/test_meson_interaction_kernel.jl` 和 `tests/unit/relaxtime/test_meson_rpa.jl`。
