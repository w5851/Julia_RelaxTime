# MesonRPA

`MesonRPA` 是完整 KMT 后端的中性 RPA 代数层。它与
`MesonInteractionKernel` 并行工作，输入固定平均场背景生成的 `FullKMTInteraction`，以及调用方已经计算好的三味 flavor-diagonal 夸克泡
`(Pi_u, Pi_d, Pi_s)`。

## 入口

```julia
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction
using Main.RelaxTime.MesonRPA

kernel = build_full_kmt_interaction((phi_u, phi_d, phi_s); G=G, K=K)
Pi = neutral_polarization_matrix((Pi_u, Pi_d, Pi_s))
D = neutral_rpa_propagator(kernel, Pi; channel=:P)
```

如果需要从当前 `PolarizationAniso` 实际计算三味同味泡，可使用诊断层
[`MesonRPAAdapter.md`](MesonRPAAdapter.md)；它不会改变本模块的纯代数入口。

`neutral_polarization_matrix` 按 Tian 等 Phys. Rev. D 114, 034012 (2026), Eq. (26) 组装 `(lambda_0, lambda_3, lambda_8)` 基底的 3×3 极化矩阵。`neutral_rpa_propagator` 严格使用

```text
D = 2 K * inv(I - 2 K * Pi)
```

的矩阵顺序；`neutral_rpa_inverse_matrix` 返回分母矩阵，`neutral_rpa_determinant` 返回其行列式。

## 边界

该模块只做有限维复矩阵代数，不计算极化泡、不寻找极点、不生成相移，也不把 `Omega_M` 写回 PNJL 平衡方程。它不修改旧的 `MesonPropagator` 2×2 接口，当前结果只能作为 diagnostic。数值接入前仍需核对 `PolarizationAniso` 的归一化、味道顺序和动量/化学势约定。

公式、来源和适用边界见 `docs/reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md`；测试见 `tests/unit/relaxtime/test_meson_rpa.jl`。
