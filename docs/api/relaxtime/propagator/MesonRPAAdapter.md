# MesonRPAAdapter

`MesonRPAAdapter` 是 Phase 3 的诊断数值桥接层。它把当前
`PolarizationAniso` 的三次同味夸克泡计算接到 `MesonRPA` 的中性
`(0, 3, 8)` 纯矩阵代数，保持两个后端的职责边界清楚。

## 入口

```julia
using Main.RelaxTime.MesonRPAAdapter

q = (
    m=(u=m_u, d=m_d, s=m_s),
    μ=(u=mu_u, d=mu_d, s=mu_s),
    A=(u=A_u, d=A_d, s=A_s),
)
t = (T=T, Φ=Phi, Φbar=PhiBar, ξ=xi)

bubbles = neutral_flavor_bubbles(:P, k0_inv_fm, k_norm_inv_fm, q, t; ensure_A=false)
Pi = Main.RelaxTime.MesonRPA.neutral_polarization_matrix(bubbles)

result = neutral_rpa_from_quark_params(
    kernel, k0_inv_fm, k_norm_inv_fm, q, t;
    channel=:P,
    ensure_A=false,
)
```

`neutral_flavor_bubbles` 返回顶层字段 `u`、`d`、`s`，因此可以直接传给
`neutral_polarization_matrix`。同时返回 `channel`、`k0_inv_fm`、
`k_norm_inv_fm`、`gamma_inv_fm`、`with_width`、`num_s_quark`、
`a_auto_built` 以及归一化后的 `quark_params`/`thermo_params`，方便保存诊断
元数据。

`neutral_rpa_from_quark_params` 另外返回 `polarization`、`inverse_matrix`、
`propagator` 和 `determinant`，并保留输入 `kernel`。它只是组合已有函数，
不改变 `FullKMTInteraction` 或 RPA 矩阵公式。

## 数值合同

- 新入口的位置、质量、化学势和温度参数均使用项目内部自然单位 `fm^-1`；
  这层不把 `MeV` 自动转换为 `fm^-1`。
- 对每个味 `f`，调用 `PolarizationAniso` 的同味对 `(f,f)`，即
  `m1=m2=m_f`、`mu1=mu2=mu_f`、`A1=A2=A_f`。`(0,3,8)` 的非对角极化
  元素来自三味泡的基底变换，而不是这里新增的非对角夸克自能。
- 默认 `num_s_quark=(u=0,d=0,s=0)`，表示明确的同味泡约定。若需复现
  旧调用方的标签，可以逐味传入其他非负整数；目前 `PolarizationAniso`
  只有 `1` 会触发 `k0` 对称平均，`2` 仅作为兼容标签保留。
- `with_width=false` 时 `gamma_inv_fm` 必须为零；使用宽度时设置
  `with_width=true`，并传入非负 `gamma_inv_fm`。
- 缺少 `A` 时默认通过 `AFieldBuilder.ensure_quark_params_has_A` 补齐；可用
  `ensure_A=false` 强制要求调用方预计算 `A`，并通过 `a_p_nodes`、
  `a_p_max`、`a_cos_nodes` 和 `a_use_aniso` 控制自动补值。

## 诊断边界

该 adapter 不运行 PNJL gap solver，不修改平均场驻点方程，不引入介子
`Omega_M` 反馈，不寻找极点/相移，也不替换旧 `MesonPropagator` 2×2 或
`MesonDensity` 生产入口。它只验证“当前泡实现能否按指定味道顺序进入完整
中性 RPA 代数”的数值接线；归一化、正则化、极点和生产级相移仍需单独的
同位旋对称/不对称以及外部文献固定点门禁。

对应实现：`src/relaxtime/MesonRPAAdapter.jl`。
对应测试：`tests/unit/relaxtime/test_meson_rpa_adapter.jl`。
