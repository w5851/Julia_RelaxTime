# ChargedPhaseBackend

`ChargedPhaseBackend` 是与旧 `MesonDensity` 并行的 strict charged-RPA/BU 诊断后端。
它不求解 PNJL 平衡、不做 second-sheet 极点延拓，也不改变任何 production 默认。
实现位于 `src/relaxtime/ChargedPhaseBackend.jl`。

## 公式约定

输入是有序 retarded 逆传播子 `Δᵣ(ω,q)`。默认相位为

```math
δ(ω,q)=-\operatorname{arg}Δ^R(ω,q),
```

其中主值相位从高能端向低能端连续 unwrap，并以有限窗口高能端为锚点。后端同时
保留 raw、unwrapped、anchored 三条 profile，并报告高能 tail span；tail 未稳定时
结果会保留数值值但标记 `accepted=false`。

单电荷严格 BU 使用

```math
n_M=\frac{d_M}{T}\int\frac{dq\,q^2}{2\pi^2}
\int\frac{d\omega}{\pi},g_B(\omega)\frac{\partial δ_M}{\partial\omega}.
```

该后端的 `domega/pi` 是默认值；`domega/(2pi)` 只能作为 legacy ratio adapter。
相位分支、阈值和 Levinson 条件必须先通过，不能用常数 anchor 替代收敛检查。

真实 bubble 可以来自两条独立诊断路径：`ChargedRPAProvider` 的
`:ordered_retarded`（有限上半平面 `eta`）或 `:ordered_pv_cut`（实轴 Cauchy
主值实部加解析 cut）。后者不等同于把 `eta` 设成一个很小的正数；两者的
`eta`、节点、端点和高能尾部都必须分别做收敛对照。

## 主要 API

### `StrictChargedPhaseSpec(; ...)`

不可变的数值合同，字段包括：

- `phase_object=:inverse_propagator` 或 `:propagator`；
- `phase_sign=-1`（逆传播子默认）或 `+1`；
- 高能目标、branch tolerance、tail 点数和 tail tolerance。

### `strict_phase_profile(omega, inverse_values; spec=...)`

返回相位和端点诊断。它只做相位构造，不声称 profile 已满足 Levinson。

### `strict_phase_gate(profile; threshold, bound_state_count, ...)`

组合阈下 inverse 实部符号变号根计数、Levinson 阈值相位和高能 tail gate。失败
时返回 `passed=false`，不静默修正相位。

### `strict_charged_bu_density(inverse_fn, mass, T; ...)`

对 `inverse_fn(omega, q)` 执行 `q`/`omega` 积分，返回 `density`、`q_profiles`、
每个 `q` 的 phase/gate 状态、Bose 支撑、测度和节点/截断配置。设置
`require_levinson=true` 时必须提供 `threshold` 和 `bound_state_count`；否则只执行
端点诊断，适合合成路径和数值定位。两者也可以是 `q -> value` callable，以表达
动量依赖的连续阈值和显式绑定态计数。负密度或非有限密度会保留原始值但返回
`status=:invalid_density`、`accepted=false`。

### `strict_charged_rpa_bu_density(spec, coupling, polarization_fn, mass, T; ...)`

这是 strict phase backend 与真实 charged-RPA provider 的组合适配层。它把

```math
\Delta^{-1,R}_{a}(\omega,q)=1-4K_a\Pi^{R}_{a}(\omega,q)
```

交给 `strict_charged_bu_density`，其中 `polarization_fn` 必须返回有序 retarded
复数极化函数。分析脚本
`scripts/analysis/relaxtime/audit_charged_phase_backend.jl` 使用
`ChargedRPAProvider` 的 `:ordered_retarded` 路径在一个固定 BQS 背景上调用该适配层，
并把每个 charged 通道的门禁失败写入本地诊断 CSV；脚本不会切换任何生产默认。

### `strict_mott_gate(before_profile, after_profile; ...)`

对 Mott 转变前后的两个 profile 分别执行端点/Levinson gate，再调用
`BUPhaseGates.mott_phase_gate` 检查绑定态数和阈值相位的同步下降。绑定态数是显式
输入，不由有限网格自动猜测。

### `strict_density_convergence_gate(coarse, refined; rtol=..., atol=...)`

同时检查两次结果是否有限、两次 strict gates 是否通过以及密度是否满足显式相对/
绝对误差门禁。节点和截断变化必须通过该函数记录，不能只比较 ratio。

## 与旧路径的边界

现有 `phase_shift_meson_number_density` 保持原有 `current`/`gbu_reference`、
`phase_anchor` 和 legacy measure 语义。新后端当前是 solver-independent diagnostic
接口；将真实 `ChargedRPAProvider` profile 接入并通过 `eta`、`omega_max`、节点、
Levinson/Mott 和 Bose-support 门禁，属于后续生产候选评审，不在本 API 中静默切换。

测试：`tests/unit/relaxtime/test_charged_phase_backend.jl`。
