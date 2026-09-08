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
n_M=d_M\int\frac{dq\,q^2}{2\pi^2}
\int\frac{d\omega}{\pi}\,g_B(\omega)\frac{\partial δ_M}{\partial\omega}.
```

该后端的 `domega/pi` 是默认值；`domega/(2pi)` 只能作为 legacy ratio adapter。
相位分支、阈值和 Levinson 条件必须先通过，不能用常数 anchor 替代收敛检查。
导数式无额外 `1/T`；该因子仅在分部积分后的 `g_B(1+g_B) delta` 形式中出现。
数密度单位为 `fm^-3`，单 q 壳层 `dn/dq` 为 `fm^-2`。

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
采样根数也必须等于显式 `bound_state_count`。这仍是全 profile 的保守旧诊断门禁，
不能替代分 cut 的物理 sheet 判定。

### `bu_phase_integral(omega, phase, T; μ=0, weight=:current)`

对有序有限采样执行 `g_B dF/pi` 的 Stieltjes 求积，`weight` 为 `:current` 或
`:gbu`。先对 GBU 权重作差，保留离散跳变，禁止输入 folded display 作为物理分支。
omega、T、μ 使用 `fm^-1`，返回无量纲谱权重。调用方必须固定同一能量坐标。
常数 anchor 不改变普通 BU 的 `d delta`，但一般会改变 GBU 的 `dF(delta)`；
仅整数倍 pi 平移对后者也不变。GBU 的有限窗口 anchor 不是无代价的归一化操作。

### `bu_phase_integral_parts(omega,phase,T; μ=0,weight=:current)`

返回 `derivative`、`bulk`、`boundary`、`lower_boundary`、`upper_boundary`、
`reconstructed` 与 `identity_residual`，均为无量纲谱权重。离散乘积恒等式为
`derivative=bulk+boundary`，其中 `bulk=-sum(mean(F)*diff(g))/pi`。
这不是对连续 `F*g*(1+g)/T` 的独立积分收敛证明。有限窗口、非零端点相位时，
不得丢弃 `boundary`；正的 `bulk` 不保证完整积分为正。

`strict_retarded_phase` 现在拒绝精确复零点；该处角度没有定义，必须用单侧极限
和独立根计数。非零小残差不被任意容差裁掉。

### `split_bu_shell(root_result, segments, q, T; μ=0, weight=:current)`

`root_result` 来自 `BUPhaseGates.certify_gap_roots`，且必须通过 gap 内认证。
`segments` 是按能量排序、互不相交的 `(omega,phase)` 列表，不能包含已认证根。
返回 `bound_shell_inv_fm2`、`continuum_shell_inv_fm2`、`total_shell_inv_fm2`；
每个正能根在 current/GBU 中都贡献一个 Bose 权重。负的连续谱/部分和原样保留。
不自动授予完整束缚态计数、Levinson 或 production 资格。

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

分析脚本 `scripts/analysis/relaxtime/audit_charged_mott_profiles.jl` 在两个显式温度
上配对真实 ordered profile，并输出质量--阈值差、独立束缚态状态、Levinson 与 Mott
gate。它不会自动寻找 Mott 温度；默认温度只是诊断候选，任何 `complex_subthreshold`
或 phase gate 失败都必须保留，不能晋升为物理转变结论。

### `strict_density_convergence_gate(coarse, refined; rtol=..., atol=...)`

同时检查两次结果是否有限、两次 strict gates 是否通过以及密度是否满足显式相对/
绝对误差门禁。节点和截断变化必须通过该函数记录，不能只比较 ratio。

## 与旧路径的边界

现有 `phase_shift_meson_number_density` 保持原有 `current`/`gbu_reference`、
`phase_anchor` 和 legacy measure 语义。新后端当前是 solver-independent diagnostic
接口；将真实 `ChargedRPAProvider` profile 接入并通过 `eta`、`omega_max`、节点、
Levinson/Mott 和 Bose-support 门禁，属于后续生产候选评审，不在本 API 中静默切换。

测试：`tests/unit/relaxtime/test_charged_phase_backend.jl`。
