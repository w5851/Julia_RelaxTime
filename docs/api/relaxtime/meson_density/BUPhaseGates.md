# BUPhaseGates

`BUPhaseGates` 提供 real-axis Beth-Uhlenbeck 路线的纯门禁，不计算传播子，也不
修改 production 默认。实现位于 `src/relaxtime/BUPhaseGates.jl`。

## 能量测度

`bu_omega_measure` 和 `bu_omega_measure_factor` 固定两种显式口径：

| canonical symbol | factor | 定位 |
|---|---:|---|
| `:single_charge_domega_over_pi` | `1/pi` | `d=1`、单 Bose 因子的 strict charged 绝对密度 |
| `:legacy_domega_over_2pi` | `1/(2pi)` | 已发表 ratio 写法兼容 adapter |

单束缚态的相移变化为 `pi`，所以 strict 口径返回一个稳定玻色子；legacy 口径的
绝对密度少一半。同一口径的 `K/pi` ratio 不受这个共同因子影响。

## 高能相位锚点

```julia
anchored = anchor_phase_high_energy(omega, phase; target=0.0)
```

函数从最大 `omega` 向低能端 unwrap，再把高能端平移到 `target`。返回值保留
raw、unwrapped、anchored 三条序列，并给出 `applied_shift`、`tail_span` 和
`max_adjacent_jump`。锚定只消除常数相位，不证明高能截断已经收敛；调用方仍需
比较更大的 `omega_max`。

## 阈下根与 Levinson 门禁

```julia
roots = count_subthreshold_roots(omega, inverse_values, threshold)
gate = levinson_phase_gate(
    omega,
    phase,
    threshold;
    bound_state_count=roots.count,
)
```

`count_subthreshold_roots` 只计数阈下 inverse propagator 实部的简单符号变号根。
若阈下虚部超过 `imag_tolerance`，返回 `status=:complex_subthreshold` 和
`passed=false`，不会把复结构静默计作束缚态。

`levinson_phase_gate` 检查：

```math
\delta(\omega_{thr})-\delta(\infty)=\pi n_B,
```

并同时要求高能 tail span 小于指定容差。返回结果包含相位残差、尾部跨度和两个
分门禁状态。

## Mott 门禁

```julia
transition = mott_phase_gate(before_gate, after_gate)
```

默认要求 Mott 前后束缚态数减少 1，阈值相位同步减少 `pi`，且两侧各自通过
Levinson 门禁。它是 fail-able acceptance contract，不替代介子质量阈值或
second-sheet pole 求解。

## Bose 支撑与数值收敛

```julia
bose_support_gate(mass, chemical_potential; omega_min, omega_max)
convergence_gate(reference, candidate; rtol=..., atol=...)
```

`bose_support_gate` 只检查正常相积分域：`m_M > μ_M` 且积分下界满足
`omega_min > μ_M`。它不会添加零动量凝聚模；返回的 `:unsafe_bose_domain` 必须在
扫描输出中保留，不能静默改用 `x_min_cut`。

`convergence_gate` 对一个标量观测量执行显式绝对/相对误差判断。节点或截断扫描应
分别记录低/高配置、参考值、候选值和门禁状态，而不能只保留最终比值。

四类可比较的密度算法由 `four_density_algorithm_labels()` 固定为：
`stable_particle_limit`、`reduced_strict_bw`、`q_pole_strict_bw` 和
`phase_shift_bu`。它们是诊断比较标签，不表示当前已有任何 production 授权。

## 当前边界

- 只提供 real-axis profile 的纯校验；
- sign-change root count 不认证偶重根或复杂 sheet 结构；
- 不处理 Bose 凝聚零模；
- 不授权冻结线或 production 切换。

测试：`tests/unit/relaxtime/test_bu_phase_gates.jl`。
