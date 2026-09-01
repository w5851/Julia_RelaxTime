# PhaseNormalization

`PhaseNormalization` 是 charged-RPA/BU 路线的纯代数约定层。它只连接散射相位、
on-shell scalar/diagonal `S`-matrix 和 BU 能量测度，不从传播子自动推断物理相位，
也不改变 `MesonDensity` 或任何 production 默认。实现位于
`src/relaxtime/PhaseNormalization.jl`。

## 标量闭合

本模块固定两个可互换的相位变量：

| 变量 | 定义 | 测度因子 |
|---|---|---:|
| `:delta` | 物理散射相位 `delta`，`S=exp(2im*delta)` | `d_delta/pi` |
| `:s_matrix_argument` | `arg(S)=2*delta` | `d_arg_s/(2pi)` |

```julia
S = phase_to_s_matrix(delta)
delta_back = s_matrix_to_phase(S)
rho = s_matrix_density_of_states(S, dS)
```

对 `S=exp(2im*delta)`，模块验证

```math
\operatorname{Im}(S^{-1}\partial_\omega S)=2\partial_\omega\delta,
\qquad
\frac{1}{2\pi}\operatorname{Im}(S^{-1}\partial_\omega S)
=\frac{1}{\pi}\partial_\omega\delta.
```

因此 `domega/pi` 与 `domega/(2pi)` 不是两种可以混用的“整体归一化”，而是分别
对应 `delta` 和完整 `arg(S)` 两种相位变量。共同因子在同一口径的 ratio 中可能抵消，
但会改变单电荷绝对密度，不能据此跳过绝对密度回归。

## 分支与传播子映射边界

`s_matrix_to_phase` 返回主值相位；连续分支和有限网格高能端点由
`BUPhaseGates.anchor_phase_high_energy` 负责。该锚定只去除相位常数并记录 tail，
不证明 `omega_max` 已足够大。

现有 `ChargedPhaseBackend` 使用
`delta=-arg(Delta^R_inverse)` 作为项目诊断映射。该映射必须由同一传播子、解析延拓
和 `S`-matrix 归一化的独立验证支持；本模块的代数恒等式不把它升级为文献唯一算法。

对多通道对角 `S`，`s_matrix_log_derivative(S,dS)` 返回
`Im tr(S^-1*dS)`，即 `d arg(det(S))`。若单独跟踪 eigenphase，仍须记录通道基底和
分支选择；非对角 coupled-channel 的生产实现不由本模块自动提供。

## 来源与测试

`S=exp(2i*delta)` 及 `d delta/pi` 来自 Dashen--Ma--Bernstein 的 S-matrix 统计力学
公式以及其现代两体 BU 重写；项目文献证据与转换边界见
`docs/analysis/relaxtime/charged_phase_literature_review_v1/README.md`。

测试：`tests/unit/relaxtime/test_phase_normalization.jl`，覆盖标量/对角矩阵、合成
连续 profile、主值分支恢复和因子 2 测度恒等式。
