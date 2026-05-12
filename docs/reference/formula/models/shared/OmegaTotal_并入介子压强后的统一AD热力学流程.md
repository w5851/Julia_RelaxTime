# Ω_total：并入介子压强后的统一 AD 热力学流程

本文档专门固定一个系统约束：

> 当介子压强被纳入 EOS 后，项目中的派生热力学量仍然必须走 `Models` 主域统一的自动微分流程；  
> 正确做法是替换总巨热力学势 `\Omega`，而不是在某个 workflow 里局部重写 `s / \epsilon / I`。

这条规则适用于当前 `meson thermo` 主线，也适用于之后任何新的相关压强修正项。

---

## 1. 起点：当前主域统一热力学骨架

当前仓库的主域热力学骨架在：

- `src/models/thermo_kernel.jl`

其统一口径是：

```math
P = -\Omega,
```

```math
\rho_i = \frac{\partial P}{\partial \mu_i}
= -\frac{\partial \Omega}{\partial \mu_i},
```

```math
s = \frac{\partial P}{\partial T}
= -\frac{\partial \Omega}{\partial T},
```

```math
\epsilon = -P + \sum_i \mu_i \rho_i + Ts.
```

在代码层对应为：

1. `model_pressure = -omega(...)`
2. `model_rho` 用 `ForwardDiff.gradient`
3. `model_thermo` 用 `ForwardDiff.derivative`

因此该层真正依赖的唯一物理输入是：

```math
\Omega(T,\mu,\text{state}).
```

---

## 2. 并入介子压强前后的巨热力学势

### 2.1 平均场基线

若没有介子相关压强，主线使用：

```math
\Omega_{\mathrm{MF}}.
```

对应：

```math
P_{\mathrm{MF}} = -\Omega_{\mathrm{MF}}.
```

### 2.2 并入介子压强后

一旦介子压强写成

```math
P_{\mathrm{meson}} = \sum_M P_M,
```

总压强就应写成

```math
P_{\mathrm{total}}
=
P_{\mathrm{MF}} + P_{\mathrm{meson}}.
```

这等价于把巨热力学势替换成：

```math
\Omega_{\mathrm{total}}
=
\Omega_{\mathrm{MF}} - P_{\mathrm{meson}}.
```

若进一步区分 QP / LD，则：

```math
P_{\mathrm{meson}}
=
\sum_M \left(P_M^{\mathrm{QP}} + P_M^{\mathrm{LD}}\right),
```

```math
\Omega_{\mathrm{total}}
=
\Omega_{\mathrm{MF}}
- \sum_M P_M^{\mathrm{QP}}
- \sum_M P_M^{\mathrm{LD}}.
```

---

## 3. 系统标准：派生量继续从 Ω_total 自动导出

一旦定义了 `\Omega_total`，其余量不应再独立手写。

### 3.1 压强

```math
P_{\mathrm{total}} = -\Omega_{\mathrm{total}}.
```

### 3.2 数密度

```math
\rho_i
=
-\frac{\partial \Omega_{\mathrm{total}}}{\partial \mu_i}.
```

### 3.3 熵密度

```math
s
=
-\frac{\partial \Omega_{\mathrm{total}}}{\partial T}.
```

### 3.4 能量密度

```math
\epsilon
=
-P_{\mathrm{total}}
+ T s
+ \sum_i \mu_i \rho_i.
```

### 3.5 Trace anomaly

```math
I(T,\mu)
=
\frac{\epsilon - 3P_{\mathrm{total}}}{T^4}.
```

这就是 `meson thermo` 主线应该服从的正式系统公式。

---

## 4. 为什么不应把 workflow 层的数值差分写成标准口径

在快速原型阶段，workflow 层可能临时用：

- `P_meson(T+\Delta T) - P_meson(T-\Delta T)`

来近似 `s_meson`。  
但这只能是：

1. 原型验证手段；
2. 对缺失 `\Omega_total` 实现时的临时桥接；
3. 数值 sanity check。

它不能成为正式系统公式，原因有三：

1. **破坏统一热力学治理**  
   项目已经有 `omega -> model_pressure/model_rho/model_thermo` 的统一主线，再为 meson thermo 单独维护一套 `entropy/energy/trace anomaly` 差分逻辑，会让不同主线的数值语义分叉。

2. **破坏可组合性**  
   后续如果再引入别的修正项，例如：
   - `sigma`
   - `eta`
   - self-consistent backreaction
   - magnetic / rotation / vector 修正  
   最稳的方式始终是更新 `\Omega_total`，而不是继续叠加更多局部差分补丁。

3. **不利于自动微分与回归治理**  
   当前项目的导数、热力学响应、susceptibility 已经围绕自动微分组织；  
   介子压强进入主线后，应该复用这套基础设施，而不是旁路出去。

---

## 5. 对 meson thermo 主线的实现约束

因此对当前 `meson thermo` 主线，应明确以下工程约束：

1. `relaxtime` 侧负责：
   - `P_M`
   - 或等价的 `\Omega_M = -P_M`
   - 以及必要时的 `QP / LD` 分拆

2. `models` 侧负责：
   - 组装 `\Omega_total`
   - 通过 `omega(model, ...)` 暴露总巨热力学势
   - 复用 `model_pressure / model_rho / model_thermo`

3. workflow / scan / CSV 层负责：
   - 调用统一主线
   - 输出字段
   - provenance / plot-review / regression 资产

换句话说，正式版本不应是：

```text
MesonThermoWorkflow 自己差分出 entropy / epsilon / trace_anomaly
```

而应是：

```text
Meson pressure -> Omega_total -> Models.model_thermo -> entropy/epsilon/trace_anomaly
```

---

## 6. 与 Maslov & Blaschke 2023 的关系

这篇文献给的是：

- 介子 off-shell correlation 如何并入 EOS；
- `QP / LD` 如何在 pressure 中进入；
- mesonic contribution 如何影响 trace anomaly。

而当前项目标准流程补充的是：

- 一旦 `P_meson` 确定，`s / \epsilon / I` 不应再在 workflow 层临时定义；
- 它们应统一从 `\Omega_total` 通过 AD 导出。

因此，把文献方向与本项目标准流程结合后的最终规范是：

```math
\Omega_{\mathrm{total}}
=
\Omega_{\mathrm{MF}}
- \sum_M P_M^{\mathrm{BU/off\text{-}shell}},
```

然后：

```math
\{P,\rho,s,\epsilon,I,\chi,\ldots\}
\Leftarrow
\Omega_{\mathrm{total}}
```

全部通过统一主域流程得到。

---

## 7. 当前阶段的文档结论

本文档为后续实现给出一个明确裁决：

1. `MesonThermodynamics` 文档负责定义 `P_M`；
2. `Omega_total` 文档负责规定系统接线方式；
3. 正式 EOS 派生量以 `\Omega_total` 为唯一源头；
4. 任何 workflow 层的局部差分热力学都只能视作过渡实现，不是最终标准口径。
