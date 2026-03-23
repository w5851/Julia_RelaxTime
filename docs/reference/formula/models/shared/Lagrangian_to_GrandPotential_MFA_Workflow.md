# 从拉氏量到巨热力学势：平均场推导统一流程（内化版）

## 1. 文档目的

本文档吸收 legacy `docs/formulas/从拉氏量到巨热力学势.md` 的核心价值，但去除对话式冗余，转化为可直接指导本仓库实现与审阅的统一流程。

适用模型：NJL / PNJL / 其可扩展变体（rPNJL、磁场、rotation 等）。

## 2. 统一七步法

### Step 1：写出模型拉氏量

以三味 PNJL 为例：

```math
\mathcal L
=\bar q(i\gamma^\mu D_\mu-\hat m+\gamma^0\hat\mu)q
+\mathcal L_{4q}
+\mathcal L_{6q}
-\mathcal U(\Phi,\bar\Phi;T).
```

其中 `D_mu` 是否含背景场/外场，决定后续分布与色散的具体形式。

### Step 2：构造配分函数

```math
Z=\int\mathcal D\bar q\,\mathcal D q\;\exp\left[\int_0^\beta d\tau\int d^3x\;\mathcal L_E\right].
```

### Step 3：平均场近似（MFA）

把高阶费米子算符线性化到序参量（凝聚）上，例如：

```math
(\bar q q)^2 \to 2\langle\bar q q\rangle(\bar q q)-\langle\bar q q\rangle^2.
```

得到“双线性夸克项 + 纯场势能项”的结构。

### Step 4：执行夸克路径积分

双线性费米积分给出泛函行列式：

```math
\Omega = -\frac{T}{V}\ln Z = \Omega_{field} - \frac{T}{V}\ln\det(S^{-1}).
```

### Step 5：Matsubara 求和与迹运算

核心技术关系：

```math
\mathrm{Tr}\ln D = \ln\det D,
```

```math
T\sum_n\ln\left(\frac{\omega_n^2+E^2}{T^2}\right)=E+2T\ln(1+e^{-E/T}).
```

PNJL 情况下色迹会把普通费米对数改写为含 `Phi/PhiBar` 的多项式对数核。

### Step 6：组装巨热力学势

最终统一结构应写成：

```math
\Omega = \Omega_{cond} + \Omega_{vac} + \Omega_{th} + \Omega_{Polyakov}(+\Omega_{vector}+\Omega_{field\,ext}).
```

其中：

- `cond`：凝聚与耦合常数项
- `vac`：真空积分（通常需截断/重整）
- `th`：有限温密热激发
- 其余为模型特有扩展项（矢量、磁场、旋转等）

### Step 7：极值条件给出平衡态

```math
\frac{\partial\Omega}{\partial X_a}=0,
```

`X_a` 为全部独立序参量（例如 `phi_u,phi_d,phi_s,Phi,PhiBar`）。

## 3. 与代码实现的映射模板

建议把每个模型都映射到同一接口语义：

- `omega_components(...) -> (cond, vac, thermal, polyakov, extra)`
- `grand_potential(...)` 只做合并
- `solve_gap(...)` / `solve_constraint(...)` 只处理极值或约束闭合

这样模型差异主要体现在 `components` 的定义，不破坏主线求解器。

## 4. 常见误区（迁移时必须避免）

- 把“模型推导流程”与“某个旧脚本控制流”绑定在一起。
- 真空项与热项的单位/截断口径混用。
- 先写大量派生量接口，再补基础 `Omega`，导致语义反转。

## 5. 本仓库内化约束

- 内部单位优先自然单位；对外 MeV 字段显式命名（`*_MeV`）。
- 公式文档必须能映射到现有或计划中的 `src/models/*` 路径。
- 以“可测试最小核”作为第一实现目标，而非一次完成全部变体。

## 6. 旋转 PNJL 的专用补充（arXiv:2307.14402）

对 rotation 变体，Step 4-6 在数学结构上有三点特化：

1) **准粒子能量带角动量模位移**

```math
\epsilon_n=\sqrt{M^2+p_t^2+p_z^2}-(n+1/2)\omega,
```

2) **积分核带柱坐标模权重**

```math
\mathcal W_n(p_tr)=J_n^2(p_tr)+J_{n+1}^2(p_tr),
```

3) **巨势仍可整理为 PNJL 双对数核 + Polyakov 势**

```math
\Omega_{rot}=\Omega_{cond}+\Omega_{int}(\epsilon_n,\mathcal W_n)+\mathcal U(\Phi,\bar\Phi;T).
```

其中 `U(Φ,Φ̄;T)` 可先采用不显含 `omega` 的多项式势，旋转影响通过夸克热项反馈到 `Phi/PhiBar`。

> 该补充用于统一“推导流程文档”与 `rotation/Rotation_PNJL_CoreEquations.md` 的符号和实现口径，避免同一仓库内出现两套不兼容写法。
