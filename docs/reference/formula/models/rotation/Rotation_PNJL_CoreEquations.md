# Rotation-PNJL 核心公式（T-mu-omega 最小可运行口径）

## 1. 适用范围与定位

本文档提炼 legacy rotation 文档中的关键内容，服务于“先适配壳、后增强”的迁移策略：

- 最小变量集：`(phi, Phi, PhiBar)`
- 最小控制参量：`(T, mu, omega)`
- 最小输出：`pressure, rho, entropy, energy`

## 2. 有效能量与旋转修正

在 rotation 近似中，单粒子能量可写为：

```math
E_{rot}(p,n;M,\omega)=\sqrt{p^2+M^2}-(n+1/2)\,\omega,
```

其中 `n` 是与旋转相关的离散模标签（具体截断与权重由数值离散决定）。

## 3. PNJL 有效势与对数项

巨热力学势沿用 PNJL 结构，但把热项中的能量替换为 `E_rot`：

```math
\Omega_{rot}
=\Omega_{cond}(\phi)
+U(\Phi,\bar\Phi;T)
\;+
\Omega_{vac}(M)
\;+
\Omega_{th}(E_{rot},\mu,T,\Phi,\bar\Phi).
```

对应热项中两类对数核：

```math
\mathcal Q_1^{rot}=\ln\!\left(1+3\Phi e^{-(E_{rot}-\mu)/T}+3\bar\Phi e^{-2(E_{rot}-\mu)/T}+e^{-3(E_{rot}-\mu)/T}\right),
```

```math
\mathcal Q_2^{rot}=\ln\!\left(1+3\bar\Phi e^{-(E_{rot}+\mu)/T}+3\Phi e^{-2(E_{rot}+\mu)/T}+e^{-3(E_{rot}+\mu)/T}\right).
```

## 4. 平衡态方程（最小）

平衡态由极值条件确定：

```math
\frac{\partial\Omega_{rot}}{\partial\phi}=0,\qquad
\frac{\partial\Omega_{rot}}{\partial\Phi}=0,\qquad
\frac{\partial\Omega_{rot}}{\partial\bar\Phi}=0.
```

若采用固定密度路径 `T-rho`，需附加：

```math
R_\rho = \rho_{calc}(T,\mu,\omega)-\rho_{target}=0.
```

因此常见为 4 变量联立系统（3 个序参量 + 1 个化学势或约束变量）。

## 5. 热力学量定义（与主线一致）

```math
P=-\Omega_{rot},\qquad
\rho=-\frac{\partial\Omega_{rot}}{\partial\mu},\qquad
s=-\frac{\partial\Omega_{rot}}{\partial T},\qquad
\epsilon=Ts+\mu\rho-P.
```

## 6. 与现有主线的边界约定

- 本文档只定义最小 rotation 核，不绑定 legacy 的具体节点实现细节。
- 迁移实现时建议先复用现有 `Models` 求解策略接口，而不是复制旧脚本控制流。
- 高阶导数（如 `dP/dT^n`, `dP/domega`）列为 Phase-B 增强，不在首批最小核强制实现。

## 7. 代码映射建议

- 公式层：本页
- 代码层（计划）：`src/models/variants/rotation/`
- 测试层（计划）：
  - `tests/unit/models/test_rotation_model.jl`
  - `tests/integration/models/test_rotation_workflow_smoke.jl`
