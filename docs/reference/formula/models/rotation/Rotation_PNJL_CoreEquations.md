# Rotation-PNJL 核心公式（基于 arXiv:2307.14402 的最小实现口径）

## 1. 适用范围与定位

本文档用于把旋转 PNJL 的关键公式内化到本仓库可实现语义，主要对应文献：

- Fei Sun, Kun Xu, Mei Huang, *Quarkyonic phase induced by Rotation*, arXiv:2307.14402v1

面向当前主线的最小口径：

- 最小序参量：`(phi, Phi, PhiBar)`
- 最小控制参量：`(T, mu, omega)`
- 最小输出：`pressure, rho, entropy, energy`

## 2. 旋转下的有效能量与质量

在文献的旋转展开与柱坐标模展开下，准粒子能量写为（对应文献 Eq. (13) 的记号）：

```math
\epsilon_n(p_t,p_z;M,\omega)=\sqrt{M^2+p_t^2+p_z^2}-(n+1/2)\,\omega,
```

其中动态质量由平均场给出：

```math
M=m-2G\langle\bar qq\rangle.
```

`n` 为 z 角动量模标签（整数）。

## 3. 巨热力学势（旋转 PNJL）

按文献 Eq. (12)-(13) 的结构，可写为：

```math
\Omega_{rot}
=G\langle\bar qq\rangle^2
-\frac{T}{4\pi^2}
\sum_{n=-\infty}^{\infty}
\int_0^{\Lambda} p_t\,dp_t
\int_{-\sqrt{\Lambda^2-p_t^2}}^{\sqrt{\Lambda^2-p_t^2}} dp_z\;
\mathcal W_n(p_t r)
\left[\mathcal Q_+ + \mathcal Q_-\right]
+\mathcal U(\Phi,\bar\Phi;T),
```

其中旋转模权重

```math
\mathcal W_n(p_t r)=J_n^2(p_t r)+J_{n+1}^2(p_t r),
```

对数核可整理为 PNJL 常见双对数形式：

```math
\mathcal Q_+=\ln\!\left(1+3\Phi e^{-(\epsilon_n-\mu)/T}+3\bar\Phi e^{-2(\epsilon_n-\mu)/T}+e^{-3(\epsilon_n-\mu)/T}\right),
```

```math
\mathcal Q_-=\ln\!\left(1+3\bar\Phi e^{-(\epsilon_n+\mu)/T}+3\Phi e^{-2(\epsilon_n+\mu)/T}+e^{-3(\epsilon_n+\mu)/T}\right).
```

说明：文献中还写出了等价的正负指数对数项组合；在实现时可统一成上面的双对数表达以减少重复积分核。

## 4. Polyakov 势的参数化

文献采用多项式势（Eq. (14)）：

```math
\frac{\mathcal U(\Phi,\bar\Phi;T)}{T^4}
=-\frac{b_2(T)}{2}\Phi\bar\Phi
-\frac{b_3}{6}(\Phi^3+\bar\Phi^3)
+\frac{b_4}{4}(\Phi\bar\Phi)^2,
```

```math
b_2(T)=a_0+a_1\frac{T_0}{T}+a_2\left(\frac{T_0}{T}\right)^2+a_3\left(\frac{T_0}{T}\right)^3.
```

在该文口径里，`U` 不显含旋转参数 `omega`；旋转影响主要通过夸克能谱与热项反馈到 `Phi/PhiBar`。

## 5. 平衡态方程（最小）

平衡态由驻点条件给出（文献 Eq. (15)）：

```math
\frac{\partial\Omega_{rot}}{\partial\langle\bar qq\rangle}=0,\qquad
\frac{\partial\Omega_{rot}}{\partial\Phi}=0,\qquad
\frac{\partial\Omega_{rot}}{\partial\bar\Phi}=0.
```

若采用固定密度路径 `T-rho`，需再联立

```math
R_\rho=\rho_{calc}(T,\mu,\omega)-\rho_{target}=0.
```

对应 4 变量闭合系统（`phi, Phi, PhiBar, mu`）。

## 6. 热力学量定义（与主线一致）

```math
P=-\Omega_{rot},\qquad
\rho=-\frac{\partial\Omega_{rot}}{\partial\mu},\qquad
s=-\frac{\partial\Omega_{rot}}{\partial T},\qquad
\epsilon=Ts+\mu\rho-P.
```

## 7. 数值实现约束（文献口径）

文献数值部分给出可直接迁移的最小约束：

- 费米子参数：`m=0.005 GeV`, `Lambda=0.65 GeV`, `G=4.93 GeV^-2`
- Polyakov 势参数：`a0=6.75`, `a1=-1.95`, `a2=2.625`, `a3=-7.44`, `b3=0.75`, `b4=7.5`, `T0=0.27 GeV`
- `n` 求和快速收敛，可先用有限截断（文中示例 `n=-5..5`）
- 旋转半径参数示例 `r=0.1 GeV^-1`
- 因果约束：`omega * r < 1`

在本仓库实现中，上述三项应作为可配置参数并带输入校验。

## 8. 与现有主线的边界约定

- 本文档给出“文献一致的最小可运行核”，不绑定 legacy 旧脚本控制流。
- 迁移实现优先复用 `Models` 统一求解接口，不复制分散式脚本状态机。
- 高阶导数（如 `dP/domega`、高阶热导）列为增强阶段，不纳入首批最小核 DoD。

## 9. 代码映射建议

- 公式层：本页
- 代码层：`src/models/variants/rotation/`
- 测试层：
  - `tests/unit/models/test_rotation_model.jl`
  - `tests/integration/models/test_rotation_workflow_smoke.jl`
