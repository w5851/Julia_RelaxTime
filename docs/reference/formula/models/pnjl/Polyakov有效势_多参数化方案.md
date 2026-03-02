# PNJL：Polyakov 有效势多参数化方案（实现规范草案）

更新日期：2026-03-01

> 目的：为后续“Polyakov 势可切换实现”提供统一公式与参数来源说明，避免代码、配置、文档口径不一致。

---

## 1. 目标与范围

- 统一记录常用 Polyakov 势形式：
  - 多项式势（Ratti 参数系）
  - Fukushima 势
- 统一记录参数符号、量纲、典型取值和实现映射。
- 本文是“公式与参数口径文档”，不直接替代代码实现。

---

## 2. 多项式势（Ratti 常用形式）

```math
\frac{\mathcal{U}_{\text{poly}}(\Phi,\bar\Phi;T)}{T^4}
= -\frac{b_2(T)}{2}\,\bar\Phi\Phi
- \frac{b_3}{6}(\Phi^3+\bar\Phi^3)
+ \frac{b_4}{4}(\bar\Phi\Phi)^2
```

```math
b_2(T)=a_0+a_1\left(\frac{T_0}{T}\right)+a_2\left(\frac{T_0}{T}\right)^2+a_3\left(\frac{T_0}{T}\right)^3
```

### 2.1 常见参数组（文献常用）

- `T0 = 270 MeV`（纯胶子）
- `a0 = 6.75`
- `a1 = -1.95`
- `a2 = 2.625`
- `a3 = -7.44`
- `b3 = 0.75`
- `b4 = 7.5`

> 备注：含动态夸克时 `T0` 常下调到约 `200~215 MeV`。

---

## 3. Fukushima 势（常用写法）

```math
\frac{\mathcal{U}_{F}(\Phi,\bar\Phi;T)}{T^4}
= -b\,T\left[54 e^{-a/T}\Phi\bar\Phi+
\ln\Big(1-6\Phi\bar\Phi-3(\Phi\bar\Phi)^2+4(\Phi^3+\bar\Phi^3)\Big)\right]
```

### 3.1 常见参数组（文献常用）

- `a = 664 MeV`
- `b = (196 MeV)^3`

---

## 4. 与当前仓库实现的对应关系

当前主线热力学实现在：
- `src/models/pnjl/core/Thermodynamics.jl`

当前常量配置在：
- `src/Constants_PNJL.jl`

当前代码使用的是“对数势结构”口径，参数项包含：
- `T0_MeV, a0, a1, a2, b3, b4`

当前默认值（仓库现状）为：
- `T0_MeV = 210.0`
- `a0 = 3.51`
- `a1 = -2.47`
- `a2 = 15.2`
- `b3 = -1.75`
- `b4 = 7.555`

> 说明：这组值与“Ratti 多项式势参数组”不是同一组口径；后续实现多参数化切换时，需显式区分“势函数形式 + 参数组来源”。

---

## 5. 建议的实现落点（供开发任务使用）

### 5.1 配置层建议

建议新增：
- `polyakov.scheme = "log" | "poly" | "fukushima"`
- `polyakov.params.<...>` 对应所选 scheme 所需参数

### 5.2 最小接口建议

建议统一接口：

```text
polyakov_potential(T, Phi, Phibar; scheme, params) -> U
```

并在热力学主流程中仅调用该统一接口，避免在多处写分支。

### 5.3 验证建议

- `scheme=log` 且参数为当前默认值时，与现有主线结果一致（回归基线）。
- `scheme=poly` 与 `scheme=fukushima` 至少提供一组固定点/扫描点 smoke 对比结果。

---

## 6. 参考说明

- 本文参数为该类 PNJL 文献中常见口径整理，实际研究请以目标论文的具体参数表为准。
- 建议后续在实现文档中补“参数来源引用清单（论文-表号-版本）”。
