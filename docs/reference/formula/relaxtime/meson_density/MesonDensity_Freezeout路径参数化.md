# 介子数密度：freeze-out 路径参数化

本文固定介子数密度 charged / freeze-out 子链中，化学冻结线参数化的最小公式口径，以及它与文献中“实际扫描路径”的关系。

相关文档：

- [MesonDensity_稳定粒子与KPi比值.md](./MesonDensity_%E7%A8%B3%E5%AE%9A%E7%B2%92%E5%AD%90%E4%B8%8EKPi%E6%AF%94%E5%80%BC.md)
- [MesonDensity_BU相移公式.md](./MesonDensity_BU%E7%9B%B8%E7%A7%BB%E5%85%AC%E5%BC%8F.md)
- [2026-05-03_charged_freezeout_validation最小契约设计.md](../../../../dev/active/2026-05-03_charged_freezeout_validation%E6%9C%80%E5%B0%8F%E5%A5%91%E7%BA%A6%E8%AE%BE%E8%AE%A1.md)

## 1. 结论先行

当前本地文献库中，与 `K^\pm/\pi^\pm`、freeze-out 现象学最直接相关的几篇文献，**确实共享一条统计模型 chemical freeze-out 参数化基线**：

```math
T_{\mathrm{fo}}(\mu_B)
= a - b \mu_B^2 - c \mu_B^4,
```

```math
\mu_B(\sqrt{s_{NN}})
= \frac{d}{1 + e \sqrt{s_{NN}}}.
```

其中默认系数通常采用：

```math
a = 0.166~\mathrm{GeV},
\quad
b = 0.139~\mathrm{GeV}^{-1},
\quad
c = 0.053~\mathrm{GeV}^{-3},
```

```math
d = 1.308~\mathrm{GeV},
\quad
e = 0.273~\mathrm{GeV}^{-1}.
```

但当前本地相关项目也已经出现**同函数形式、不同系数组**的变体，例如：

```math
a = 0.164,\quad
b = 0.132848,\quad
c = 0.032,\quad
d = 1.336,\quad
e = 0.265.
```

这说明当前更稳定的统一结论应是：

1. **函数形式高度统一**；
2. **系数组并不天然唯一**；
3. 因此项目侧必须把“freeze-out parameterization family”与“default coefficient set”分开治理。

但要特别强调：

1. 这条式子只定义**统计模型 freeze-out 基线**；
2. 它**不自动等于**文献最终用于出图的模型扫描路径；
3. 在 `Blaschke:2016sqn`、`Blaschke:2019col`、`Blaschke:2020bzh` 一类工作中，常常还会进一步：
   - 改用 pseudo-critical line；
   - 采用 fixed `\mu/T` 映射；
   - 或用 critical-line + constant-`T` 的拼接代理路径。

因此，freeze-out 参数化与“最终扫描路径”必须在公式层明确拆开。

## 2. 引用口径与文献分工

当前这条链路里，文献角色应分成三类：

1. `Cleymans:2006qe` / 同系列 thermal-statistical freeze-out 拟合
   - 角色：**freeze-out 参数化的直接来源**
   - 作用：给出 `T_{\mathrm{fo}}(\mu_B)` 与 `\mu_B(\sqrt{s_{NN}})` 的经验拟合
   - 这是当前应固定进公式文档层的“默认 baseline”。
2. `Friesen:2019ojk`、`Blaschke:2021yml`
   - 角色：**当前项目最接近的 `K/\pi` 现象学使用范例**
   - 作用：明确说明 PNJL / EPNJL 中如何把 `K/\pi` 沿 freeze-out line 看作 `\sqrt{s_{NN}}` 的函数；
   - 这些文献直接复用了上面的 Cleymans 型参数化。
3. `Blaschke:2016sqn`、`Blaschke:2019col`、`Blaschke:2020bzh`
   - 角色：**freeze-out 基线如何被映射/替换为模型代理扫描路径的范例**
   - 作用：说明在 BU / horn 解释里，实际出图路径往往不是“原始 freeze-out 拟合曲线”本身，而是：
     - fixed `\mu/T`
     - pseudo-critical line
     - 或拼接扫描线。

因此当前公式治理建议是：

- **公式 baseline**：固定 Cleymans 型 freeze-out 参数化；
- **模型扫描路径**：作为独立上层对象，不直接塞进同一条公式定义里。

## 3. 默认 freeze-out 参数化

### 3.1 `T_{\mathrm{fo}}(\mu_B)`

```math
T_{\mathrm{fo}}(\mu_B)
= a - b \mu_B^2 - c \mu_B^4.
```

这里：

- `T_{\mathrm{fo}}` 与 `\mu_B` 都以 `\mathrm{GeV}` 计；
- `a,b,c` 使用上文默认系数。

### 3.2 `\mu_B(\sqrt{s_{NN}})`

```math
\mu_B(\sqrt{s_{NN}})
= \frac{d}{1 + e \sqrt{s_{NN}}}.
```

这里：

- `\sqrt{s_{NN}}` 以 `\mathrm{GeV}` 计；
- `\mu_B` 以 `\mathrm{GeV}` 计；
- `d,e` 使用上文默认系数。

### 3.3 组合映射

把二者组合后，可得：

```math
\sqrt{s_{NN}}
\rightarrow
\mu_B(\sqrt{s_{NN}})
\rightarrow
T_{\mathrm{fo}}(\mu_B(\sqrt{s_{NN}})).
```

这就是当前最小的 freeze-out 路径 baseline。

## 4. 文献中的实际使用方式

### 4.1 `Friesen:2019ojk`

这篇文献明确说明：

1. 实验横轴是 `\sqrt{s_{NN}}`；
2. 有效模型本身并不直接以 `\sqrt{s_{NN}}` 为输入；
3. 因此要先用 Cleymans 型参数化把 `\sqrt{s_{NN}}` 映射成 freeze-out 下的 `(T,\mu_B)`。

这属于“直接把 freeze-out 参数化当 baseline 输入”的做法。

### 4.2 `Blaschke:2021yml`

这篇文献同样明确采用：

1. freeze-out baryon chemical potential 参数化；
2. 沿 freeze-out line 考察 `K/\pi` 比值；
3. 再叠加 `\mu_\pi`、`\mu_s` 的现象学口径。

因此它支持把 Cleymans 型参数化作为 charged / freeze-out 子链的默认基线。

### 4.3 `Blaschke:2016sqn` 与 `Blaschke:2019col`

这两篇文献的重要区别在于：

1. 它们并不满足于“直接沿统计模型 freeze-out 曲线扫”；
2. 而是引入了：
   - fixed `\mu/T` 映射；
   - pseudo-critical line 代理；
   - 或 critical-line 与 constant-`T` 拼接路径。

这意味着：

- 文献中的 Figure 4 结果曲线，不能简单等同于“原始 freeze-out 参数化曲线跑出来的结果”；
- 若要严格复现，必须额外定义“路径代理规则”。

## 5. 当前项目中应如何落地

当前建议把 freeze-out 相关对象拆成两层：

### 5.1 第一层：参数化 baseline

只固定：

```math
\sqrt{s_{NN}} \mapsto \mu_B \mapsto T_{\mathrm{fo}}.
```

这层的职责是：

1. 为 charged / freeze-out workflow 提供统一输入；
2. 让 Figure 1 中的 freeze-out 背景层可由项目自绘；
3. 提供最小一致的 `\sqrt{s_{NN}} \leftrightarrow (T,\mu_B)` 数据表。

### 5.2 第二层：模型扫描路径

单独定义：

1. `actual_freezeout`
2. `critical_proxy`
3. `fixed_mu_over_T_mapping`
4. `stitched_critical_plus_constT`

这些属于“路径策略”，不应混在 baseline 公式定义中。

## 6. 与 solver 设计的关系

如果后续在 solver 层加入 freeze-out 模式，更合理的对象应是：

1. **baseline 路径模式**
   - 输入：`x = \sqrt{s_{NN}}` 或等价一维路径参数
   - 内部映射：`x -> (T,\mu_B)`
   - 底层求解：复用 `FixedMu` 主链
2. **代理路径模式**
   - 若文献需要 pseudo-critical / stitched path
   - 则应在 baseline 之上显式命名为另一种 path spec

不建议把“所有 freeze-out 及其代理路径”直接揉成一个不透明的新公式对象。

## 7. 当前统一性验证

当前仓库外的 `D:\Desktop\PNJL_Simulation` 中，已经存在同一组基线公式定义：

- `src/Gas_Liquid/Advanced_KappaforS.jl`
- `examples/Gas_Liquid/scan_tmu_for_s.jl`

其中明确写出：

```julia
T_from_μB(μB; a=0.166, b=0.139, c=0.053)
μB_from_S(s; d=1.308, e=0.273)
```

与此同时，`examples/Gas_Liquid/scan_tmu_for_s.jl` 的主调用又显式覆盖为：

```julia
main(["2.0", "200.0", "0.1"]; a=0.164, b=0.132848, c=0.032, d=1.336, e=0.265)
```

因此当前统一性验证给出的更准确结论是：

1. `PNJL_Simulation` 与本仓库都适合采用同一**函数族**；
2. 但应允许：
   - `default_cleymans_2006_like`
   - `updated_project_fit`
   这样的显式系数组标签；
3. 在未进一步追溯替代系数组文献来源前，本仓库公式层默认仍应优先固定 Cleymans 型默认组，而把替代组视为可选 profile。

## 8. 当前结论

当前应固定的 freeze-out 公式层结论是：

1. **有统一的 baseline 参数化**
   - 即 Cleymans 型 `T(\mu_B)` 与 `\mu_B(\sqrt{s_{NN}})` 函数形式；
2. **没有统一的最终扫描路径**
   - BU / horn 文献常在 baseline 之上再做路径代理；
3. **实现层应显式区分 baseline 与 path strategy**
   - 这既适合公式治理，也适合 solver / workflow 设计。
