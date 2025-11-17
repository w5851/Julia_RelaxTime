# 总传播子计算 - 实现文档

## 1. 公式标识
- **程序实现**: `total_propagator`
- **对应文件**: `src/relaxtime/TotalPropagator.jl`（待实现）

## 2. 物理意义
- **物理背景**: 在PNJL模型中，夸克弹性散射通过多种介子交换实现，同一散射道中不同介子的贡献需按味因子加权求和
- **在计算链中的作用**: 将单个介子传播子组合为特定散射过程的总传播子，用于计算散射矩阵元
- **相关物理量**: 介子传播子、味因子、散射道(s/t/u)

## 3. 数学表达式

### 3.1 总传播子公式

对于给定散射道，总传播子的一般形式：

```math
\mathcal{D}_{\text{channel}}^{P,S} = T_1 \times \sum_{\text{meson}} \mathcal{D}_{\text{meson}}^{P(S)} \times T_2
```

其中：
- $\mathcal{D}_{\text{channel}}^{P,S}$：包含味因子的总传播子（复数，单位：fm²）
- $T_1, T_2$：两个顶角的味因子（实数，见下表5.3）
- $\mathcal{D}_{\text{meson}}^{P(S)}$：单个介子传播子（调用`MesonPropagator`计算）
- $\sum_{\text{meson}}$：对该散射道允许的所有介子求和

### 3.2 混合介子的展开形式

对于η/η'和σ/σ'混合介子，已在`meson_propagator_mixed`中实现，总传播子直接调用该函数：

```math
\mathcal{D} = 2\frac{\det K}{\det M} J^T M J'
```

展开形式（已在MesonPropagator.jl中实现）：
```math
\mathcal{D} = 2\frac{\det K}{M_{00}M_{88}-M_{08}^2}(M_{00}\bar\psi\lambda_0\psi\cdot\bar\psi'\lambda_0\psi' + M_{08}\bar\psi\lambda_0\psi\cdot\bar\psi'\lambda_8\psi' + M_{08}\bar\psi\lambda_8\psi\cdot\bar\psi'\lambda_0\psi' + M_{88}\bar\psi\lambda_8\psi\cdot\bar\psi'\lambda_8\psi')
```

**注意**：混合介子传播子已包含完整的味结构，**不需要再乘以额外的味因子**。

## 4. 味因子查询表

### 表5.3 各散射道的味因子

| 夸克类型 | u | d | s |
|----------|---|---|---|
| **u** | 1 | √2 | √2 |
| **d** | √2 | -1 | √2 |
| **s** | √2 | √2 | 2 |

**使用规则**：
- 行对应第一个夸克，列对应第二个夸克
- 表格对称：$T(i,j) = T(j,i)$
- 正反夸克使用相同味因子

**各散射道的T₁和T₂确定规则**：
- **t道**：T₁对应入射夸克1→出射夸克3，T₂对应入射夸克2→出射夸克4
- **u道**：T₁对应入射夸克1→出射夸克4，T₂对应入射夸克2→出射夸克3
- **s道**：T₁对应入射夸克1-反入射夸克2，T₂对应出射夸克3-反出射夸克4

**示例**：
```julia
# us → us 过程的t道
# T₁ = T(u,u) = 1，T₂ = T(s,s) = 2
# 总传播子 = 1 × (D_K + D_σ_K) × 2

# ud → ud 过程的u道  
# T₁ = T(u,d) = √2，T₂ = T(d,u) = √2
# 总传播子 = √2 × (介子和) × √2 = 2 × (介子和)
```

## 5. 计算流程

### 5.1 一般介子传播子求和（π, K, σ_π, σ_K）
1. **确定允许的介子列表**：根据散射过程和道类型（参考文献表5.1和表5.2）
2. **计算单个传播子**：对每个介子调用`meson_propagator_simple`
3. **应用味因子**：$\mathcal{D}_{\text{total}} = T_1 \times \sum_i \mathcal{D}_i \times T_2$

### 5.2 混合介子传播子（η/η', σ/σ'）
1. **直接调用**：使用`meson_propagator_mixed`函数
2. **注意**：混合介子已包含完整味结构，无需额外乘味因子

## 6. 依赖关系

### 6.1 输入依赖
- `MesonPropagator.meson_propagator_simple`：一般介子传播子
- `MesonPropagator.meson_propagator_mixed`：混合介子传播子
- 极化函数Π：通过`Polarization`或`PolarizationAniso`模块计算
- K系数：通过`EffectiveCouplings`模块计算

### 6.2 输出用途
- 散射矩阵元计算
- 驰豫时间计算
- 输运系数计算

## 7. 实现要点

1. **复数运算**：所有传播子为ComplexF64类型
2. **味因子对称性**：利用T(i,j) = T(j,i)简化计算
3. **混合介子特殊处理**：η/η'和σ/σ'不能简单与一般介子求和
4. **性能优化策略**：
   - K系数可在批量计算中复用
   - **极化函数缓存**：使用`PolarizationCache`模块避免重复计算相同参数的极化函数（详见`doc/formula/PolarizationCache_缓存使用文档.md`）
   - 典型收益：缓存命中率30%-70%，可节省约50%计算时间