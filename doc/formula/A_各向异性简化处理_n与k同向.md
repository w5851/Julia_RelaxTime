# [有限温度有限密度场论中的单圈积分A（one-line integral）各向异性简化处理] - 实现文档

## 1. 公式标识
- **程序实现**: `A_aniso`
- **对应文件**: `../../src/relaxtime/OneLoopIntegrals.jl`

## 2. 各向异性修正原理

在考虑各向异性夸克分布函数的一阶修正后，原A积分需要引入PNJL模型中的各向异性分布函数。假设各向异性方向 $\mathbf{n}$ 和动量 $\vec{k}$ 均沿z轴。

### 2.1 原始各向同性表达式
```math
A(m, \mu, \beta, \Lambda) = 4\int_0^{\Lambda} p^2 dp \frac{1}{E} \left(f(E - \mu) + f(E + \mu) - 1\right)
```

### 2.2 各向异性分布函数修正

在PNJL模型中，各向异性分布函数的一阶修正为：

- **对于夸克分布**：
  \[
  f_q^{\text{aniso}} = f_q^+ + \frac{\xi \beta}{2} \frac{p^2 x^2}{E} \cdot G(u)
  \]
  其中 $u = \beta(E - \mu)$，$G(u) = \frac{N'(u)D(u) - N(u)D'(u)}{[D(u)]^2}$

- **对于反夸克分布**：
  \[
  f_{\bar{q}}^{\text{aniso}} = f_q^- + \frac{\xi \beta}{2} \frac{p^2 x^2}{E} \cdot H(y)
  \]
  其中 $y = \beta(E + \mu)$，$H(y) = \frac{C'(y)D(y) - C(y)D'(y)}{[D(y)]^2}$

其中Polyakov圈相关函数定义与B0文档中相同。

## 3. 各向异性修正后的A积分

### 3.1 分解形式
将A积分分解为各向同性部分和各向异性修正部分：
\[
A = A^{\text{iso}} + A^{\text{aniso}}
\]

### 3.2 各向同性部分
使用PNJL的各向同性分布函数：
\[
A^{\text{iso}} = 4\int_0^{\Lambda} p^2 dp \frac{1}{E} \left[f_q^+(u) + f_q^-(y) - 1\right]
\]
其中：
- $u = \beta(E - \mu)$
- $y = \beta(E + \mu)$

### 3.3 各向异性修正部分

由于各向异性修正项包含 $x^2 = \cos^2\theta$，需要对立体角进行积分：
\[
A^{\text{aniso}} = 4\int_0^{\Lambda} p^2 dp \frac{1}{E} \left[\delta f_q + \delta f_{\bar{q}}\right]
\]

其中各向异性修正项为：
\[
\delta f_q = \frac{\xi \beta}{2} \frac{p^2}{E} G(u) \cdot \frac{1}{4\pi} \int d\Omega\, x^2
\]
\[
\delta f_{\bar{q}} = \frac{\xi \beta}{2} \frac{p^2}{E} H(y) \cdot \frac{1}{4\pi} \int d\Omega\, x^2
\]

对立体角积分：
\[
\frac{1}{4\pi} \int d\Omega\, x^2 = \frac{1}{4\pi} \int_0^{2\pi} d\phi \int_{-1}^1 x^2 dx = \frac{1}{3}
\]

因此各向异性修正部分简化为：
\[
A^{\text{aniso}} = \frac{4\xi \beta}{3} \int_0^{\Lambda} p^4 dp \frac{1}{E^2} \left[G(u) + H(y)\right]
\]

### 3.4 最终表达式

**完整的各向异性A积分**：
\[
A(m, \mu, T, \Lambda, \xi, \Phi, \bar{\Phi}) = 4\int_0^{\Lambda} p^2 dp \frac{1}{E} \left[f_q^+(u) + f_q^-(y) - 1\right] + \frac{4\xi \beta}{3} \int_0^{\Lambda} p^4 dp \frac{1}{E^2} \left[G(u) + H(y)\right]
\]

## 4. 常数项的处理

对于常数项"-1"的积分，使用解析结果：
\[
\int_0^{\Lambda} \frac{p^2}{E} dp = \frac{\Lambda}{2} \sqrt{\Lambda^2 + m^2} - \frac{m^2}{2} \ln\left(\frac{\Lambda + \sqrt{\Lambda^2 + m^2}}{m}\right)
\]

## 5. 参数说明表

| 参数 | 符号 | 类型 | 单位 | 物理意义 | 取值范围 |
|------|------|------|------|----------|----------|
| m | m | 输入 | fm⁻¹ | 夸克质量 | [0.005, 2.5] fm⁻¹ |
| μ | μ | 输入 | fm⁻¹ | 夸克化学势 | [-5, 5] fm⁻¹ |
| T | T | 输入 | fm⁻¹ | 温度 | [0.005, 2.0] fm⁻¹ |
| Λ | Λ | 输入 | fm⁻¹ | 动量截断 | [2.0, 10.0] fm⁻¹ |
| ξ | ξ | 输入 | 无量纲 | 各向异性参数 | [-0.5, 0.5] |
| Φ | Φ | 输入 | 无量纲 | Polyakov圈 | [0, 1] |
| $\bar{\Phi}$ | $\bar{\Phi}$ | 输入 | 无量纲 | 共轭Polyakov圈 | [0, 1] |
| A | A | 输出 | fm⁻² | 单圈积分值 | 实数 |

## 6. 实现细节

### 6.1 数值积分策略
- 对分布函数部分使用高斯-勒让德积分
- 对常数项使用解析表达式
- 各向异性修正项单独积分

### 6.2 代码结构
```julia
function A_aniso(m::Float64, μ::Float64, T::Float64, Λ::Float64, 
                ξ::Float64, Φ::Float64, Φbar::Float64;
                n_points::Int=64)
    
    # 1. 计算各向同性部分
    A_iso = compute_A_iso(m, μ, T, Λ, Φ, Φbar, n_points)
    
    # 2. 计算各向异性修正部分
    A_aniso_corr = compute_A_aniso_correction(m, μ, T, Λ, ξ, Φ, Φbar, n_points)
    
    return A_iso + A_aniso_corr
end
```

### 6.3 分布函数实现
- `f_q_plus(u, Φ, Φbar)`: 夸克分布函数
- `f_q_minus(y, Φ, Φbar)`: 反夸克分布函数  
- `G_function(u, Φ, Φbar)`: 各向异性修正函数G
- `H_function(y, Φ, Φbar)`: 各向异性修正函数H

## 7. 极限情况验证

### 7.1 各向同性极限
当 $\xi = 0$ 时，$A^{\text{aniso}} = 0$，回归标准PNJL模型。

### 7.2 高温极限
当 $T \to \infty$，$\Phi, \bar{\Phi} \to 1$，分布函数退化为费米分布。

### 7.3 零化学势情况
当 $\mu = 0$ 时，$f_q^+ = f_q^-$，表达式对称化。

## 8. 物理意义

各向异性修正反映了动量空间分布的各向异性对单圈积分的贡献：
- **实部**：影响能谱和有效质量
- **虚部**：在A积分中通常为零（无传播子极点）

这种修正在研究重离子碰撞产生的各向异性夸克胶子等离子体时尤为重要。

## 9. 与B0积分的对比

| 特征 | A积分 | B0积分 |
|------|-------|--------|
| 传播子数 | 1个 | 2个 |
| 角度依赖 | 各向同性+各向异性修正 | 显式角度积分 |
| 虚部 | 通常为零 | 非零（有传播子极点）|
| 发散性 | 紫外发散 | 紫外发散 |
| 各向异性修正形式 | 乘以1/3因子 | 复杂角度积分 |