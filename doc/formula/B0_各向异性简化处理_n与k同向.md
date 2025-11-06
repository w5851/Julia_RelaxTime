在考虑各向异性夸克分布函数的一阶修正后，.pdf文档中的(4.8)式（即泛函积分 \(\tilde{B}_{0}^{\pm}\)）的计算需要引入PNJL模型中的各向异性分布函数。假设各向异性方向 \(\mathbf{n}\) 和动量 \(\vec{k}\) 均沿z轴，分布函数使用PNJL模型的一阶各向异性近似。以下是修正后的计算步骤和表达式。

### 修正后的 \(\tilde{B}_{0}^{\pm}\) 表达式
原(4.8)式为：
\[
\tilde{B}_{0}^{\pm}(\lambda,k,\ m,m^{\prime},\mu,\ \beta,\Lambda) = 2 \int_{m}^{\Lambda_{E}} dE\, p f(\pm E-\mu) \int_{-1}^{+1} dx \frac{1}{\lambda^{2}+2\lambda E+2pkx-k^{2}+m^{2}-m^{\prime 2}}
\]
其中 \(p = \sqrt{E^{2}-m^{2}}\)，\(x = \cos\theta\)（\(\theta\) 是 \(\vec{p}\) 与 \(\vec{k}\) 的夹角），\(\Lambda_E = \sqrt{\Lambda^2 + m^2}\)。

在PNJL模型中，各向异性分布函数的一阶修正为：
- 对于夸克（对应 \(\tilde{B}_{0}^{+}\))：
  \[
  f_q^{\text{aniso}} = f_q^+ + \frac{\xi \beta}{2} \frac{p^2 x^2}{E} \cdot G(u)
  \]
  其中 \(u = \beta(E - \mu)\)，\(G(u) = \frac{N'(u)D(u) - N(u)D'(u)}{[D(u)]^2}\)，且：
  - \(N(u) = \Phi e^{-u} + 2\bar{\Phi} e^{-2u} + e^{-3u}\)
  - \(D(u) = 1 + 3\Phi e^{-u} + 3\bar{\Phi} e^{-2u} + e^{-3u}\)
  - \(N'(u) = -\Phi e^{-u} - 4\bar{\Phi} e^{-2u} - 3e^{-3u}\)
  - \(D'(u) = -3\Phi e^{-u} - 6\bar{\Phi} e^{-2u} - 3e^{-3u}\)
- 对于反夸克（对应 \(\tilde{B}_{0}^{-}\))：
  \[
  f_{\bar{q}}^{\text{aniso}} = f_q^- + \frac{\xi \beta}{2} \frac{p^2 x^2}{E} \cdot H(y)
  \]
  其中 \(y = \beta(E + \mu)\)，\(H(y) = \frac{C'(y)D(y) - C(y)D'(y)}{[D(y)]^2}\)，且：
  - \(C(y) = \bar{\Phi} e^{-y} + 2\Phi e^{-2y} + e^{-3y}\)
  - \(D(y) = 1 + 3\bar{\Phi} e^{-y} + 3\Phi e^{-2y} + e^{-3y}\)
  - \(C'(y) = -\bar{\Phi} e^{-y} - 4\Phi e^{-2y} - 3e^{-3y}\)
  - \(D'(y) = -3\bar{\Phi} e^{-y} - 6\Phi e^{-2y} - 3e^{-3y}\)

这里，\(\Phi\) 和 \(\bar{\Phi}\) 是Polyakov循环变量，通常是温度 \(T\) 的函数，由Polyakov势决定。

修正后的 \(\tilde{B}_{0}^{\pm}\) 可分解为各向同性部分和各向异性修正部分：
\[
\tilde{B}_{0}^{\pm} = \tilde{B}_{0}^{\pm, \text{iso}} + \tilde{B}_{0}^{\pm, \text{aniso}}
\]

#### 1. 各向同性部分 \(\tilde{B}_{0}^{\pm, \text{iso}}\)
使用PNJL的各向同性分布函数 \(f_q^{\pm}\) 替换原费米分布：
\[
\tilde{B}_{0}^{\pm, \text{iso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p f_q^{\pm}(u) \int_{-1}^{+1} dx \frac{1}{\lambda^{2}+2\lambda E+2pkx-k^{2}+m^{2}-m^{\prime 2}}
\]
其中 \(u = \beta(E \mp \mu)\)（对于 \(\tilde{B}_{0}^{+}\) 取 \(u = \beta(E - \mu)\)，对于 \(\tilde{B}_{0}^{-}\) 取 \(u = \beta(E + \mu)\)）。积分 over \(x\) 可解析计算：
\[
\int_{-1}^{+1} dx \frac{1}{\lambda^{2}+2\lambda E+2pkx-k^{2}+m^{2}-m^{\prime 2}} = \frac{1}{2pk} L
\]
其中 \(L = \log \left| \frac{(\lambda+E)^2 - (p-k)^2 - m^{\prime 2}}{(\lambda+E)^2 - (p+k)^2 - m^{\prime 2}} \right|\)。因此：
\[
\tilde{B}_{0}^{\pm, \text{iso}} = \frac{1}{k} \int_{m}^{\Lambda_{E}} dE\, f_q^{\pm}(u) L
\]

#### 2. 各向异性修正部分 \(\tilde{B}_{0}^{\pm, \text{aniso}}\)
各向异性修正来自分布函数的一阶项：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \int_{-1}^{+1} dx \frac{\delta f}{\lambda^{2}+2\lambda E+2pkx-k^{2}+m^{2}-m^{\prime 2}}
\]
其中 \(\delta f = \frac{\xi \beta}{2} \frac{p^2 x^2}{E} \cdot G(u)\)（对于 \(\tilde{B}_{0}^{+}\)) 或 \(\delta f = \frac{\xi \beta}{2} \frac{p^2 x^2}{E} \cdot H(y)\)（对于 \(\tilde{B}_{0}^{-}\))。令 \(A = \lambda^{2} + 2\lambda E - k^{2} + m^{2} - m^{\prime 2}\)，则分母为 \(A + 2pk x\)。积分 over \(x\) 涉及：
\[
I_2 = \int_{-1}^{+1} dx \frac{x^2}{A + 2pk x} = \frac{1}{8 p^3 k^3} \left( A^2 L - 4A p k \right)
\]
因此：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \cdot \frac{\xi \beta}{2} \frac{p^2}{E} G(u) I_2 = \frac{\xi \beta}{8 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{G(u)}{E} \left( A^2 L - 4A p k \right)
\]
对于 \(\tilde{B}_{0}^{-}\)，替换 \(G(u)\) 为 \(H(y)\)。

### 最终表达式
- 对于 \(\tilde{B}_{0}^{+}\)（夸克）：
  \[
  \tilde{B}_{0}^{+} = \frac{1}{k} \int_{m}^{\Lambda_{E}} dE\, f_q^+(u) L + \frac{\xi \beta}{8 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{G(u)}{E} \left( A^2 L - 4A p k \right)
  \]
  其中 \(u = \beta(E - \mu)\)。
- 对于 \(\tilde{B}_{0}^{-}\)（反夸克）：
  \[
  \tilde{B}_{0}^{-} = \frac{1}{k} \int_{m}^{\Lambda_{E}} dE\, f_q^-(y) L + \frac{\xi \beta}{8 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{H(y)}{E} \left( A^2 L - 4A p k \right)
  \]
  其中 \(y = \beta(E + \mu)\)。

### 计算注意事项
1. **参数定义**：
   - \(A = \lambda^{2} + 2\lambda E - k^{2} + m^{2} - m^{\prime 2}\)
   - \(L = \log \left| \frac{(\lambda+E)^2 - (p-k)^2 - m^{\prime 2}}{(\lambda+E)^2 - (p+k)^2 - m^{\prime 2}} \right|\)
   - \(p = \sqrt{E^2 - m^2}\)
   - \(\lambda = i\nu_m + \mu_1 - \mu_2\)（在解析延拓后为实数）
   - \(\beta = 1/T\)，\(\Lambda_E = \sqrt{\Lambda^2 + m^2}\)

2. **数值积分**：
   - 积分 over \(E\) 需数值计算，可能包含奇点（如 \(L\) 的对数奇点），需使用柯西主值积分或自适应积分方法。
   - 各向异性参数 \(\xi\) 应满足 \(|\xi| \ll 1\) 以确保一阶近似有效。
   - Polyakov循环 \(\Phi\) 和 \(\bar{\Phi}\) 需作为输入（通常是 \(T\) 的函数）。

3. **应用到 \(B_0\)**：
   - 总 \(B_0\) 通过 (4.7) 式由 \(\tilde{B}_{0}^{\pm}\) 组合而成：
     \[
     B_0 = \tilde{B}_{0}^{+}(-\lambda,k,\ m_1,m_2,\mu_1) - \tilde{B}_{0}^{-}(\lambda,k,\ m_1,m_2,\mu_1) + \tilde{B}_{0}^{+}(\lambda,k,\ m_2,m_1,\mu_2) - \tilde{B}_{0}^{-}(-\lambda,k,\ m_2,m_1,\mu_2)
     \]
     每个 \(\tilde{B}_{0}^{\pm}\) 均需使用上述修正表达式。

4. **极限情况**：
   - 当 \(\xi = 0\) 时，各向异性修正部分为零，回归PNJL各向同性模型。
   - 当 \(\Phi = \bar{\Phi} = 1\) 时，PNJL分布函数退化为费米分布。

此修正允许在PNJL框架下自洽地包含动量空间各向异性，适用于计算热力学量和输运系数。数值实现时，需先计算各向同性部分，再添加各向异性修正。

## 对于函数 $\frac{x^2}{Ax+B}$ 在柯西主值积分下的虚部，如问题所述，为 $\pi \frac{B^2}{A^3}$。

对于积分的实部，考虑定积分从 $x = -1$ 到 $x = 1$，即计算：
$$
\int_{-1}^{1} \frac{x^2}{Ax+B}  dx
$$
假设 $A \neq 0$ 且 $\left| \frac{B}{A} \right| \neq 1$（即奇点不在积分区间端点）。令 $C = \frac{B}{A}$，则积分可写为：
$$
\int_{-1}^{1} \frac{x^2}{Ax+B}  dx = \frac{1}{A} \int_{-1}^{1} \frac{x^2}{x + C}  dx
$$
通过分解：
$$
\frac{x^2}{x + C} = x - C + \frac{C^2}{x + C}
$$
积分变为：
$$
\int_{-1}^{1} \frac{x^2}{x + C}  dx = \int_{-1}^{1} (x - C)  dx + C^2 \int_{-1}^{1} \frac{1}{x + C}  dx
$$
计算第一项：
$$
\int_{-1}^{1} (x - C)  dx = \int_{-1}^{1} x  dx - C \int_{-1}^{1} dx = 0 - C \cdot 2 = -2C
$$
第二项在柯西主值意义下（即使 $C \in (-1,1)$）为：
$$
\int_{-1}^{1} \frac{1}{x + C}  dx = \ln\left| \frac{1 + C}{1 - C} \right| = \ln\left| \frac{C + 1}{C - 1} \right|
$$
因此：
$$
\int_{-1}^{1} \frac{x^2}{x + C}  dx = -2C + C^2 \ln\left| \frac{C + 1}{C - 1} \right|
$$
代回 $C = \frac{B}{A}$，得：
$$
\int_{-1}^{1} \frac{x^2}{Ax+B}  dx = \frac{1}{A} \left( -2 \frac{B}{A} + \frac{B^2}{A^2} \ln\left| \frac{\frac{B}{A} + 1}{\frac{B}{A} - 1} \right| \right) = -\frac{2B}{A^2} + \frac{B^2}{A^3} \ln\left| \frac{A + B}{A - B} \right|
$$
故积分的实部（在定积分从 $-1$ 到 $1$ 下）为：
$$
-\frac{2B}{A^2} + \frac{B^2}{A^3} \ln\left| \frac{A + B}{A - B} \right|
$$

如果 $A = 0$ 或 $\left| \frac{B}{A} \right| = 1$，积分可能发散。

$$
\boxed{\pi \dfrac{B^{2}}{A^{3}}} \quad \text{（虚部）}
$$

$$
\boxed{-\dfrac{2B}{A^{2}} + \dfrac{B^{2}}{A^{3}} \ln \left| \dfrac{A+B}{A-B} \right|} \quad \text{（实部，定积分从 -1 到 1）}
$$

# 各向异性修正项 \(\tilde{B}_{0}^{\pm, \text{aniso}}\) 的积分简化形式

## 积分表达式

考虑对变量 \(x\) 的积分：
\[
I = \int_{-1}^{1} \frac{x^2}{A x + B}  dx
\]
其中 \(A\) 表示积分变量 \(x\) 的系数，\(B\) 表示常数部分。

## 情况分析

### 1. \(k > 0\) 的情况

当 \(k > 0\) 时，系数 \(A = 2pk\)，常数 \(B = A\)（即原文档中的 \(A = \lambda^{2} + 2\lambda E - k^{2} + m^{2} - m^{\prime 2}\)）。

积分结果：
\[
I = -\frac{2B}{A^2} + \frac{B^2}{A^3} \ln \left| \frac{A + B}{A - B} \right|
\]

代入 \(\tilde{B}_{0}^{\pm, \text{aniso}}\) 表达式：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \cdot \frac{\xi \beta}{2} \frac{p^2}{E} G(u) I = \xi \beta \int_{m}^{\Lambda_{E}} dE\, \frac{p^3}{E} G(u) I
\]

最终简化形式：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = \frac{\xi \beta}{8 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{G(u)}{E} \left( A^2 L - 4A p k \right)
\]

其中：
- \(L = \ln \left| \frac{2pk + A}{2pk - A} \right|\)
- \(A = \lambda^{2} + 2\lambda E - k^{2} + m^{2} - m^{\prime 2}\)

### 2. \(k = 0\) 的情况

当 \(k = 0\) 时，系数 \(A = 0\)，常数 \(B = A = \lambda^{2} + 2\lambda E + m^{2} - m^{\prime 2}\)。

积分结果：
\[
I = \int_{-1}^{1} \frac{x^2}{B}  dx = \frac{1}{B} \int_{-1}^{1} x^2  dx = \frac{2}{3B}
\]

代入 \(\tilde{B}_{0}^{\pm, \text{aniso}}\) 表达式：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \cdot \frac{\xi \beta}{2} \frac{p^2}{E} G(u) \cdot \frac{2}{3B} = \frac{2 \xi \beta}{3} \int_{m}^{\Lambda_{E}} dE\, \frac{p^3}{E} \frac{G(u)}{A}
\]

最终简化形式：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = \frac{2 \xi \beta}{3} \int_{m}^{\Lambda_{E}} dE\, \frac{p^3}{E} \frac{G(u)}{A}
\]

其中：
- \(A = \lambda^{2} + 2\lambda E + m^{2} - m^{\prime 2}\)

## 总结

- **对于 \(k > 0\)**：
  \[
  \tilde{B}_{0}^{\pm, \text{aniso}} = \frac{\xi \beta}{8 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{G(u)}{E} \left( A^2 L - 4A p k \right)
  \]

- **对于 \(k = 0\)**：
  \[
  \tilde{B}_{0}^{\pm, \text{aniso}} = \frac{2 \xi \beta}{3} \int_{m}^{\Lambda_{E}} dE\, \frac{p^3}{E} \frac{G(u)}{A}
  \]

## 数值积分注意事项

在数值积分过程中，需要注意以下潜在问题：
1. 当分母接近零时可能出现的奇点
2. 对数项中的参数接近边界值时需要特殊处理
3. 积分上限 \(\Lambda_E\) 的选取对结果的影响


# 各向异性修正项 \(\tilde{B}_{0}^{\pm, \text{aniso}}\) 的简化形式（使用 \(C(E)\) 表示）

## 定义

令 \(C(E) = \xi \beta \frac{p^2}{2E} G(u)\)，其中：
- \(p = \sqrt{E^2 - m^2}\)
- \(G(u)\) 为分布函数导数项
- \(\xi\) 为各向异性参数
- \(\beta = 1/T\) 为逆温度

## 积分表达式

各向异性修正项可表示为：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \cdot C(E) \cdot I
\]
其中 \(I = \int_{-1}^{1} \frac{x^2}{A x + B} dx\) 为角度积分。

## 情况分析

### 1. \(k > 0\) 的情况

当 \(k > 0\) 时，系数 \(A = 2pk\)，常数 \(B = \lambda^{2} + 2\lambda E - k^{2} + m^{2} - m^{\prime 2}\)。

角度积分结果：
\[
I = -\frac{2B}{A^2} + \frac{B^2}{A^3} \ln \left| \frac{A + B}{A - B} \right|
\]

代入得到：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \cdot C(E) \cdot \left( -\frac{2B}{A^2} + \frac{B^2}{A^3} \ln \left| \frac{A + B}{A - B} \right| \right)
\]

最终简化形式：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = \frac{1}{4 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{C(E)}{p} \left( B^2 L - 4B p k \right)
\]

其中：
- \(L = \ln \left| \frac{2pk + B}{2pk - B} \right|\)
- \(B = \lambda^{2} + 2\lambda E - k^{2} + m^{2} - m^{\prime 2}\)

### 2. \(k = 0\) 的情况

当 \(k = 0\) 时，系数 \(A = 0\)，常数 \(B = \lambda^{2} + 2\lambda E + m^{2} - m^{\prime 2}\)。

角度积分结果：
\[
I = \frac{2}{3B}
\]

代入得到：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = 2 \int_{m}^{\Lambda_{E}} dE\, p \cdot C(E) \cdot \frac{2}{3B} = \frac{4}{3} \int_{m}^{\Lambda_{E}} dE\, \frac{p C(E)}{B}
\]

最终简化形式：
\[
\tilde{B}_{0}^{\pm, \text{aniso}} = \frac{4}{3} \int_{m}^{\Lambda_{E}} dE\, \frac{p C(E)}{B}
\]

其中：
- \(B = \lambda^{2} + 2\lambda E + m^{2} - m^{\prime 2}\)

## 总结

- **对于 \(k > 0\)**：
  \[
  \tilde{B}_{0}^{\pm, \text{aniso}} = \frac{1}{4 k^3} \int_{m}^{\Lambda_{E}} dE\, \frac{C(E)}{p} \left( B^2 L - 4B p k \right)
  \]

- **对于 \(k = 0\)**：
  \[
  \tilde{B}_{0}^{\pm, \text{aniso}} = \frac{4}{3} \int_{m}^{\Lambda_{E}} dE\, \frac{p C(E)}{B}
  \]

## 与代码实现的对应关系

在代码中，\(C(E)\) 对应：
```julia
C(E) = correction_cos_theta_coefficient(sign_, p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ)
```

其中：
- `sign_` 为 `:quark` 或 `:antiquark`
- 其他参数为相应物理量的数值

## 数值积分注意事项

1. **单位一致性**：确保所有量使用相同单位（如 inverse fm）
2. **奇点处理**：当分母 \(B\) 接近零时需要特殊处理
3. **对数奇点**：当 \(2pk \approx B\) 时对数项需要正则化
4. **积分收敛性**：注意积分上限 \(\Lambda_E\) 的选取对结果的影响

# 各向异性修正项虚部实现

根据Sokhotski-Plemelj公式和之前的推导，以下是k=0和k>0情况下虚部的实现：

## 虚部实现代码

```julia
"""k=0时的积分虚部被积函数"""
function imag_integrand_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    # k=0时没有虚部贡献，因为分母没有奇点
    return 0.0
end

"""k>0时的积分虚部被积函数"""
function imag_integrand_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    coeff_x, denominator_const = compute_coefficients(λ, k, m, m_prime, E)
    p = internal_momentum(E,m)
    
    # 计算奇点位置
    x_pole = -denominator_const / coeff_x
    
    # 检查奇点是否在积分区间[-1,1]内
    if -1.0 ≤ x_pole ≤ 1.0
        # 虚部贡献
        value = 2.0 * p * imag_integral_tool(coeff_x, denominator_const) *
                correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
        return value
    else
        # 奇点不在积分区间内，无虚部贡献
        return 0.0
    end
end
```

## 虚部计算的物理解释

### k>0 情况

对于k>0的情况，虚部来源于角度积分中的极点：

1. **极点位置**：当分母 $Ax + B = 0$ 时，即 $x = -B/A$ 时出现极点
2. **虚部贡献条件**：只有当极点位于物理积分区间 $[-1, 1]$ 内时才贡献虚部
3. **虚部大小**：由 `imag_integral_tool` 函数计算，对应Sokhotski-Plemelj公式中的 $\pi$ 乘以留数

### k=0 情况

对于k=0的情况：
- 分母简化为常数 $B$，没有$x$依赖
- 因此没有极点，虚部为零

## 完整的各向异性修正项计算

结合实部和虚部，完整的各向异性修正项为：

```julia
"""计算完整的各向异性修正项（实部 + i×虚部）"""
function compute_full_B0_aniso_correction(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64, E_min::Float64, E_max::Float64)
    
    # 计算实部
    if k > EPS_K
        real_part, real_err = quadgk(E -> real_integrand_k_positive(sign_, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar), 
                                    E_min, E_max, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL)
        imag_part, imag_err = quadgk(E -> imag_integrand_k_positive(sign_, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar), 
                                    E_min, E_max, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL)
    else
        real_part, real_err = quadgk(E -> real_integrand_k_zero(sign_, λ, m, m_prime, E, ξ, T, μ, Φ, Φbar), 
                                    E_min, E_max, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL)
        imag_part = 0.0
        imag_err = 0.0
    end
    
    return complex(real_part, imag_part), (real_err, imag_err)
end
```

## 注意事项

1. **数值稳定性**：当极点非常接近积分边界时，可能需要特殊处理
2. **收敛性检查**：虚部积分应该在奇点附近有较好的收敛性，因为它是通过解析延拓得到的
3. **物理意义**：虚部通常对应物理过程中的衰减或共振现象

这样的实现确保了各向异性修正项在数学上的完整性，同时考虑了物理过程中的虚部贡献。