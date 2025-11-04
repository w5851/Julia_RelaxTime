# [有限温度有限密度场论中的单圈积分A（one-line integral）] - 实现文档

## 1. 公式标识
- **程序实现**: `A`
- **对应文件**: `../../src/relaxtime/OneLoopIntegrals.jl`

## 2. 数学表达式
### 2.1 原始公式
```math
A(m, \mu, \beta, \Lambda) = \frac{16\pi^2}{\beta} \sum_n \exp(i\omega_n \eta) \int \frac{d^3p}{(2\pi)^3} \frac{1}{( (i\omega_n + \mu)^2 - E^2 )}
```

### 2.2 数学处理(松原频率求和)后的形式
```math
A(m, \mu, \beta, \Lambda) = 16\pi^2 \int \frac{d^3p}{(2\pi)^3} \frac{1}{2E} \left(f(E - \mu) + f(E + \mu) -1
\right)
```
```math
A(m, \mu, \beta, \Lambda) = 4\int_0^{\Lambda} p^2 dp \frac{1}{E} \left(f(E - \mu) + f(E + \mu) -1
\right)
```
## 3. 额外说明
- PNJL模型下，参数需要加上Φ和Φbar；f(E-μ)和f(E+μ)分别表示夸克和反夸克的分布函数，具体形式见**PNJL_夸克有效分布函数_补充详细说明.md**。
- 计算时对动量积分即可，传入预生成的高斯-勒让德节点计算，**对常数1使用三维动量积分的数值截断Λ，对分布函数部分使用无限积分**。

### 3.1. 对常数项的数学处理
# 不定积分 ∫x²/√(x²+m²) dx 的初等函数表示

## 问题描述
考虑积分：
\[
\int \frac{x^2}{\sqrt{x^2 + m^2}}  dx
\]
其中 \(m\) 是参数。

## 解答
该积分存在初等函数表示，结果为：
\[
\int \frac{x^2}{\sqrt{x^2 + m^2}}  dx = \frac{x}{2} \sqrt{x^2 + m^2} - \frac{m^2}{2} \ln \left| x + \sqrt{x^2 + m^2} \right| + C
\]
其中 \(C\) 为积分常数。

## 推导过程

该结果可以通过分部积分或三角代换等方法推导得出。

### 方法一：代数变形与分部积分

设 \(I = \int \frac{x^2}{\sqrt{x^2 + m^2}}  dx\)，通过变形 \(x^2 = (x^2 + m^2) - m^2\)，可得：

\[
I = \int \frac{(x^2 + m^2) - m^2}{\sqrt{x^2 + m^2}}  dx = \int \sqrt{x^2 + m^2}  dx - m^2 \int \frac{1}{\sqrt{x^2 + m^2}}  dx
\]

其中：
- \(\int \frac{1}{\sqrt{x^2 + m^2}}  dx = \ln \left| x + \sqrt{x^2 + m^2} \right| + C\)
- \(\int \sqrt{x^2 + m^2}  dx\) 可通过分部积分求解

最终得到上述结果。

### 方法二：三角代换

令 \(x = m \tan \theta\)，则 \(dx = m \sec^2 \theta  d\theta\)，且：
\[
\sqrt{x^2 + m^2} = \sqrt{m^2 \tan^2 \theta + m^2} = m \sec \theta
\]

代入原积分：
\[
\int \frac{x^2}{\sqrt{x^2 + m^2}}  dx = \int \frac{m^2 \tan^2 \theta}{m \sec \theta} \cdot m \sec^2 \theta  d\theta = m^2 \int \tan^2 \theta \sec \theta  d\theta
\]

利用三角恒等式 \(\tan^2 \theta = \sec^2 \theta - 1\)，继续化简可得相同结果。

## 结论
因此，该积分可以表示为初等函数。
# 定积分 ∫₀^{x_max} x²/√(x²+m²) dx 的计算结果

## 定积分表达式

考虑定积分：
\[
\int_{0}^{x_{\text{max}}} \frac{x^2}{\sqrt{x^2 + m^2}}  dx
\]

## 计算结果

根据不定积分的结果，计算该定积分得到：
\[
\int_{0}^{x_{\text{max}}} \frac{x^2}{\sqrt{x^2 + m^2}}  dx = \frac{x_{\text{max}}}{2} \sqrt{x_{\text{max}}^2 + m^2} - \frac{m^2}{2} \ln \left( \frac{x_{\text{max}} + \sqrt{x_{\text{max}}^2 + m^2}}{m} \right)
\]

### 4. A的单位和量纲
- 由于积分变量 \(x\) 的单位为动量（fm⁻¹），被积函数的分母包含能量项（√(x² + m²)），因此整个积分的结果 \(A\) 的单位为 fm⁻²。
- 量纲为[能量]²,即对应单位 fm⁻²。