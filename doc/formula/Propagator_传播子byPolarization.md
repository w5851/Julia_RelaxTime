# 介子传播子计算 - 实现文档

## 1. 公式标识
- **程序实现**: `meson_propagator`
- **对应文件**: `src/relaxtime/MesonPropagator.jl`

## 2. 物理意义
- **物理背景**: 在PNJL模型中，介子被视为夸克-反夸克相互作用的中间玻色子，通过随机相位近似(RPA)构建介子传播子，用于计算夸克弹性散射过程中的散射截面
- **在计算链中的作用**: 作为夸克两体散射计算的核心输入，直接影响驰豫时间和输运系数的计算结果
- **相关物理量**: 介子质量、衰变宽度、极化函数、散射矩阵元

## 3. 数学表达式

### 3.1 一般介子传播子(π, K, σ_K, σ_π)(σ_K和σ_π是标量介子S，而π和K是赝标量介子P)
```math
\mathcal{D}_{\pi} = \frac{2K_3^{\pm}}{1 - 4K_3^{\pm} \prod_{u\bar u}^{P(S)}(k_0, k)}
```
```math
\mathcal{D}_K = \frac{2K_4^{\pm}}{1 - 4K_4^{\pm} \prod_{u\bar u}^{P(S)}(k_0, k)}
```

### 3.2 混合介子传播子（η/η', σ/σ'）(η/η'是赝标量介子P，σ/σ'是标量介子S)
```math
\mathcal{D} = 2\frac{\det K}{M_{00}M_{88}-M_{08}^2}(M_{00}\bar\psi\lambda_0\psi\cdot\bar\psi'\lambda_0\psi' + M_{08}\bar\psi\lambda_0\psi\cdot\bar\psi'\lambda_8\psi' + M_{80}\bar\psi\lambda_8\psi\cdot\bar\psi'\lambda_0\psi' + M_{88}\bar\psi\lambda_8\psi\cdot\bar\psi'\lambda_8\psi')
```
- 其中**λ是盖尔曼矩阵**,**ψ与其变体是夸克味空间波函数**,详情见**3.3小节**

其中：
```math
M_{00} = K_0^+ - \frac{4}{3} \det K (\Pi_{u\bar u}^{P(S)} + 2 \Pi_{s\bar s}^{P(S)})
```
```math
M_{08} = K_0^+ - \frac{4}{3} \sqrt{2} \det K (\Pi_{u\bar u}^{P(S)} - \Pi_{s\bar s}^{P(S)})
```
```math
M_{88} = K_0^+ - \frac{4}{3} \det K (2 \Pi_{u\bar u}^{P(S)} + \Pi_{s\bar s}^{P(S)})
```
```math
\det K = K_0^+ K_8^+ - K_{08}^2
```
- 这个传播子可以写成**矩阵乘法**的形式,推导如下:

**定义以下矩阵和向量**：

#### 3.2.1 耦合矩阵
$$
M = \begin{pmatrix}
M_{00} & M_{08} \\
M_{80} & M_{88}
\end{pmatrix}
$$

#### 3.2.2 流算符向量
$$
J = \begin{pmatrix}
\bar{\psi}\lambda_0\psi \\
\bar{\psi}\lambda_8\psi
\end{pmatrix}, \quad 
J' = \begin{pmatrix}
\bar{\psi}'\lambda_0\psi' \\
\bar{\psi}'\lambda_8\psi'
\end{pmatrix}
$$

#### 3.2.3 矩阵形式的传播子
$$
\mathcal{D} = 2\frac{\det K}{\det M} J^T M J'
$$

其中 $\det M = M_{00}M_{88} - M_{08}M_{80}$。

#### 3.2.4 详细展开验证

将矩阵乘法展开：
$$
J^T M J' = \begin{pmatrix}
\bar{\psi}\lambda_0\psi & \bar{\psi}\lambda_8\psi
\end{pmatrix}
\begin{pmatrix}
M_{00} & M_{08} \\
M_{80} & M_{88}
\end{pmatrix}
\begin{pmatrix}
\bar{\psi}'\lambda_0\psi' \\
\bar{\psi}'\lambda_8\psi'
\end{pmatrix}
$$

计算中间步骤：
$$
\begin{pmatrix}
\bar{\psi}\lambda_0\psi & \bar{\psi}\lambda_8\psi
\end{pmatrix}
\begin{pmatrix}
M_{00} & M_{08} \\
M_{80} & M_{88}
\end{pmatrix}
= \begin{pmatrix}
M_{00}\bar{\psi}\lambda_0\psi + M_{80}\bar{\psi}\lambda_8\psi & 
M_{08}\bar{\psi}\lambda_0\psi + M_{88}\bar{\psi}\lambda_8\psi
\end{pmatrix}
$$

再与 $J'$ 相乘：
$$
[M_{00}\bar{\psi}\lambda_0\psi + M_{80}\bar{\psi}\lambda_8\psi] \bar{\psi}'\lambda_0\psi' + 
[M_{08}\bar{\psi}\lambda_0\psi + M_{88}\bar{\psi}\lambda_8\psi] \bar{\psi}'\lambda_8\psi'
$$

整理后得到：
$$
M_{00}\bar{\psi}\lambda_0\psi\cdot\bar{\psi}'\lambda_0\psi' + 
M_{80}\bar{\psi}\lambda_8\psi\cdot\bar{\psi}'\lambda_0\psi' + 
M_{08}\bar{\psi}\lambda_0\psi\cdot\bar{\psi}'\lambda_8\psi' + 
M_{88}\bar{\psi}\lambda_8\psi\cdot\bar{\psi}'\lambda_8\psi'
$$

这与原始表达式完全一致。
- **(参考文献Hadronization in the SU(3) Nambu–Jona-Lasinio model和Elastic scattering and transport coefficients for a quark plasma in SUf(3) at finite temperatures)**

### 3.3 盖尔曼矩阵与夸克味空间波函数
#### 3.3.1 **Gell-Mann 矩阵 λₐ（3×3）**

在 SU(3) 味空间中，Gell-Mann 矩阵是 3×3 的厄米矩阵，满足：

\[
\text{Tr}(\lambda_a \lambda_b) = 2\delta_{ab}
\]

常用的 λₐ 包括：

- **λ₀（味单态，归一化）**：
\[
\lambda_0 = \sqrt{\frac{2}{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
\]

- **λ₃（同位旋第三分量）**：
\[
\lambda_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
\]

- **λ₈（超荷）**：
\[
\lambda_8 = \frac{1}{\sqrt{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -2 \end{pmatrix}
\]

- **λ₄ 和 λ₅（K⁺ 介子相关）**：
\[
\lambda_4 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \quad
\lambda_5 = \begin{pmatrix} 0 & 0 & -i \\ 0 & 0 & 0 \\ i & 0 & 0 \end{pmatrix}
\]

- **λ₆ 和 λ₇（K⁰ 介子相关）**：
\[
\lambda_6 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}, \quad
\lambda_7 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -i \\ 0 & i & 0 \end{pmatrix}
\]

---

#### 3.3.2 **夸克与反夸克的味空间波函数（3×1 或 1×3）**

在味空间中，每个夸克用一个 **3 维列向量** 表示，反夸克用 **3 维行向量** 表示：

##### 3.3.2.1 夸克（列向量）：
- **u 夸克**：
\[
\psi_u = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}
\]

- **d 夸克**：
\[
\psi_d = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}
\]

- **s 夸克**：
\[
\psi_s = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
\]

##### 3.3.2.2 反夸克（行向量）：
- **ū 反夸克**：
\[
\bar{\psi}_u = \begin{pmatrix} 1 & 0 & 0 \end{pmatrix}
\]

- **d̄ 反夸克**：
\[
\bar{\psi}_d = \begin{pmatrix} 0 & 1 & 0 \end{pmatrix}
\]

- **s̄ 反夸克**：
\[
\bar{\psi}_s = \begin{pmatrix} 0 & 0 & 1 \end{pmatrix}
\]

---

#### 3.3.3 **散射过程中的 ψ, ψ̄, ψ′, ψ̄′**

在夸克-反夸克散射过程中：

- **ψ**：入射夸克的味波函数（列向量）
- **ψ̄**：入射反夸克的味波函数（行向量）
- **ψ′**：出射夸克的味波函数（列向量）
- **ψ̄′**：出射反夸克的味波函数（行向量）

##### 示例：\( u\bar{s} \to u\bar{s} \) 过程

- **入射**：u 夸克 + s̄ 反夸克
  \[
  \psi = \psi_u = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \quad
  \bar{\psi} = \bar{\psi}_s = \begin{pmatrix} 0 & 0 & 1 \end{pmatrix}
  \]

- **出射**：u 夸克 + s̄ 反夸克
  \[
  \psi' = \psi_u = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \quad
  \bar{\psi}' = \bar{\psi}_s = \begin{pmatrix} 0 & 0 & 1 \end{pmatrix}
  \]

##### 总结

- **λₐ** 是 3×3 的 Gell-Mann 矩阵，描述夸克在味空间的变换。
- **ψ** 是夸克的味波函数（列向量），**ψ̄** 是反夸克的味波函数（行向量）。
- **ψ′** 和 **ψ̄′** 是出射态的味波函数。

### 3.4 极化函数
```math
\Pi_{f}^{P,S}(p_0,k) = -\frac{N_c}{8\pi^2}\left\{A(m_f,\mu_f,T) + A(m_{f'},\mu_{f'},T) + [(m_f \mp m_{f'})^2 - (p_0 + \mu_f - \mu_{f'})^2 + p^2] \times B_0(|\vec{p}|, m_f,\mu_f,m_{f'},\mu_{f'},p_0,T)\right\}
```

## 4. 参数说明表

| 参数 | 符号 | 类型 | 单位 | 物理意义 | 取值范围 |
|------|------|------|------|----------|----------|
| meson_type | - | 输入 | - | 介子类型标识 | [:pi, :sigma_pi, :K, :sigma_K, :eta, :eta_prime, :sigma, :sigma_prime] |
| k0 | $k_0$ | 输入 | fm⁻¹ | 能量分量 | [0, ∞) |
| k | $\vec{k}$ | 输入 | fm⁻¹ | 动量 | [0, ∞) |
| T | $T$ | 输入 | fm⁻¹ | 温度 | [0, 2.03] (对应0-400 MeV) |
| μ | $\mu_q$ | 输入 | fm⁻¹ | 夸克化学势 | [0, 5.07] (对应0-1000 MeV) |
| m_q | $m_q$ | 输入 | fm⁻¹ | 夸克质量 | 模型相关 |
| K | $K_\alpha^+$ | 输入 | fm² | 耦合常数 | 模型参数 |

## 5. 输入参数详细说明

### 5.1 必需参数
- **meson_type**: 指定计算的介子类型，影响耦合常数和味因子的选择
- **k0, k**: 四动量的能量和动量分量，决定传播子的运动学行为
- **T, μ**: 热力学参数，影响夸克分布函数和有效质量
- **m_q**: 夸克质量数组，包含u,d,s夸克在当前温度化学势下的有效质量

### 5.2 可选参数
- **Γ**: 介子衰变宽度，默认通过自洽求解获得，可手动指定用于测试
- **cutoff_Λ**: 动量截断，默认使用模型参数中的截断值

## 6. 输出结果说明

### 6.1 主输出
- **propagator**: 复数类型的介子传播子，包含实部和虚部，单位fm²

### 6.2 辅助输出
- **mass**: 介子物理质量，通过求解$Re[\mathcal{D}^{-1}(m,\vec{0})]=0$获得，单位fm⁻¹
- **width**: 介子衰变宽度，反映传播子虚部的大小，单位fm⁻¹
- **Mott_flag**: Mott相变标志，当$m_{meson} = m_{q1} + m_{q2}$时触发

## 7. 依赖关系

### 7.1 输入依赖
- **polarization_function**: 计算夸克-反夸克极化函数，包含主值积分和虚部计算
- **quark_mass**: 提供当前热力学条件下的夸克有效质量
- **coupling_constant**: 根据介子类型返回相应的耦合常数
- **distribution_function**: 费米-狄拉克分布函数，用于极化函数计算

### 7.2 输出用途
- **scattering_amplitude**: 被散射矩阵元计算函数调用，用于计算微分散射截面
- **relaxation_time**: 间接影响驰豫时间的计算
- **transport_coefficients**: 最终影响剪切粘滞系数等输运系数的计算结果

### 7.3 物理常数
- **N_c**: 色数，$N_c = 3$
- **hbar_c**: 自然单位制转换常数，hbar_c = 197.327 MeV·fm

**单位转换说明**：
- 能量单位转换：1 MeV = 1/197.327 fm⁻¹
- 耦合常数单位转换：1 MeV⁻² = (197.327)² fm²

## 8. 程序实现细节计划
- 计算一般介子传播子时，将极化函数的值作为参数传入而非在函数内计算，以避免重复计算极化函数
- 计算混合介子传播子时，直接使用矩阵乘法计算，避免手动展开；另外，由于耦合矩阵M只与极化函数有关，因此将耦合矩阵M也作为参数传入，在批量调用函数前预计算耦合矩阵即可
- 混合介子传播子的矩阵乘法计算与散射过程中的夸克类型有关，计划将u,d,s三种夸克的味波函数预定义为常量，直接在散射过程中调用，避免每次计算时重复创建向量；另外，这样需要将散射过程涉及到的夸克类型作为参数传入函数，以选择正确的味波函数进行计算