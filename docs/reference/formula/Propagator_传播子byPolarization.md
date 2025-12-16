# 介子传播子计算 - 实现文档

## 1. 公式标识
- **程序实现**: `meson_propagator`
- **对应文件**: `src/relaxtime/MesonPropagator.jl`

## 2. 物理意义
- **物理背景**: 在PNJL模型中，介子被视为夸克-反夸克相互作用的中间玻色子，通过随机相位近似(RPA)构建介子传播子，用于计算夸克弹性散射过程中的散射截面
- **在计算链中的作用**: 作为夸克两体散射计算的核心输入，直接影响驰豫时间和输运系数的计算结果
- **相关物理量**: 介子质量、衰变宽度、极化函数、散射矩阵元

## 3. 数学表达式

### 3.1 一般介子传播子(π, K, σ_π, σ_K)

**介子类型与通道说明**：
- **赝标量介子P**（使用K⁺系数）：π介子（uū/d̄d组合）、K介子（us̄/ds̄组合）
- **标量介子S**（使用K⁻系数）：σ_π介子（uū/d̄d组合）、σ_K介子（us̄/ds̄组合）

**传播子公式**（使用复数除法）：

**π介子传播子**（赝标量P通道）：
```math
\mathcal{D}_{\pi} = \frac{2K_3^{+}}{1 - 4K_3^{+} \Pi_{u\bar u}^{P}(k_0, k)}
```
其中 $\Pi_{u\bar u}^{P} = \Pi_{u\bar u}^{P,\text{real}} + i\Pi_{u\bar u}^{P,\text{imag}}$ 是复数。

**K介子传播子**（赝标量P通道）：
```math
\mathcal{D}_{K} = \frac{2K_4^{+}}{1 - 4K_4^{+} \Pi_{u\bar s}^{P}(k_0, k)}
```
注意：K介子使用 $\Pi_{u\bar s}^{P}$ 而非 $\Pi_{u\bar u}^{P}$。

**σ_π标量介子传播子**（标量S通道）：
```math
\mathcal{D}_{\sigma_\pi} = \frac{2K_3^{-}}{1 - 4K_3^{-} \Pi_{u\bar u}^{S}(k_0, k)}
```

**σ_K标量介子传播子**（标量S通道）：
```math
\mathcal{D}_{\sigma_K} = \frac{2K_4^{-}}{1 - 4K_4^{-} \Pi_{u\bar s}^{S}(k_0, k)}
```

**返回值格式**：
- 传播子 $\mathcal{D}$ 是复数类型 `ComplexF64`
- 可通过 `real(D)` 和 `imag(D)` 提取实部和虚部
- 单位：fm²

### 3.2 混合介子传播子（η/η', σ/σ'）

**介子类型与通道说明**：
- **赝标量混合介子P**（使用K⁺系数）：η/η'混合（uū、s̄s通过λ₀、λ₈耦合）
- **标量混合介子S**（使用K⁻系数）：σ/σ'混合

**传播子公式**（使用复数矩阵运算）：
```math
\mathcal{D} = 2\frac{\det K}{\det M} J^T M J'
```

展开形式：
```math
\mathcal{D} = 2\frac{\det K}{M_{00}M_{88}-M_{08}^2}(M_{00}\bar\psi\lambda_0\psi\cdot\bar\psi'\lambda_0\psi' + M_{08}\bar\psi\lambda_0\psi\cdot\bar\psi'\lambda_8\psi' + M_{08}\bar\psi\lambda_8\psi\cdot\bar\psi'\lambda_0\psi' + M_{88}\bar\psi\lambda_8\psi\cdot\bar\psi'\lambda_8\psi')
```
- 其中**λ是盖尔曼矩阵**，**ψ与其变体是夸克味空间波函数**，详情见**3.3小节**
- **M矩阵是对称的**：$M_{08} = M_{80}$
- **所有量均为复数**：$M_{ij}$、$\det K$、$\det M$、$\mathcal{D}$ 都是 `ComplexF64` 类型

**耦合矩阵 $M$ 的定义**（对赝标量P使用K⁺，标量S使用K⁻）：

对于**赝标量通道P**（η/η'）：
```math
M_{00}^P = K_0^{+} - \frac{4}{3} \det K^{+} (\Pi_{u\bar u}^{P} + 2 \Pi_{s\bar s}^{P})
```
```math
M_{08}^P = M_{80}^P = K_{08}^{+} + \frac{4\sqrt{2}}{3} \det K^{+} (\Pi_{u\bar u}^{P} - \Pi_{s\bar s}^{P})
```
```math
M_{88}^P = K_8^{+} - \frac{4}{3} \det K^{+} (2 \Pi_{u\bar u}^{P} + \Pi_{s\bar s}^{P})
```

对于**标量通道S**（σ/σ'）：
```math
M_{00}^S = K_0^{-} - \frac{4}{3} \det K^{-} (\Pi_{u\bar u}^{S} + 2 \Pi_{s\bar s}^{S})
```
```math
M_{08}^S = M_{80}^S = K_{08}^{-} + \frac{4\sqrt{2}}{3} \det K^{-} (\Pi_{u\bar u}^{S} - \Pi_{s\bar s}^{S})
```
```math
M_{88}^S = K_8^{-} - \frac{4}{3} \det K^{-} (2 \Pi_{u\bar u}^{S} + \Pi_{s\bar s}^{S})
```

其中：
```math
\det K^{\pm} = K_0^{\pm} K_8^{\pm} - (K_{08}^{\pm})^2
```

**注意**：$M_{08}$ 公式中的系数是 $\frac{4\sqrt{2}}{3}$（即 $\frac{4}{3}\sqrt{2}$），不是 $\frac{4}{3\sqrt{2}}$。
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

**本模块使用的Gell-Mann矩阵**（需在Constants_PNJL.jl中定义为常量）：

- **λ₀（味单态，归一化）**：
\[
\lambda_0 = \sqrt{\frac{2}{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
\]

- **λ₈（超荷）**：
\[
\lambda_8 = \frac{1}{\sqrt{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -2 \end{pmatrix}
\]

**其他Gell-Mann矩阵**（供参考，本模块暂不需要）：

- **λ₃（同位旋第三分量）**：
\[
\lambda_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
\]

- **λ₄ 和 λ₅（K⁺ 介子相关）**：
\[
\lambda_4 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}, \quad
\lambda_5 = \begin{pmatrix} 0 & 0 & -i \\ 0 & 0 & 0 \\ i & 0 & 0 \end{pmatrix}
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

##### 流算符 J 的计算示例

对于 u 夸克和 s̄ 反夸克的混合介子传播子：
\[
J = \begin{pmatrix}
\bar{\psi}_s \lambda_0 \psi_u \\
\bar{\psi}_s \lambda_8 \psi_u
\end{pmatrix}
\]

具体计算（以 u 夸克为例）：
- $\bar{\psi}_u \lambda_0 \psi_u = [1\;0\;0] \cdot \sqrt{\frac{2}{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \cdot \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} = \sqrt{\frac{2}{3}}$
- $\bar{\psi}_u \lambda_8 \psi_u = [1\;0\;0] \cdot \frac{1}{\sqrt{3}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -2 \end{pmatrix} \cdot \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} = \frac{1}{\sqrt{3}}$

对于 s 夸克：
- $\bar{\psi}_s \lambda_0 \psi_s = \sqrt{\frac{2}{3}}$
- $\bar{\psi}_s \lambda_8 \psi_s = \frac{1}{\sqrt{3}} \times (-2) = -\frac{2}{\sqrt{3}}$

**性能注记**：这些矩阵乘法可以预计算为常量，但由于计算开销相对极化函数和传播子可忽略，也可直接在函数内部计算。

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
- **流算符 J** 是通过矩阵乘法 $\bar{\psi}\lambda\psi$ 计算得到的2×1列向量。

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

## 8. 程序实现细节说明

### 8.1 返回值格式与数据类型
- **所有传播子函数返回 `ComplexF64` 类型**
- 可通过 `real(D)` 和 `imag(D)` 提取实部和虚部
- 单位：fm²

### 8.2 一般介子传播子实现
- **函数签名**：`meson_propagator_simple(meson_type::Symbol, channel::Symbol, K_coeffs::NamedTuple, Π::ComplexF64) -> ComplexF64`
- **参数说明**：
  - `meson_type`: 介子类型(`:pi`, `:K`, `:sigma_pi`, `:sigma_K`)
  - `channel`: 通道类型(`:P`赝标量或`:S`标量)
  - `K_coeffs`: 通过`EffectiveCouplings.calculate_effective_couplings`预先计算的K系数
  - `Π`: 预计算的极化函数(ComplexF64)
- **K系数自动选择**：函数内部根据`meson_type`和`channel`自动选择正确的K系数，实现逻辑如下：
  ```julia
  if meson_type == :pi
      K_coeff = (channel == :P) ? K_coeffs.K123_minus : K_coeffs.K123_plus
  elseif meson_type == :K
      K_coeff = (channel == :P) ? K_coeffs.K4567_minus : K_coeffs.K4567_plus
  elseif meson_type == :sigma_pi
      K_coeff = (channel == :S) ? K_coeffs.K123_plus : K_coeffs.K123_minus
  elseif meson_type == :sigma_K
      K_coeff = (channel == :S) ? K_coeffs.K4567_plus : K_coeffs.K4567_minus
  end
  ```
- **性能优化**：将极化函数的复数值 `Π = Π_real + im*Π_imag` 作为参数传入，避免重复计算；K系数选择的条件判断开销可忽略不计
- **介子类型与K系数映射**（**关键修正**）：
  | 介子 | 通道 | 夸克组合 | 使用的Π | K系数 |
  |------|------|----------|---------|-------|
  | π | P | uū/d̄d | Π_{uu}^P | K₃⁺ |
  | K | P | us̄/ds̄ | Π_{us}^P | K₄⁺ |
  | σ_π | S | uū/d̄d | Π_{uu}^S | K₃⁻ |
  | σ_K | S | us̄/ds̄ | Π_{us}^S | K₄⁻ |

### 8.3 混合介子传播子实现
- **函数签名**：`meson_propagator_mixed(det_K::Float64, M_matrix::Matrix{ComplexF64}, q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol, channel::Symbol) -> ComplexF64`
- **参数说明**：
  - `q1, q2, q3, q4`: 散射过程q1+q2→q3+q4中的粒子类型(`:u`, `:d`, `:s`, `:ubar`, `:dbar`, `:sbar`)
  - `channel`: 散射道类型(`:t`, `:s`, `:u`)
- **散射过程约束**：只有q2和q4可能是反夸克，且若有反夸克则q2和q4同时为反夸克
- **场算符映射**：函数内部根据channel自动将q1,q2,q3,q4映射到ψ,ψ̄,ψ',ψ̄'场算符，调用者无需处理复杂的场算符对应关系
- **矩阵运算**：使用Julia内置的复数矩阵运算 `J' * M * J'`
- **耦合矩阵M**：2×2复数矩阵，作为参数传入（在批量调用前预计算）
- **味波函数常量**：在 `Constants_PNJL.jl` 中预定义（使用ASCII bar命名）：
  ```julia
  const ψ_u = [1.0, 0.0, 0.0]  # 列向量
  const ψ_d = [0.0, 1.0, 0.0]
  const ψ_s = [0.0, 0.0, 1.0]
  const ψbar_u = [1.0 0.0 0.0]  # 行向量（1×3矩阵）
  const ψbar_d = [0.0 1.0 0.0]
  const ψbar_s = [0.0 0.0 1.0]
  ```
- **Gell-Mann矩阵常量**：在 `Constants_PNJL.jl` 中定义：
  ```julia
  const λ₀ = sqrt(2/3) * [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
  const λ₈ = (1/sqrt(3)) * [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -2.0]
  ```
- **辅助函数**：
  - `extract_flavor(q::Symbol)`: 从Symbol中提取味类型(去除bar标记)，如`:ubar→:u`
  - `get_quark_wavefunction(flavor::Symbol, is_bar::Bool)`: 根据味类型和是否为bar返回对应的列向量或行向量波函数

### 8.4 辅助函数
- **`calculate_coupling_matrix`**:根据 `(Π_uu::ComplexF64, Π_ss::ComplexF64)`、预计算的K系数 `K_coeffs::NamedTuple` 和通道类型 `channel::Symbol` 计算复数矩阵M(2×2)
- **函数签名**:`calculate_coupling_matrix(Π_uu::ComplexF64, Π_ss::ComplexF64, K_coeffs::NamedTuple, channel::Symbol) -> Matrix{ComplexF64}`
- **采用方案A**:K系数通过`EffectiveCouplings.calculate_effective_couplings(G, K, G_u, G_s)`预先计算并作为参数传入,函数内部根据channel自动选择对应的K系数(`:P`通道用K⁺系数,`:S`通道用K⁻系数)
- **性能优化**:批量计算多个传播子时,可以只计算一次K系数并复用,避免重复调用`calculate_effective_couplings`
- **M矩阵对称性**：利用 `M₀₈ = M₈₀` 减少计算量
- **det_K计算**：`det_K`为Float64类型的实数，已在EffectiveCouplings.jl中定义
- **det_M计算**：直接使用Julia标准库的`det()`函数计算M矩阵行列式，无需额外定义函数

### 8.5 关键公式修正总结
1. **K介子极化函数**：使用 `Π_{us}` 而非 `Π_{uu}` ✓
2. **通道与K系数映射**（**已修正**）：
   - 赝标量通道P（π、K、η/η'）→ 使用 K⁺ 系数
   - 标量通道S（σ_π、σ_K、σ/σ'）→ 使用 K⁻ 系数
3. **M₀₈系数**：$\frac{4\sqrt{2}}{3}$（即 $\frac{4}{3}\sqrt{2}$）✓
4. **数据类型**：
   - 所有传播子和M矩阵使用 `ComplexF64` 类型 ✓
   - 极化函数Π作为 `ComplexF64` 传入 ✓
   - K系数和det_K是 `Float64` 类型的实数 ✓
5. **M矩阵对称性**：$M_{08} = M_{80}$ ✓