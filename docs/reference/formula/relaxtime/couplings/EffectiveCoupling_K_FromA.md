# 有效耦合系数 - 实现文档

## 1. 公式标识
- **程序实现**: `K_effective_couplings`
- **对应文件**: `../../src/relaxtime/EffectiveCouplings.jl`

## 2. 物理意义
- **物理背景**: 在3味PNJL模型中，有效耦合系数描述了夸克-反夸克相互作用通过介子传播的有效强度，是随机相位近似下介子传播子的关键参数
- **在计算链中的作用**: 作为介子传播子和散射截面的输入，直接影响驰豫时间和输运系数的计算
- **相关物理量**: 原始耦合常数G、K，夸克凝聚 \(\phi_u,\phi_d,\phi_s\)，温度T，化学势μ

## 3. 数学表达式

### 3.1 原始公式

为与当前物理凝聚口径区分，下面先记历史 helper
`H_u=-\phi_u`、`H_s=-\phi_s`；它们对应代码参数名 `G_u`、`G_s`。

```math
\begin{aligned}
K_0^{\pm} &= G \mp \frac{1}{3} K(2H_u + H_s) \\
K_1^{\pm} = K_2^{\pm} = K_3^{\pm} &= G \pm \frac{1}{2} KH_s \\
K_4^{\pm} = K_5^{\pm} = K_6^{\pm} = K_7^{\pm} &= G \pm \frac{1}{2} KH_u \\
K_8^{\pm} &= G \pm \frac{1}{6} K(4H_u - H_s) \\
K_{08}^{\pm} &= \pm \frac{1}{6} \sqrt{2} K(H_u - H_s)
\end{aligned}
```

其中，物理上应优先从夸克凝聚
```math
\begin{aligned}
\phi_f &= \frac{N_c}{4\pi^2} \; m_f\; A_f(T, \mu), \\
A_f(T, \mu) &= -16\pi^2 \int_0^\Lambda \frac{d^3 p}{(2\pi)^3} \frac{1}{2E}\left[ 1 - n_f^+ - n_f^- \right].
\end{aligned}
```
出发，其中
```math
\phi_f \equiv \langle \bar q_f q_f \rangle.
```

等价地，把角积分做掉后，
```math
A_f(T,\mu)= -4\int_0^\Lambda dp\,\frac{p^2}{E}
\left[1-n_q(E;\mu,\Phi,\bar\Phi)-n_{\bar q}(E;\mu,\Phi,\bar\Phi)\right].
```

若改写成当前项目代码直接实现的积分核，则同一个 `A_f` 也可写成
```math
A_f(T,\mu)= 4\int_0^\Lambda dp\,\frac{p^2}{E}
\left[n_q(E;\mu,\Phi,\bar\Phi)+n_{\bar q}(E;\mu,\Phi,\bar\Phi)-1\right].
```

在项目的自然单位约定下，`p,m_f,T,mu_f` 的单位为 `fm^-1`，因此
```math
[A_f]=\mathrm{fm}^{-2},\qquad
[\phi_f]=\mathrm{fm}^{-3},\qquad
[H_f]=\mathrm{fm}^{-3},\qquad
[K H_f]=\mathrm{fm}^{2}.
```

这两式之所以等价，是因为第二式把第一式中整体的负号吸收到括号里，改写成了 `n_q+n_{\bar q}-1`；这里的物理凝聚主口径始终是 `\phi_f`，而不是 `-\phi_f`。

- 需要注意：当前实现中的辅助函数 `calculate_G_from_A(A_f, m_f)` 返回的是
```math
G^f_{\text{helper}} = -\phi_f,
```
这是历史实现中为适配后续有效耦合公式而保留的中间变量口径；在文档解释层，应优先把 `\phi_f` 视为基本物理量。

### 3.2 程序实现形式
```julia
function calculate_effective_couplings(G, K, G_u, G_s)
    K0_plus = G - (1/3)*K*(2*G_u + G_s)
    K0_minus = G + (1/3)*K*(2*G_u + G_s)
    
    K123_plus = G + (1/2)*K*G_s
    K123_minus = G - (1/2)*K*G_s
    
    K4567_plus = G + (1/2)*K*G_u
    K4567_minus = G - (1/2)*K*G_u
    
    K8_plus = G + (1/6)*K*(4*G_u - G_s)
    K8_minus = G - (1/6)*K*(4*G_u - G_s)
    
    K08_plus = (1/6)*√2*K*(G_u - G_s)
    K08_minus = -(1/6)*√2*K*(G_u - G_s)
    
    return (K0_plus, K0_minus, K123_plus, K123_minus, 
            K4567_plus, K4567_minus, K8_plus, K8_minus,
            K08_plus, K08_minus)
end
```

## 4. 参数说明表

| 参数 | 符号 | 类型 | 单位 | 物理意义 | 取值范围 |
|------|------|------|------|----------|----------|
| G | G | 输入 | fm² | 四夸克相互作用耦合常数 | 正实数 |
| K | K | 输入 | fm⁵ | 't Hooft相互作用耦合常数 | 实数 |
| G_u | \(G^u_{\text{helper}}=-\phi_u\) | 输入 | fm⁻³ | 当前实现中进入耦合公式的历史中间变量 | 实数 |
| G_s | \(G^s_{\text{helper}}=-\phi_s\) | 输入 | fm⁻³ | 当前实现中进入耦合公式的历史中间变量 | 实数 |
| Nc | N_c | 输入 | 无量纲 | 色数 | 3(固定) |
| K0± | K₀± | 输出 | fm² | 单态道有效耦合系数 | 实数 |
| K123± | K₁₂₃± | 输出 | fm² | π介子道有效耦合系数 | 实数 |
| K4567± | K₄₅₆₇± | 输出 | fm² | K介子道有效耦合系数 | 实数 |
| K8± | K₈± | 输出 | fm² | 八重态道有效耦合系数 | 实数 |
| K08± | K₀₈± | 输出 | fm² | 混合道有效耦合系数 | 实数 |

## 5. 输入参数详细说明

### 5.1 必需参数
- **G**: 基础四夸克相互作用耦合常数，典型值约为5.0×10⁻⁶ MeV⁻²，决定手征对称性自发破缺强度,在src/Constants/Constants_PNJL.jl中定义
- **K**: 't Hooft六夸克相互作用耦合常数，典型值约为1.0×10⁻¹³ MeV⁻⁵，描述U_A(1)反常,在src/Constants/Constants_PNJL.jl中定义
- **\(\phi_u,\phi_s\)**: 物理上的夸克凝聚，通过 A 积分确定
- **G_u,G_s**: 当前实现中传给 `calculate_effective_couplings` 的 helper 变量，满足 `G_u=-\phi_u, G_s=-\phi_s`

`G_u/G_s` 的命名是兼容旧接口的历史遗留。后续完整 KMT 路径应直接复用
平衡解中已得到的 `\phi_u,\phi_d,\phi_s`，避免再次从 `A_f` 重积分或把
`H_f` 当作独立物理输入。

### 5.2 非对称背景下的旧/完整通道映射

由旧 helper 定义 `H_f=-\phi_f` 可直接得到：

```math
K_{123}^{\pm}=G\pm\frac12 K H_s
                  =K_{12}^{P/S},
\qquad
K_{4567}^{\pm}=G\pm\frac12 K H_u
                  =K_{67}^{P/S}.
```

第二个等式说明旧 `K4567` 在 `\phi_u\ne\phi_d` 时是 `K67`（u spectator）
的代数形式，而不是 charged kaon 的耦合。物理通道应按

```text
pi^±                 -> K12
K^+=u bar(s), K^-=s bar(u) -> K45  (d spectator)
K^0=d bar(s), bar(K)^0=s bar(d)  -> K67  (u spectator)
```

在 `\phi_u=\phi_d` 时 `K45=K67=K4567`，这才恢复旧接口把四个 kaon 味道
通道合并的极限。该映射是纯代数兼容关系，不等同于已经完成非对称背景下的
完整 RPA 极化或介子热力学反馈。


### 5.3 可选参数
- **Nc**: 色数，默认值为3，对应QCD的SU(3)规范群

## 6. 输出结果说明

### 6.1 主输出
- **有效耦合系数组**: 包含10个耦合系数的元组，分别对应不同味道通道的标量(+)和赝标量(-)相互作用强度

### 6.2 辅助输出
- **耦合矩阵行列式**: det K = K₀⁺K₈⁺ - (K₀₈⁺)²，用于混合介子传播子计算

## 7. 依赖关系

### 7.1 输入依赖
- **夸克凝聚计算**: 物理口径上先计算
```math
\phi_f = \frac{N_c}{4\pi^2}\, m_f\, A_f(T,\mu),
```
当前 helper `calculate_G_from_A(A_f, m_f)` 返回的则是 `-\phi_f`
- **分布函数**: n_f±依赖于Polyakov环Φ、温度T和化学势μ_f
- **能隙方程求解**: 夸克质量m_f需要通过能隙方程自洽求解

### 7.2 输出用途
- **介子传播子**: 直接用于计算π、K、η、η'等介子的传播子
- **散射截面**: 作为夸克-夸克散射振幅的输入参数
- **驰豫时间**: 通过散射截面影响输运系数的计算
