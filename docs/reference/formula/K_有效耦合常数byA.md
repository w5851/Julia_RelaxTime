# 有效耦合系数 - 实现文档

## 1. 公式标识
- **程序实现**: `K_effective_couplings`
- **对应文件**: `../../src/relaxtime/EffectiveCouplings.jl`

## 2. 物理意义
- **物理背景**: 在3味PNJL模型中，有效耦合系数描述了夸克-反夸克相互作用通过介子传播的有效强度，是随机相位近似下介子传播子的关键参数
- **在计算链中的作用**: 作为介子传播子和散射截面的输入，直接影响驰豫时间和输运系数的计算
- **相关物理量**: 原始耦合常数G、K，夸克凝聚G^μ 、G^s，温度T，化学势μ

## 3. 数学表达式

### 3.1 原始公式
```math
\begin{aligned}
K_0^{\pm} &= G \mp \frac{1}{3} K(2G^\mu + G^s) \\
K_1^{\pm} = K_2^{\pm} = K_3^{\pm} &= G \pm \frac{1}{2} KG^s \\
K_4^{\pm} = K_5^{\pm} = K_6^{\pm} = K_7^{\pm} &= G \pm \frac{1}{2} KG^\mu \\
K_8^{\pm} &= G \pm \frac{1}{6} K(4G^\mu - G^s) \\
K_{08}^{\pm} &= \pm \frac{1}{6} \sqrt{2} K(G^\mu - G^s)
\end{aligned}
```

其中：
```math
\begin{aligned}
G^f &= -\frac{N_c}{4\pi^2} \; m_f\; A_f(T, \mu) \\
A_f(T, \mu) &= 16\pi^2 \int_0^\Lambda \frac{d^3 p}{(2\pi)^3} \left[ 1 - n_f^+ \frac{1}{2E} - n_f^- \frac{1}{2E} \right]
\end{aligned}
```

### 3.2 程序实现形式
```julia
function calculate_effective_couplings(G, K, G_u, G_s, Nc=3)
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
| G_u | G^μ | 输入 | 无量纲 | u,d夸克凝聚相关函数 | 实数 |
| G_s | G^s | 输入 | 无量纲 | s夸克凝聚相关函数 | 实数 |
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
- **G_u,G_d**: u,d夸克凝聚相关函数，通过积分计算得到，依赖温度T和化学势μ,依赖函数src/relaxtime/OneLoopIntegrals.jl中的A函数


### 5.2 可选参数
- **Nc**: 色数，默认值为3，对应QCD的SU(3)规范群

## 6. 输出结果说明

### 6.1 主输出
- **有效耦合系数组**: 包含10个耦合系数的元组，分别对应不同味道通道的标量(+)和赝标量(-)相互作用强度

### 6.2 辅助输出
- **耦合矩阵行列式**: det K = K₀⁺K₈⁺ - (K₀₈⁺)²，用于混合介子传播子计算

## 7. 依赖关系

### 7.1 输入依赖
- **夸克凝聚计算**: 需要先计算 $G^f = -N_c/(4\pi^2)\, m_f\, A_f(T,\mu)$，对应函数为 `calculate_G_from_A(A_f, m_f)`
- **分布函数**: n_f±依赖于Polyakov环Φ、温度T和化学势μ_f
- **能隙方程求解**: 夸克质量m_f需要通过能隙方程自洽求解

### 7.2 输出用途
- **介子传播子**: 直接用于计算π、K、η、η'等介子的传播子
- **散射截面**: 作为夸克-夸克散射振幅的输入参数
- **驰豫时间**: 通过散射截面影响输运系数的计算