# 散射矩阵元公式文档（修订版）

## 1. 公式标识
- **程序实现**: `scattering_amplitude`
- **对应文件**: `src/relaxtime/ScatteringAmplitude.jl`

## 2. 物理意义
- **物理背景**: 在三味PNJL模型中，夸克-夸克、夸克-反夸克弹性散射通过标量和赝标量介子交换传递相互作用
- **在计算链中的作用**: 连接介子传播子和散射截面，是计算弛豫时间和输运系数的核心输入
- **重要说明**: **每个散射过程独立计算其散射矩阵元，不同过程之间没有耦合**

## 3. 计算结构说明

### 3.1 散射过程的独立性
在程序中，散射矩阵元函数应该是**从散射过程到振幅的映射**：

```julia
function scattering_amplitude(process::Symbol, s, t, T, μ_q, masses)
    # 每个散射过程独立计算
    # process ∈ [:uu_to_uu, :ud_to_ud, :us_to_us, :uubar_to_uubar, ...]
end
```

### 3.2 同一过程中的散射道求和
**重要澄清**: 公式5.16中提到的"所有可能的散射道求和"指的是**同一散射过程中**的不同费曼图道（如t道、u道、s道），而不是不同散射过程之间的求和。

例如：
- 对于 `uu → uu` 过程：`M_total = M_t + M_u` （t道 + u道）
- 对于 `uū → uū` 过程：`M_total = M_s + M_t` （s道 + t道）

## 4. 数学表达式

### 4.1 Mandelstam变量定义
```
s = (p₁ + p₂)²
t = (p₁ - p₃)²  
u = (p₁ - p₄)²
约束条件: s + t + u = m₁² + m₂² + m₃² + m₄²
```

### 4.2 夸克-夸克散射（表5.1中的独立过程）

**过程1: uu → uu**
```
M = M_t - M_u  （注意符号约定）
1/(4N_c²) ∑|M|² = 1/(4N_c²) ∑|M_t|² + 1/(4N_c²) ∑|M_u|² - 2/(4N_c²) ∑(M_t M_u^*)
```

**过程2: ss → ss**
```
M = M_t - M_u
涉及介子: t道和u道的η, η', σ, σ'
```

**过程3: ud → ud**  
**过程4: us → us**

### 4.3 夸克-反夸克散射（表5.2中的独立过程）

**过程1: uū → uū**
```
M = M_s - M_t
1/(4N_c²) ∑|M|² = 1/(4N_c²) ∑|M_s|² + 1/(4N_c²) ∑|M_t|² - 2/(4N_c²) ∑(M_s M_t^*)
```

**过程2: sš → sš**  
**过程3: uū → dđ**  
**过程4: uū → sš**  
**过程5: sš → uū**

## 5. 参数说明表

| 参数 | 类型 | 物理意义 | 单位 | 备注 |
|------|------|----------|------|------|
| `process` | Symbol | 散射过程标识 | 无 | 从预定义列表中选择 |
| `s` | Float64 | Mandelstam s变量 | MeV² | 质心系能量平方 |
| `t` | Float64 | Mandelstam t变量 | MeV² | 动量转移平方 |
| `T` | Float64 | 温度 | MeV | 热力学参数 |
| `μ_q` | Float64 | 夸克化学势 | MeV | 化学势参数 |
| `masses` | Dict | 夸克质量字典 | MeV | `{:u => m_u, :d => m_d, :s => m_s}` |

## 6. 程序接口设计

### 6.1 散射过程枚举
```julia
const QQ_PROCESSES = [
    :uu_to_uu,     # u + u → u + u
    :ss_to_ss,     # s + s → s + s  
    :ud_to_ud,     # u + d → u + d
    :us_to_us      # u + s → u + s
]

const QBARQ_PROCESSES = [
    :uubar_to_uubar,  # u + ū → u + ū
    :ssbar_to_ssbar,  # s + š → s + š
    :uubar_to_ddbar,  # u + ū → d + đ
    :uubar_to_ssbar,  # u + ū → s + š
    :ssbar_to_uubar   # s + š → u + ū
]
```

### 6.2 主要函数签名
```julia
function scattering_amplitude(process::Symbol, s::Float64, t::Float64, 
                             T::Float64, μ_q::Float64, masses::Dict)::Complex{Float64}
    # 根据process选择对应的计算路径
    # 返回该过程的复数散射振幅
end

function scattering_amplitude_squared(process::Symbol, s::Float64, t::Float64,
                                    T::Float64, μ_q::Float64, masses::Dict)::Float64
    # 计算 |M|²，包括色和自旋平均
end
```

## 7. 依赖关系

- **输入依赖**: 
  - `MesonPropagator.get_propagators(process, s, t, T, μ_q)` - 获取特定过程所需的介子传播子
  - `FlavorFactors.get_vertex_factors(process)` - 获取顶点味因子

- **输出用途**:
  - `CrossSection.dsigma_dt(process, s, t, ...)` - 计算微分散射截面
  - `RelaxationTime.scattering_rate(process, ...)` - 计算散射率

## 8. 重要注意事项

1. **过程独立性**: 每个散射过程的矩阵元计算完全独立，不同过程之间没有耦合关系
2. **道求和规则**: 同一过程中的不同散射道（t、u、s道）需要按规则求和
3. **符号约定**: 注意夸克-夸克散射和夸克-反夸克散射的振幅符号差异
4. **味因子处理**: 不同过程的顶点味因子不同，需要根据表5.3正确处理

这样的设计确保了每个散射过程的计算逻辑清晰，便于调试和验证，同时也符合物理上的独立性原则。