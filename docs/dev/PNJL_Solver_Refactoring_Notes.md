# PNJL 求解器重构笔记

> 临时文档：记录 Julia_test 实验代码的可借鉴之处及待解决问题

---

## 一、Julia_test (Pnjl_minimal) 可借鉴的技术点

### 1.1 数值稳定性：Log-Sum-Exp 技巧

**问题**：在计算 Polyakov loop 修正的对数项时，指数函数可能导致数值溢出。

**Julia_test 实现** (`_log_polyakov_series`):
```julia
function _log_polyakov_series(a, Phi1, Phi2)
    m = max.(0.0, a)
    m = max.(m, 2 .* a)
    m = max.(m, 3 .* a)
    term = exp.(-m) .+ 3.0 .* Phi1 .* exp.(a .- m) .+ 
           3.0 .* Phi2 .* exp.(2 .* a .- m) .+ exp.(3 .* a .- m)
    clamped = max.(term, POLYAKOV_EPS)
    return m .+ log.(clamped)
end
```

**当前 AnisoGapSolver 实现** (`calculate_log_term`):
```julia
exp1 = exp(-x_i)
exp2 = exp1 * exp1
exp3 = exp1 * exp2
f1_val = 1.0 + 3.0 * Φ * exp1 + 3.0 * Φ̄ * exp2 + exp3
return safe_log(f1_val) + safe_log(f2_val)
```

**建议**：在极端参数（高温、高化学势）下，当前实现可能出现数值问题。考虑引入 log-sum-exp 技巧。

### 1.2 自适应动量截断策略

**Julia_test 实现** (`analytic_p_max`):
- 使用离散节点策略 [250, 500, 1000, 2000, 4000, 6000] MeV
- 通过 `ForwardDiff.value` 冻结参数，保证 AD 稳定性
- 低温快捷路径：T_mev ≤ (μ - m)/10 时直接返回 1.3 × pF

**当前 AnisoGapSolver**：使用固定上限 (默认 10 fm⁻¹)

**权衡**：
- 固定上限：简单、AD 友好，但可能在低温时浪费计算资源
- 自适应截断：更高效，但需要处理 AD 冻结带来的导数误差（~3-5%）

### 1.3 节点缓存策略

**Julia_test**：双层缓存
- `MINIMAL_LEGENDRE_CACHE`: 标准 [-1, 1] 节点
- `MINIMAL_BOUNDARY_CACHE`: 预定义边界的映射节点

**当前 AnisoGapSolver**：单层缓存 `NODE_CACHE`

**建议**：当前实现已足够，但若引入自适应截断，可参考双层缓存设计。

### 1.4 调试输出与组件分解

**Julia_test** 提供 `Omega_minimal_components` 和 `print_omega_components`：
- 返回各组成部分（chiral, polyakov, vacuum, thermal）
- 便于调试和验证

**建议**：在 AnisoGapSolver 中添加类似的诊断函数。

---

## 二、待解决的架构问题

### 2.1 ImplicitDifferentiation.jl 集成与求解器重构

**现状**：
- `ThermoDerivatives.jl` 已成功使用 ImplicitDifferentiation.jl
- 需要将求解器分离为 `conditions` 和 `forward_solve` 两部分

**ImplicitDifferentiation.jl 要求的接口**：
```julia
# 前向求解：给定参数 θ，返回解 x
function forward_solve(θ::AbstractVector)
    # ... 求解 F(x; θ) = 0
    return (x, nothing)  # 或 (x, auxiliary_data)
end

# 条件函数：返回残差 F(θ, x)
function conditions(θ::AbstractVector, x::AbstractVector, z)
    # ... 计算残差
    return F
end

# 创建隐函数
implicit_solver = ImplicitFunction(forward_solve, conditions)
```

**重构方向**：
1. 将 `residual_mu!` 和 `residual_rho!` 重构为符合 ImplicitDifferentiation 接口的形式
2. 分离"条件定义"和"求解执行"
3. 统一 fixed-μ 和 fixed-ρ 两种模式的接口

**建议的新结构**：
```
src/pnjl/
├── core/                      # 核心计算（无求解逻辑）
│   ├── Thermodynamics.jl      # Ω, P, s, ε 计算
│   ├── GapEquations.jl        # 能隙方程残差定义
│   └── MassCalculation.jl     # 质量计算
├── solvers/                   # 求解器
│   ├── Conditions.jl          # conditions 函数定义
│   ├── ForwardSolve.jl        # forward_solve 实现
│   ├── ImplicitSolver.jl      # ImplicitFunction 封装
│   └── LegacySolver.jl        # 兼容旧接口（可选）
├── analysis/                  # 分析工具
│   ├── ThermoDerivatives.jl   # 热力学导数
│   ├── PhaseTransition.jl     # 相变分析
│   └── CEPFinder.jl           # 临界点搜索
├── scans/                     # 扫描工具
│   └── ...
└── PNJL.jl                    # 主模块入口
```

### 2.2 一阶相变处理

**问题描述**（参见 `PNJL模型一阶相变处理策略.md`）：
- 在一阶相变区域，能隙方程存在多值解
- 需要选择热力学势最小的物理稳定解
- 相变线位置需要通过 Ω 比较确定

**当前状态**：代码未处理一阶相变问题

**策略：双分支逆向扫描**
1. 低化学势分支：从 μ_min 开始，使用强子相初值，沿 μ 增加方向扫描
2. 高化学势分支：从 μ_max 开始，使用夸克相初值，沿 μ 减小方向扫描
3. 对每个 μ，比较两分支的 Ω，选择较小者
4. 相变点：Ω_low(μ_c) = Ω_high(μ_c)

**实现要点**：
- 需要可靠的分支初值（强子相 vs 夸克相）
- 需要连续性跟踪（前一点的解作为下一点的初值）
- 需要 Ω 计算的一致性

### 2.3 初解选择问题

**问题**：给定任意 (T, μ) 点，如何选择合适的初解？

**挑战**：
1. 不知道该点与一阶相变线的相对位置
2. 错误的初解可能导致收敛到亚稳态或不收敛
3. 单点求解时无法使用连续性跟踪

**可能的解决方案**：

#### 方案 A：多初值尝试
```julia
function solve_with_multiple_seeds(T, μ)
    seeds = [
        hadron_phase_seed(T, μ),   # 强子相典型初值
        quark_phase_seed(T, μ),    # 夸克相典型初值
        interpolated_seed(T, μ),   # 基于已知解的插值
    ]
    
    results = []
    for seed in seeds
        res = try_solve(T, μ, seed)
        res.converged && push!(results, res)
    end
    
    # 选择 Ω 最小的解
    return argmin(r -> r.Omega, results)
end
```

#### 方案 B：种子缓存与插值
- 维护已求解点的缓存
- 对新点，从最近的已知解插值得到初值
- 参见 `SeedCache.jl` 的现有实现

#### 方案 C：相图预计算
- 预先计算粗网格上的相图
- 确定一阶相变线的大致位置
- 根据 (T, μ) 与相变线的关系选择初值

#### 方案 D：连续性路径
- 从已知收敛点出发
- 沿参数空间路径逐步移动到目标点
- 每步使用前一点的解作为初值

**建议**：结合方案 A 和 B，实现"多初值 + 缓存"策略。

---

## 三、优先级排序

### 高优先级
1. **求解器重构**：分离 conditions/forward_solve，统一 ImplicitDifferentiation 接口
2. **初解策略**：实现多初值尝试 + 种子缓存

### 中优先级
3. **一阶相变处理**：实现双分支扫描策略
4. **数值稳定性**：引入 log-sum-exp 技巧

### 低优先级
5. **自适应截断**：评估是否值得引入
6. **诊断工具**：添加 Ω 组件分解函数

---

## 四、具体实施计划

### Phase 1：求解器重构（预计 2-3 天）

1. 创建 `src/pnjl/core/` 目录，提取核心计算函数
2. 定义统一的 `conditions` 接口
3. 实现 `forward_solve` 封装
4. 创建 `ImplicitSolver.jl`，统一 fixed-μ 和 fixed-ρ 模式
5. 更新 `ThermoDerivatives.jl` 使用新接口
6. 添加测试确保行为一致

### Phase 2：初解策略（预计 1-2 天）

1. 定义强子相/夸克相典型初值函数
2. 实现多初值尝试逻辑
3. 增强 `SeedCache.jl` 的插值能力
4. 添加 Ω 比较选择逻辑

### Phase 3：一阶相变处理（预计 2-3 天）

1. 实现双分支扫描函数
2. 添加相变点检测逻辑
3. 集成到 `TmuScan.jl` 和 `TrhoScan.jl`
4. 添加相变线输出

---

## 五、重构设计讨论（2025-12-28 更新）

### 5.1 目录结构修订

经讨论后的目录结构：

```
src/pnjl/
├── core/
│   ├── Thermodynamics.jl      # Ω, P, s, ε, ρ, M 计算（含质量计算）
│   └── Integrals.jl           # 真空项、热项积分
├── solver/
│   ├── ConstraintModes.jl     # FixedMu, FixedRho, FixedEntropy 等类型定义
│   ├── Conditions.jl          # 条件定义（gap_conditions + 各模式约束）
│   └── ImplicitSolver.jl      # ImplicitFunction 封装 + 初值策略
├── derivatives/
│   └── ThermoDerivatives.jl   # 热力学导数（使用 ImplicitDifferentiation）
├── analysis/
│   ├── PhaseTransition.jl     # 相变分析
│   ├── CEPFinder.jl           # 临界点搜索
│   └── MaxwellRhoMu.jl        # Maxwell 构造
├── scans/
│   └── ...
├── SeedCache.jl               # 种子缓存（保留现有功能）
└── PNJL.jl                    # 主模块入口
```

**关键设计决策**：
- `ThermoDerivatives.jl` 从 `analysis/` 移到独立的 `derivatives/` 目录
- 质量计算合并到 `Thermodynamics.jl`，不单独成模块
- 初值策略暂时放在 `ImplicitSolver.jl` 中，复杂化后再抽取

### 5.2 多求解模式的条件设计

**共同部分**：5 个能隙方程
```math
\frac{\partial\Omega}{\partial \phi_u} = \frac{\partial\Omega}{\partial \phi_d} = \frac{\partial\Omega}{\partial \phi_s} = \frac{\partial\Omega}{\partial \Phi} = \frac{\partial\Omega}{\partial \bar\Phi} = 0
```

**不同模式的约束**：

| 模式 | 已知参数 | 未知参数 | 额外约束 |
|------|----------|----------|----------|
| FixedMu | T, μ | φ(3), Φ, Φ̄ | 无（5方程5未知） |
| FixedRho | T, ρ | φ(3), Φ, Φ̄, μ(3) | μ_u=μ_d=μ_s, ρ(μ)=ρ_target |
| FixedEntropy | T, s | φ(3), Φ, Φ̄, μ(3) | μ_u=μ_d=μ_s, s(μ)=s_target |
| FixedSigma | T, σ=s/ρ | φ(3), Φ, Φ̄, μ(3) | μ_u=μ_d=μ_s, s/ρ=σ_target |

**实现模式**：组合式条件构建

```julia
# ConstraintModes.jl
abstract type ConstraintMode end
struct FixedMu <: ConstraintMode end
struct FixedRho <: ConstraintMode 
    rho_target::Float64
end
struct FixedEntropy <: ConstraintMode
    s_target::Float64
end
struct FixedSigma <: ConstraintMode
    sigma_target::Float64
end

# Conditions.jl
# 核心能隙条件（所有模式共用）
function gap_conditions(x_state, params)
    # 返回 5 个能隙方程的残差
end

# 根据模式构建完整条件
function build_conditions(mode::FixedMu)
    return (θ, x) -> gap_conditions(x, make_params(θ, mode))
end

function build_conditions(mode::FixedRho)
    return (θ, x) -> begin
        gap = gap_conditions(x[1:5], make_params(θ, x[6:8], mode))
        mu_eq1 = x[6] - x[7]
        mu_eq2 = x[7] - x[8]
        rho_constraint = compute_rho(x, θ) - mode.rho_target
        return [gap..., mu_eq1, mu_eq2, rho_constraint]
    end
end
```

### 5.3 forward_solve 与 conditions 的关系

`forward_solve` 和 `conditions` 共享同一个条件定义：

```julia
# 条件定义（核心）
function make_conditions(mode::ConstraintMode)
    return (θ, x) -> ... # 返回残差向量
end

# forward_solve 使用相同的条件
function make_forward_solve(mode::ConstraintMode)
    cond_fn = make_conditions(mode)
    return θ -> begin
        f! = (F, x) -> (F .= cond_fn(θ, x))
        res = nlsolve(f!, initial_guess(mode, θ); ...)
        return (res.zero, nothing)
    end
end

# ImplicitFunction 的 conditions 包装
function make_implicit_conditions(mode::ConstraintMode)
    cond_fn = make_conditions(mode)
    return (θ, x, z) -> cond_fn(θ, x)
end
```

---

## 六、参考资料

- `docs/reference/formula/pnjl/Omega_各向同性.md` - 各向同性 PNJL 公式
- `docs/reference/formula/pnjl/Omega_RS各向异性.md` - 各向异性 PNJL 公式
- `docs/reference/formula/pnjl/PNJL模型一阶相变处理策略.md` - 一阶相变策略
- `Julia_test/integral_test/Pnjl_minimal.jl` - 实验性实现参考
- ImplicitDifferentiation.jl 文档：https://gdalle.github.io/ImplicitDifferentiation.jl/

---

*创建日期：2025-12-28*
*最后更新：2025-12-28*
*状态：设计讨论中*
