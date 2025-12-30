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
│   ├── SeedStrategies.jl      # 初值策略（多种策略的统一接口）
│   └── ImplicitSolver.jl      # ImplicitFunction 封装
├── derivatives/
│   └── ThermoDerivatives.jl   # 热力学导数（使用 ImplicitDifferentiation）
├── analysis/
│   ├── PhaseTransition.jl     # 相变分析
│   ├── CEPFinder.jl           # 临界点搜索
│   └── MaxwellRhoMu.jl        # Maxwell 构造
├── scans/
│   └── ...
└── PNJL.jl                    # 主模块入口
```

**关键设计决策**：
- `ThermoDerivatives.jl` 从 `analysis/` 移到独立的 `derivatives/` 目录
- 质量计算合并到 `Thermodynamics.jl`，不单独成模块
- `SeedStrategies.jl` 作为独立模块，整合现有 `SeedCache.jl` 功能

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

### 5.4 SeedStrategies.jl 初值策略模块设计

**目标**：统一管理"如何为求解器选择初值"的问题，支持多种策略及其组合。

#### 关于 SeedCache.jl 的决策

经分析，现有的 `SeedCache.jl` 存在以下问题：
- 依赖外部预计算数据文件，增加维护成本
- 200+ 行代码实现相对简单的功能（过度工程）
- 实际使用场景有限：扫描用连续性跟踪更有效，单点用多初值尝试更可靠

**决策**：移除 `SeedCache.jl`，用更简单的策略替代。

#### 简化后的策略设计

```julia
# 初值策略的抽象类型
abstract type SeedStrategy end

# 统一接口
function get_seed(strategy::SeedStrategy, θ::AbstractVector, mode::ConstraintMode)::Vector{Float64}
```

#### 具体策略实现

```julia
# 策略1：固定默认值（基于物理直觉）
struct DefaultSeed <: SeedStrategy
    hadron_seed::Vector{Float64}  # 强子相典型值
    quark_seed::Vector{Float64}   # 夸克相典型值
    phase_hint::Symbol            # :hadron, :quark, 或 :auto
end

# 内置默认值
const HADRON_SEED_5 = [-0.07, -0.07, -1.87, 0.02, 0.02]
const QUARK_SEED_5 = [-1.8, -1.8, -2.1, 0.8, 0.8]

# :auto 模式的启发式规则
function auto_phase_hint(T_fm, μ_fm)
    T_mev = T_fm * 197.327
    μ_mev = μ_fm * 197.327
    return (T_mev > 150 || μ_mev > 300) ? :quark : :hadron
end

# 策略2：多初值尝试（用于处理多值解）
struct MultiSeed <: SeedStrategy
    candidates::Vector{SeedStrategy}
    selector::Function  # (results) -> best_result
end

# 默认选择器：选 Ω 最小的收敛解
default_omega_selector(results) = argmin(r -> r.converged ? r.omega : Inf, results)

# 策略3：连续性跟踪（用于参数扫描）
mutable struct ContinuitySeed <: SeedStrategy
    previous_solution::Union{Nothing, Vector{Float64}}
    fallback::SeedStrategy
end

# 策略4：基于相图的智能选择（用于一阶相变区域）
struct PhaseAwareSeed <: SeedStrategy
    phase_boundary::Function  # (T) -> μ_c
    hadron_strategy::SeedStrategy
    quark_strategy::SeedStrategy
end
```

#### 辅助函数

```julia
# 根据求解模式扩展基础种子
function extend_seed(base_seed::Vector{Float64}, mode::FixedMu)
    return base_seed  # 5 维
end

function extend_seed(base_seed::Vector{Float64}, mode::FixedRho)
    # 扩展为 8 维，μ 初值基于 ρ 估计
    μ_guess = estimate_mu_from_rho(mode.rho_target)
    return [base_seed..., μ_guess, μ_guess, μ_guess]
end

# 基于 ρ 估计 μ 的简单公式（自由费米气体近似）
function estimate_mu_from_rho(rho_norm)
    # ρ/ρ₀ ≈ (μ/μ₀)³ 的粗略近似
    μ₀ = 1.5  # fm⁻¹，约 300 MeV
    return μ₀ * cbrt(max(rho_norm, 0.1))
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
*最后更新：2025-12-30*
*状态：✅ 重构完成*

> **注意**：后续开发任务已迁移到 [`PNJL_Phase_Diagram_TODO.md`](./PNJL_Phase_Diagram_TODO.md)

## 七、重构完成记录

### 已完成（2025-12-29）

1. ✅ 创建目录结构
   - `src/pnjl/core/` - 核心计算
   - `src/pnjl/solver/` - 求解器
   - `src/pnjl/derivatives/` - 导数计算

2. ✅ Core 模块
   - `Integrals.jl` - 积分节点缓存、真空项、热项计算
   - `Thermodynamics.jl` - 质量、势能、压强、热力学量计算

3. ✅ Solver 模块
   - `ConstraintModes.jl` - FixedMu, FixedRho, FixedEntropy, FixedSigma 类型
   - `SeedStrategies.jl` - DefaultSeed, MultiSeed, ContinuitySeed, PhaseAwareSeed
   - `Conditions.jl` - gap_conditions, build_conditions, build_residual!
   - `ImplicitSolver.jl` - solve(), create_implicit_solver()

4. ✅ Derivatives 模块
   - `ThermoDerivatives.jl` - mass_derivatives, thermo_derivatives, bulk_viscosity_coefficients
   - 修复了 bulk_viscosity_coefficients 的 Dual 类型泄漏问题（使用有限差分替代 ForwardDiff.jacobian）

5. ✅ PNJL.jl 主模块更新
   - 新接口导出
   - 旧接口兼容（AnisoGapSolver）
   - SeedCache.jl 已弃用

6. ✅ 初值常量更新（基于 trho_seed_table.csv 数据分析）
   - 数据来源：`data/raw/pnjl/seeds/trho_seed_table.csv`
   - 分析脚本：`scripts/analyze_seeds.jl`

### 待完成

> 已迁移到 [`PNJL_Phase_Diagram_TODO.md`](./PNJL_Phase_Diagram_TODO.md)

### 已完成（2025-12-29 续）

7. ✅ 扫描模块迁移
   - `TmuScan.jl` - T-μ 参数空间扫描（使用新求解器架构）
   - `TrhoScan.jl` - T-ρ 参数空间扫描（使用新求解器架构）
   - 特性：
     - 连续性跟踪初值策略
     - 多初值尝试处理收敛困难点
     - 断点续扫支持
     - 输出增加有效质量列 (M_u, M_d, M_s)
   - 测试脚本：`scripts/test_scan_migration.jl`

8. ✅ 双分支扫描与一阶相变处理
   - `DualBranchScan.jl` - 双分支扫描模块
   - 功能：
     - 强子分支：从低 μ 向高 μ 扫描
     - 夸克分支：从高 μ 向低 μ 扫描
     - 自动选择物理分支（Ω 最小）
     - 输出三个关键化学势：
       - μ_c: 相变点（Ω 交叉点）
       - μ_hadron_spinodal: 强子相亚稳态边界
       - μ_quark_spinodal: 夸克相亚稳态边界
     - 解跳跃检测（防止追踪到非物理解）
     - 同解检测（区分真正的共存区）
   - 导出函数：
     - `run_dual_branch_scan()` - 单温度双分支扫描
     - `find_phase_transition()` - 查找相变点和 spinodal 边界
     - `merge_branches()` - 合并输出物理解
     - `scan_phase_diagram()` - 多温度相变线扫描
   - 测试结果（xi=0.0）：
     - T = 50 MeV: μ_c ≈ 356 MeV, spinodal: [350, 365] MeV
     - T = 80 MeV: μ_c ≈ 345 MeV, spinodal: [340, 350] MeV
     - T = 100 MeV: μ_c ≈ 331 MeV, spinodal: [330, 335] MeV
     - T = 120 MeV: 接近 CEP，共存区极窄
   - 文档：`docs/notes/pnjl/一阶相变双分支扫描策略.md`
   - 测试脚本：`scripts/test_dual_branch_scan.jl`

9. ✅ 一阶相变方法对比文档
   - 文档：`docs/notes/pnjl/一阶相变方法对比.md`
   - 对比内容：
     - 方法 A：T-ρ 扫描 + Maxwell 等面积构造
     - 方法 B：双分支 T-μ 扫描 + Ω 比较
   - 精度对比结论：
     - 两种方法在 μ_c 计算上精度一致（误差 < 0.01 MeV）
     - T-ρ 方法直接输出密度边界（ρ_gas, ρ_liquid）
     - T-μ 方法直接输出化学势边界（μ_spinodal）
   - 效率对比：
     - T-ρ 方法：8 维求解，计算量较大
     - T-μ 方法：5 维求解 × 2，效率更高
   - 适用场景建议：
     - T-μ 相图 → 双分支扫描
     - T-ρ 相图 → Maxwell 构造
     - 快速扫描 → 双分支扫描

---

## 八、初值常量说明

### 数据来源

初值常量基于 `data/raw/pnjl/seeds/trho_seed_table.csv` 中的收敛解数据分析得出。
该数据集包含 303 个点，覆盖：
- T: 50 - 200 MeV
- μ: 0 - 376 MeV
- ρ/ρ₀: 0 - 3.0
- ξ = 0.0（各向同性）

### 初值常量定义

```julia
# 强子相（T=50 MeV, μ≈0, ρ=0）
# 特征：手征凝聚完整，Polyakov loop 接近零（禁闭相）
const HADRON_SEED_5 = [-1.84329, -1.84329, -2.22701, 1.0e-5, 4.0e-5]

# 中等密度（T=100 MeV, ρ=1.0）
# 特征：部分手征恢复，Polyakov loop 中等
const MEDIUM_SEED_5 = [-1.3647, -1.3647, -2.14502, 0.10594, 0.15569]

# 高密度（T=100 MeV, ρ=3.0）
# 特征：手征对称性大部分恢复，Polyakov loop 较高
const HIGH_DENSITY_SEED_5 = [-0.21695, -0.21695, -2.01372, 0.18601, 0.22333]

# 高温（T=200 MeV, ρ=0）
# 特征：高温解禁闭相，Polyakov loop 接近 0.6
const HIGH_TEMP_SEED_5 = [-0.73192, -0.73192, -1.79539, 0.60532, 0.60532]

# 8 维版本（含化学势，单位 fm⁻¹）
const HADRON_SEED_8 = [HADRON_SEED_5..., 0.22367, 0.22367, 0.22367]
const MEDIUM_SEED_8 = [MEDIUM_SEED_5..., 1.70267, 1.70267, 1.70267]
const HIGH_DENSITY_SEED_8 = [HIGH_DENSITY_SEED_5..., 1.7516, 1.7516, 1.7516]
```

### 使用建议

1. **低温低密度区域**：使用 `HADRON_SEED_5`
2. **中等密度区域**：使用 `MEDIUM_SEED_5`
3. **高密度区域**：使用 `HIGH_DENSITY_SEED_5`
4. **高温区域**：使用 `HIGH_TEMP_SEED_5`（或 `QUARK_SEED_5`）

对于未知区域，建议使用 `MultiSeed` 策略尝试多个初值。


---

## 九、计算工作流规划（2025-12-29）

### 9.1 方法选择决策

经过对比分析，确定以下技术路线：

| 任务 | 选择的方法 | 理由 |
|------|-----------|------|
| 一阶相变线计算 | T-ρ 扫描 + Maxwell 构造 | 直接输出密度边界，与相图绘制一致 |
| T-μ 扫描初解选择 | 基于预生成的 boundary.csv | 效率高，单次扫描，逻辑简单 |
| 双分支扫描 | 保留作为备用 | 用于验证或无预生成数据时 |

### 9.2 数据生成工作流

```
Step 1: 相结构计算
├── CEP 搜索
├── 一阶相变线 (T-ρ + Maxwell)
├── Crossover 线
└── 输出: cep.csv, boundary.csv, crossover.csv
        ↓
Step 2: 相图绘制
└── 输出: phase_diagram_*.png
        ↓
Step 3: 参数空间扫描
├── T-μ 扫描 (依赖 boundary.csv)
├── T-ρ 扫描
└── 输出: tmu_scan.csv, trho_scan.csv
        ↓
Step 4: 热力学导数计算
└── 输出: derivatives.csv, bulk_viscosity.csv
```

**关键依赖**：T-μ 扫描必须在相结构计算之后，因为需要 boundary.csv 来选择正确的初解。

### 9.3 GitHub Actions 自动化计划

**目标**：
- 用户可在 GitHub 网站手动触发计算（workflow_dispatch）
- 通过 UI 输入扫描参数（ξ, T 范围, μ/ρ 范围等）
- 使用 GitHub 托管的运行器执行计算
- 结果作为 Artifacts 上传，可下载

**优势**：
- 无需本地 Julia 环境
- 随时随地可运行
- 计算资源由 GitHub 提供
- 结果自动保存和版本化

**配置文件**：`.github/workflows/pnjl-pipeline.yml`（待创建）

### 9.4 待实现任务

1. **PhaseAwareSeed 策略**：在 `SeedStrategies.jl` 中实现基于 boundary.csv 的初解选择
2. **TmuScan 更新**：集成 PhaseAwareSeed，支持加载相变线数据
3. **相结构计算脚本**：`scripts/calculate_phase_structure.jl`
4. **GitHub Actions 配置**：完整的 pipeline workflow

### 9.5 相关文档

- 详细工作流说明：`docs/notes/pnjl/PNJL计算工作流.md`
- 一阶相变方法对比：`docs/notes/pnjl/一阶相变方法对比.md`
- 双分支策略（备用）：`docs/notes/pnjl/一阶相变双分支扫描策略.md`


### 9.6 已完成（2025-12-29）

10. ✅ 数据目录重组
    - 创建 `data/reference/pnjl/` 目录存放代码依赖的预计算数据
    - 合并不同 ξ 值的数据到单个文件（通过 xi 列区分）
    - 新文件：
      - `data/reference/pnjl/boundary.csv` - 相变线数据（xi=0.0, 0.2, 0.4）
      - `data/reference/pnjl/cep.csv` - CEP 数据（xi=0.0, 0.2, 0.4）

11. ✅ PhaseAwareSeed 策略实现
    - 位置：`src/pnjl/solver/SeedStrategies.jl`
    - 新增类型：
      - `PhaseBoundaryData` - 相变线数据结构
      - `PhaseAwareSeed` - 单点相变感知策略
      - `PhaseAwareContinuitySeed` - 连续性 + 相变感知策略
    - 新增函数：
      - `load_phase_boundary(xi)` - 从 CSV 加载相变线数据
      - `interpolate_mu_c(data, T)` - 线性插值获取 μ_c(T)
      - `get_phase_hint(strategy, T, μ)` - 获取相位提示
    - 功能：
      - 自动加载 `data/reference/pnjl/boundary.csv` 和 `cep.csv`
      - 连续性跟踪 + 跨越相变线时自动切换初值
      - μ < μ_c(T) → hadron 初值
      - μ > μ_c(T) → quark 初值
      - T > T_CEP → crossover 策略
    - 测试脚本：`scripts/test_phase_aware_seed.jl`

12. ✅ TmuScan 集成 PhaseAwareContinuitySeed
    - 更新 `src/pnjl/scans/TmuScan.jl`
    - 新增参数：`use_phase_aware=true`（默认启用）
    - 功能：
      - 为每个 xi 值创建独立的 PhaseAwareContinuitySeed 跟踪器
      - 每个新温度重置跟踪器
      - 跨越相变线时自动切换初值
      - 保留普通连续性种子作为回退
    - 测试脚本：`scripts/test_tmu_scan_phase_aware.jl`
    - 测试结果：
      - T=100 MeV, μ_c≈331 MeV
      - μ < 331 MeV: 6/6 点正确识别为 hadron 相
      - μ ≥ 331 MeV: 5/5 点正确识别为 quark 相

### 9.7 下一步任务

1. **改进相结构计算脚本**
   - 改进 S 形检测算法
   - 改进 Maxwell 构造算法
   - 或者集成现有的 analysis 模块

2. **GitHub Actions Pipeline**
   - 创建 `.github/workflows/pnjl-pipeline.yml`

### 9.8 已完成（2025-12-29 续）

13. ✅ 相结构计算脚本框架
    - 位置：`scripts/calculate_phase_structure.jl`
    - 功能：
      - Step 1: T-ρ 扫描
      - Step 2: CEP 搜索（简化版）
      - Step 3: Maxwell 构造（简化版）
      - Step 4: 保存结果到 data/reference/pnjl/
    - 命令行参数：
      - `--xi`: 各向异性参数
      - `--T_min/T_max/T_step`: 温度范围
      - `--rho_max/rho_step`: 密度范围
      - `--skip_trho`: 跳过扫描使用已有数据
      - `--verbose`: 详细输出
    - 状态：基本框架完成，S 形检测和 Maxwell 构造算法需要改进

14. ✅ S 形检测和 Maxwell 构造算法改进
    - 更新 `scripts/calculate_phase_structure.jl`
    - 改进内容：
      - S 形检测：基于 `src/pnjl/analysis/PhaseTransition.jl` 的完善实现
        - 精确的 spinodal 点定位（斜率变号位置）
        - 符号变化计数
        - 更准确的 +1 → -1 → +1 模式检测
      - Maxwell 构造：基于 `src/pnjl/analysis/MaxwellRhoMu.jl` 的完善实现
        - 从 spinodal 点估计 μ 搜索区间
        - 区间收缩避免边界问题
        - 符号变化点搜索
        - 二分法精确求解等面积点
        - 详细的失败原因报告
    - 测试结果（xi=0.0, T=50-140 MeV）：
      - 7/10 温度点成功（T=50-80, 110-130 MeV）
      - 3/10 失败（T=90-100, 140 MeV）
        - T=90-100 MeV：求解器在低密度区域失败（NaN），导致曲线不完整
        - T=140 MeV：已超过 CEP，无 S 形（正常）
    - 输出示例：
      ```
      T=50.0 MeV: μ_c=356.46 MeV, ρ_gas=0.085, ρ_liquid=2.591
      T=80.0 MeV: μ_c=345.24 MeV, ρ_gas=0.242, ρ_liquid=2.505
      T=130.0 MeV: μ_c=292.94 MeV, ρ_gas=1.652, ρ_liquid=2.198
      ```

### 9.9 已知问题

1. ~~**T=90-105 MeV 区域求解失败**~~ ✅ 已解决
   - 原因：ρ=0 是奇异点（μ=0 的真空解），从这里开始连续性跟踪会失败
   - 解决方案：使用反向扫描（ρ_max → 0），在 `TrhoScan.jl` 中添加 `reverse_rho=true` 参数
   - 测试结果：T=90 MeV 从 8/16 成功提升到 16/16 成功

2. ~~**CEP 估计精度**~~ ✅ 已解决
   - 实现：二分法细化 CEP 位置
   - 通过线性插值获取中间温度的曲线
   - 测试结果：
     - 初始区间 [130.0, 135.0] MeV
     - 3 次二分后 [130.62, 131.25] MeV
     - T_CEP ≈ 130.94 MeV, μ_CEP ≈ 291.79 MeV
     - 与文献值（T_CEP ≈ 131 MeV, μ_CEP ≈ 291 MeV）吻合

### 9.10 后续任务

> 已迁移到 [`PNJL_Phase_Diagram_TODO.md`](./PNJL_Phase_Diagram_TODO.md)，包括：
> - Crossover 线计算
> - T-μ 扫描脚本
> - 热力学导数计算脚本
> - GitHub Actions Pipeline
> - 完整测试套件
> - 数值稳定性改进（Log-Sum-Exp）
> - 自适应动量截断
> - 诊断工具（Ω 组件分解）

### 9.11 已完成（2025-12-29 续）

15. ✅ TrhoScan 反向扫描
    - 添加 `reverse_rho=true` 参数（默认启用）
    - 解决了 ρ=0 奇异点导致的连续性跟踪失败问题
    - 扫描成功率从 95% 提升到 100%

16. ✅ CEP 二分法细化
    - 默认精度 0.01 MeV
    - 通过线性插值获取中间温度的曲线
    - 结果：T_CEP ≈ 131.04 MeV, μ_CEP ≈ 291.14 MeV (xi=0.0)

17. ✅ 脚本重组和绘图功能
    - 移动 `calculate_phase_structure.jl` 到 `scripts/pnjl/`
    - 创建 `scripts/pnjl/plot_phase_diagram.py`（Python 绘图）
    - 支持 T-μ、T-ρ、组合相图
    - 支持多 ξ 值叠加显示
    - 输出格式：png/pdf/svg
