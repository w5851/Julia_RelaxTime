# 错误处理指南

## 概述

本文档定义QCD模型库的统一错误处理策略，包括异常类型层次结构、各层的异常处理责任、日志记录规范和降级策略。

**适用范围**：所有QCD模型实现（PNJL、rPNJL等）及相关工具模块

---

## 1. 异常类型层次结构

### 1.1 自定义异常类型

```julia
# src/core/exceptions.jl

"""所有QCD模型库异常的根类型"""
abstract type QCDModelException <: Exception end

"""参数验证失败"""
struct ModelParameterError <: QCDModelException
    parameter::Symbol
    value::Any
    reason::String
end

"""数值计算收敛失败"""
struct NumericalConvergenceError <: QCDModelException
    operation::String
    iterations::Int
    tolerance::Float64
    achieved::Float64
end

"""物理约束违反"""
struct PhysicalConstraintError <: QCDModelException
    constraint::String
    values::Dict{Symbol, Any}
end

"""模型配置错误"""
struct ModelConfigurationError <: QCDModelException
    model_type::Type
    message::String
end

"""数值不稳定"""
struct NumericalInstabilityError <: QCDModelException
    location::String
    problematic_value::Any
    context::String
end
```

### 1.2 异常消息格式

所有异常消息应包含：
1. **问题描述**：简洁说明发生了什么
2. **上下文信息**：相关参数值、计算状态
3. **建议操作**：用户可以采取的修正措施

**示例**：
```julia
function Base.showerror(io::IO, e::ModelParameterError)
    print(io, "ModelParameterError: Invalid parameter '$(e.parameter)'\n")
    print(io, "  Value: $(e.value)\n")
    print(io, "  Reason: $(e.reason)\n")
    print(io, "  Suggestion: Check parameter constraints in model documentation")
end
```

---

## 2. 参数验证策略

### 2.1 验证时机

**构造时验证**（推荐）：
```julia
struct PNJLModel <: AbstractIsotropicModel
    params::Dict{Symbol, Any}
    
    function PNJLModel(params::Dict{Symbol, Any})
        validate_pnjl_parameters(params)  # 抛出 ModelParameterError
        new(params)
    end
end
```

**计算前验证**（备选）：
```julia
function calculate_mass_vec(model::PNJLModel, φ::SVector{3, T}) where {T}
    validate_chiral_condensate(φ)  # 快速检查
    # ... 计算逻辑
end
```

### 2.2 参数约束定义

```julia
# src/models/base/validation.jl

"""PNJL模型参数约束"""
const PNJL_CONSTRAINTS = Dict(
    :T => (min=0.0, max=1.0, unit="fm⁻¹", description="Temperature"),
    :mu => (min=0.0, max=2.0, unit="fm⁻¹", description="Chemical potential"),
    :xi => (min=-1.0, max=1.0, unit="dimensionless", description="Anisotropy parameter"),
    :G => (min=0.0, max=Inf, unit="fm²", description="Four-quark coupling"),
    :K => (min=0.0, max=Inf, unit="fm⁵", description="Six-quark coupling"),
)

function validate_parameter(name::Symbol, value, constraints)
    c = get(constraints, name, nothing)
    if c === nothing
        @warn "No constraints defined for parameter :$name"
        return
    end
    
    if value < c.min || value > c.max
        throw(ModelParameterError(
            name, value,
            "Value must be in range [$(c.min), $(c.max)] $(c.unit)"
        ))
    end
end
```

### 2.3 物理约束检查

```julia
"""检查Polyakov loop物理约束：0 ≤ Φ, Φ̄ ≤ 1"""
function validate_polyakov_loops(Φ, Φ̄)
    if !(0 ≤ Φ ≤ 1) || !(0 ≤ Φ̄ ≤ 1)
        throw(PhysicalConstraintError(
            "Polyakov loop bounds",
            Dict(:Φ => Φ, :Φ̄ => Φ̄)
        ))
    end
end

"""检查质量的物理合理性"""
function validate_quark_masses(masses::SVector{3, T}) where {T}
    for (i, m) in enumerate(masses)
        if m < 0
            throw(PhysicalConstraintError(
                "Positive quark mass",
                Dict(:flavor => i, :mass => m)
            ))
        end
        if m > 2.0  # fm⁻¹, ~400 MeV
            @warn "Unusually large quark mass" flavor=i mass=m
        end
    end
end
```

---

## 3. 各层异常处理责任

### 3.1 分层架构

```
用户脚本
    ↓
高层API (workflows, scans)
    ↓
模型接口 (AbstractQCDModel)
    ↓
基础工具 (base utilities)
    ↓
数值库 (Julia stdlib)
```

### 3.2 异常处理责任矩阵

| 层级 | 捕获异常 | 转换异常 | 记录日志 | 向上传播 |
|------|---------|---------|---------|---------|
| **基础工具** | 数值库异常 | → NumericalInstabilityError | DEBUG级别 | ✅ |
| **模型接口** | 基础工具异常 | → ModelConfigurationError | INFO级别 | ✅ |
| **高层API** | 模型接口异常 | → 用户友好消息 | WARN级别 | ✅ |
| **用户脚本** | 所有异常 | - | ERROR级别 | 决定是否继续 |

### 3.3 异常处理示例

**基础工具层**：
```julia
# src/models/base/integrals.jl
function safe_log(x::T) where {T<:Real}
    if x ≤ 0
        @debug "safe_log: non-positive input" x
        if x == 0
            return T(-Inf)  # 数学上合理的回退
        else
            throw(NumericalInstabilityError(
                "safe_log",
                x,
                "Logarithm of negative number"
            ))
        end
    end
    return log(x)
end
```

**模型接口层**：
```julia
# src/models/isotropic/pnjl.jl
function calculate_mass_vec(model::PNJLModel, φ::SVector{3, T}) where {T}
    try
        validate_chiral_condensate(φ)
        # ... 计算逻辑
    catch e
        if e isa PhysicalConstraintError
            @info "Invalid chiral condensate" φ exception=e
            throw(ModelConfigurationError(
                PNJLModel,
                "Failed to calculate masses: $(e.constraint)"
            ))
        else
            rethrow()
        end
    end
end
```

**高层API层**：
```julia
# src/models/scans/TmuScan.jl
function run_tmu_scan(model, T_range, mu_range; kwargs...)
    results = []
    for T in T_range, mu in mu_range
        try
            result = solve_gap_equations(model, T, mu; kwargs...)
            push!(results, result)
        catch e
            if e isa NumericalConvergenceError
                @warn "Convergence failed at (T, μ)" T mu iterations=e.iterations
                push!(results, missing)  # 标记失败点
            else
                rethrow()  # 其他异常向上传播
            end
        end
    end
    return results
end
```

---

## 4. 日志记录规范

### 4.1 日志级别使用指南

| 级别 | 使用场景 | 示例 |
|------|---------|------|
| **DEBUG** | 详细的内部状态、中间计算结果 | "Iteration 5: residual = 1e-8" |
| **INFO** | 正常操作的关键步骤 | "Starting PNJL scan: T ∈ [0.1, 0.3]" |
| **WARN** | 可恢复的异常情况、性能问题 | "Vandermonde term negative: -0.05" |
| **ERROR** | 严重错误，但程序可继续 | "Failed to solve at (T=0.15, μ=0.8)" |

### 4.2 日志消息格式

**推荐格式**：
```julia
@info "Operation description" key1=value1 key2=value2 exception=e
```

**示例**：
```julia
@debug "Calculating thermal contribution" T=T_fm mu=mu_vec xi=xi
@info "Gap equations converged" iterations=iter residual=norm(F)
@warn "Slow convergence detected" iterations=iter max_iterations=max_iter
@error "Numerical instability" location="calculate_U" value=problematic_value
```

### 4.3 敏感数据处理

**不要记录**：
- 完整的大型数组（使用统计信息代替）
- 用户的私有配置路径

**推荐做法**：
```julia
# ❌ 不好：记录完整数组
@debug "Mesh values" mesh=thermal_p_mesh

# ✅ 好：记录统计信息
@debug "Mesh statistics" n=length(thermal_p_mesh) min=minimum(thermal_p_mesh) max=maximum(thermal_p_mesh)
```

---

## 5. 数值计算降级策略

### 5.1 积分不收敛

**策略**：增加节点数或切换积分方法

```julia
function calculate_thermal_contribution(model, masses, T, mu; 
                                       nodes=64, max_nodes=256)
    current_nodes = nodes
    while current_nodes ≤ max_nodes
        try
            result = integrate_with_nodes(model, masses, T, mu, current_nodes)
            return result
        catch e
            if e isa NumericalConvergenceError
                @warn "Integration failed, increasing nodes" current_nodes next_nodes=current_nodes*2
                current_nodes *= 2
            else
                rethrow()
            end
        end
    end
    throw(NumericalConvergenceError(
        "thermal_contribution",
        max_nodes,
        1e-6,
        Inf
    ))
end
```

### 5.2 求解器不收敛

**策略**：调整初值或放宽容差

```julia
function solve_with_fallback(F!, x0, tol=1e-10; max_attempts=3)
    tolerances = [tol, tol*10, tol*100]
    
    for (attempt, current_tol) in enumerate(tolerances)
        try
            result = nlsolve(F!, x0; ftol=current_tol)
            if converged(result)
                if attempt > 1
                    @warn "Converged with relaxed tolerance" attempt tolerance=current_tol
                end
                return result
            end
        catch e
            @debug "Solve attempt failed" attempt exception=e
        end
    end
    
    throw(NumericalConvergenceError(
        "gap_equations",
        max_attempts,
        tol,
        Inf
    ))
end
```

### 5.3 Vandermonde项负值

**策略**：记录警告但继续计算（物理上可能合理）

```julia
function safe_vandermonde(Φ, Φ̄)
    discriminant = 1 - 6Φ*Φ̄ + 4(Φ^3 + Φ̄^3) - 3(Φ*Φ̄)^2
    
    if discriminant < 0
        @warn "Vandermonde discriminant negative" Φ Φ̄ discriminant
        # 不修正值，让调用者决定如何处理
    end
    
    return discriminant
end
```

---

## 6. 用户指南

### 6.1 捕获和处理异常

**基本模式**：
```julia
using QCDModels

try
    model = create_model(:PNJL, params=my_params)
    result = calculate_thermodynamics(model, T, mu)
catch e
    if e isa ModelParameterError
        println("Parameter error: $(e.reason)")
        println("Please check: $(e.parameter) = $(e.value)")
    elseif e isa NumericalConvergenceError
        println("Calculation did not converge")
        println("Try adjusting initial conditions or tolerances")
    else
        rethrow()  # 未预期的异常
    end
end
```

### 6.2 配置日志级别

```julia
using Logging

# 开发/调试：显示所有日志
global_logger(ConsoleLogger(stderr, Logging.Debug))

# 生产：只显示警告和错误
global_logger(ConsoleLogger(stderr, Logging.Warn))

# 保存到文件
io = open("qcd_model.log", "w")
global_logger(SimpleLogger(io, Logging.Info))
```

### 6.3 自定义错误处理

```julia
# 忽略特定警告
function my_scan(...)
    Logging.with_logger(Logging.SimpleLogger(stderr, Logging.Error)) do
        # 这里的WARN级别日志会被抑制
        run_scan(...)
    end
end

# 收集所有错误
errors = []
for point in scan_points
    try
        result = calculate(point)
    catch e
        push!(errors, (point=point, error=e))
    end
end
```

---

## 7. 测试要求

### 7.1 异常测试

每个可能抛出异常的函数都应有测试：

```julia
@testset "Parameter validation" begin
    # 测试无效参数
    @test_throws ModelParameterError PNJLModel(Dict(:T => -0.1))
    @test_throws ModelParameterError PNJLModel(Dict(:mu => 5.0))
    
    # 测试边界值
    @test_nowarn PNJLModel(Dict(:T => 0.0))  # 边界应该允许
    @test_nowarn PNJLModel(Dict(:T => 1.0))
end

@testset "Numerical stability" begin
    # 测试极端输入
    @test_throws NumericalInstabilityError safe_log(-1.0)
    @test safe_log(0.0) == -Inf
    @test safe_log(1e-100) ≈ -230.26  # 应该能处理极小值
end
```

### 7.2 降级策略测试

```julia
@testset "Fallback mechanisms" begin
    # 测试积分降级
    result = calculate_thermal_contribution(model, masses, T, mu; nodes=8)
    @test !isnan(result)
    
    # 测试求解器降级
    result = solve_with_fallback(F!, x0)
    @test converged(result)
end
```

---

## 8. 检查清单

在实现新功能时，请确认：

- [ ] 所有用户输入参数都经过验证
- [ ] 关键计算步骤有适当的日志记录
- [ ] 数值不稳定情况有明确的异常或警告
- [ ] 异常消息包含足够的上下文信息
- [ ] 实现了合理的降级策略
- [ ] 编写了异常情况的测试用例
- [ ] 文档说明了可能的异常和处理方法

---

## 9. 参考资料

- Julia异常处理：https://docs.julialang.org/en/v1/manual/control-flow/#Exception-Handling
- Julia日志系统：https://docs.julialang.org/en/v1/stdlib/Logging/
- 项目代码规范：`docs/dev/开发规范与思考清单.md`
