# 项目级 Solver 时序图（workflow -> solver -> physics -> result）

更新日期：2026-04-06

## 1. 主时序（当前实现视角）

```mermaid
sequenceDiagram
    autonumber
    participant W as Workflow or API
    participant E as Models Entry
    participant S as Solver Entry
    participant P as ProblemSpec Builder
    participant SS as Seed and Policy
    participant R as Root Engine
    participant C as Conditions and Residual
    participant G as Gap Solver
    participant M as Physics Model
    participant GV as Governance
    participant O as SolverResult

    W->>E: call solve* at point theta
    E->>S: solve or solve_constraint
    S->>P: build problem contract by mode or spec
    S->>SS: resolve seed strategy and root policy

    loop attempt plan
        R->>C: request residual function
        C->>G: evaluate gap residual or constraint residual
        G->>M: query omega thermo rho mass
        M-->>G: return physics quantities
        G-->>C: return residual values
        C-->>R: return F(x)
        R-->>S: candidate with diagnostics
    end

    S->>GV: evaluate hard rules and select candidate
    GV-->>S: selected candidate and failure reasons
    S->>O: normalize payload to SolverResult
    O-->>E: return standard result
    E-->>W: deliver equilibrium state and thermo outputs
```

## 2. 分层责任映射（与时序对应）

- `Workflow/API`：给定控制参数，发起点求解或扫描。
- `Models Entry`：统一外部调用入口，屏蔽内部实现差异。
- `Solver Entry`：参数归一化、调度链路组织、结果标准化。
- `ProblemSpec`：构造约束问题契约（当前为 mode 主导，后续可 spec-first）。
- `Seed and Policy`：设置初值策略、容差与 fallback 规则。
- `Root Engine`：执行主方法与回退策略，产出尝试诊断。
- `Conditions and Residual`：将物理约束转译为可求解残差。
- `Gap Solver + Physics Model`：完成具体物理量计算与方程评估。
- `Governance`：做硬约束过滤、候选择优和失败归因。
- `SolverResult`：输出统一结构供上层 workflow 复用。

## 3. 可用于评审的简化口径

一句话：

- Solver 在项目中的作用是把“物理约束定义”转换成“可治理的数值求解流程”，并稳定交付统一格式的平衡态结果给上层业务链路。

## 4. 异常分支时序（主方法失败路径）

```mermaid
sequenceDiagram
    autonumber
    participant W as Workflow or API
    participant S as Solver Entry
    participant PS as Primary Strategy
    participant R as Root Engine
    participant GV as Candidate Governance
    participant LEG as Legacy Adapter Plugin
    participant O as SolverResult

    W->>S: solve request
    S->>PS: build strategy bundle
    PS->>R: run primary method
    R-->>PS: failed or degraded candidate

    alt strategy allows fallback
        PS->>R: run fallback method or extra seeds
        R-->>PS: fallback candidates
    else no fallback
        PS-->>S: primary-only candidate
    end

    PS-->>S: candidate set with diagnostics
    S->>GV: hard rules and selector

    alt governance rejected and plugin enabled
        GV-->>S: rejected with reasons
        S->>LEG: optional rescue solve
        LEG-->>S: rescue candidate or failure
        S->>GV: re-evaluate rescue candidate
    else governance accepted
        GV-->>S: selected candidate
    end

    S->>O: finalize result and diagnostics
    O-->>W: return
```
