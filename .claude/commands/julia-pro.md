---
description: Julia 科学计算/性能优化/多重派发设计指导
---

# Julia Pro

## 核心流程
1. 明确数值目标、不变量、预期精度/容差。
2. 在边界处显式约定数据契约（输入/输出类型与单位）。
3. 保持热路径 kernel 类型稳定、分配轻量。
4. 先验证正确性，再优化性能。
5. Profile 热路径，只优化实测瓶颈。

## 实现规则
- 热循环中优先使用具体容器元素类型。
- 计算 kernel 中避免全局可变状态。
- IO/配置解析与数值 kernel 隔离。
- 遇到抽象输入时使用函数屏障（function barrier）。
- 避免类型盗用（type piracy），在自己的命名空间中扩展方法。
- API 保持可组合、可预测。

## 性能检查清单
- 对代表性工作负载运行分配检查（`@allocated`）。
- 确认关键 kernel 返回类型稳定（`@code_warntype`）。
- 每次优化前后做 benchmark 对比。
- 用确定性 smoke case 防止回退。

## 本仓库约定
- 内部物理量使用自然单位制；fm^-1 用 `_inv_fm` 后缀。
- 外部 MeV 量用显式命名（`T_MeV`、`mu_MeV`）。
- 高阶单位遵循本地惯例（`sigma_fm4`、`coupling_inv_fm4`）。
- `using`/`import` 置顶；优先显式导入；相对导入正常使用（`using ..GaussLegendre`）。
- 保留 `if !isdefined(Main, :ModuleName)` 守卫。
- 用多重派发而非大标志函数扩展行为。
- `const` 声明全局绑定；`@inline`/`@inbounds` 仅热路径使用。
