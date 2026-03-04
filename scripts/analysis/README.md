# scripts/analysis/

探索性分析与诊断脚本，从 `tests/analysis/` 迁入（2026-03-04）。

## 说明

这些脚本用于数值探索、收敛分析、缓存调试等场景。它们**不纳入 CI**，仅供开发者手动运行。

## 目录结构

- `cache/` — 缓存命中率/开销分析
- `convergence/` — 数值收敛性测试（节点数、角度、sigma 网格）
- `pnjl/` — PNJL 求解器参数扫描分析
- `relaxtime/` — 输运计算数值探索
- `relaxtime_diagnostics/` — 积分诊断工具
- `relaxtime_validation/` — 非正式验证脚本（宽松断言）
- `struct_migration/` — 迁移阶段文档归档
- `unit_summaries/` — 单元测试摘要文档

## 运行方式

```bash
julia --project=. scripts/analysis/relaxtime/analyze_1_over_p.jl
julia --project=. scripts/analysis/convergence/run_convergence.jl
```
