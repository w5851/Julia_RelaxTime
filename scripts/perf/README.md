# scripts/perf

本目录存放性能观测与瓶颈诊断脚本，不承担 correctness 测试职责。

## 目录分工

- `scripts/perf/pnjl/`
  - PNJL / `Models` 主链相关的性能探针与基线脚本
- `scripts/perf/relaxtime/`
  - relaxtime / scattering / transport 相关的性能探针
  - `bulk_derivative_context_probe.jl` 对比旧式五组独立 Taylor series 与 workflow equilibrium 锁支、三方向共享线性化路径；只报告测量结果，不设置 correctness 或固定加速比门槛
  - `bench_meson_conserved_charge_partial_feedback.jl` 在完成一次 warm-up 后，以新 evaluator/cache 重复测量 charged meson partial-feedback 单点，并输出 baseline、gap、BU、outer 与总耗时；不承担 correctness 或 production gate
  - 示例：`julia --project=. scripts/perf/relaxtime/bulk_derivative_context_probe.jl --repeats=3`

## 稳定入口性能基线

P1-F 起，稳定 CLI 入口的最小基线脚本为：

- `scripts/perf/pnjl/stable_entrypoint_baseline.jl`
- `scripts/perf/pnjl/compile_time_breakdown.jl`
- `scripts/perf/pnjl_phase_thermal_quadrature_probe.jl`

用途：

- 跑真实稳定入口命令，而不是只测内部函数
- 记录端到端 wall-clock 时间
- 固定输出 JSON + Markdown 报告，便于后续复测比较

当前默认覆盖：

- `scripts/models/run_unified_scan.jl scan tmu`
- `scripts/models/run_unified_scan.jl scan trho`

其中：

- `stable_entrypoint_baseline.jl` 记录冷启动端到端 wall-clock 基线
- `compile_time_breakdown.jl` 对比冷启动 CLI、同进程热态调用、以及可选 sysimage CLI，用于判断编译/JIT 占比
- `pnjl_phase_thermal_quadrature_probe.jl` 在多温度和 `xi` 固定态上比较 `24/8`、`32/10` tensor 与 RS-reduced adaptive 热势的耗时、分配和数值；它只提供性能证据，不替代 correctness/phase 收敛 gate

示例：

```powershell
julia --project=. scripts/perf/pnjl/stable_entrypoint_baseline.jl --samples=1 --warmups=0
```

输出目录：

- `data/outputs/perf/stable_entrypoints/<tag>/baseline_summary.json`
- `data/outputs/perf/stable_entrypoints/<tag>/baseline_summary.md`
- `data/outputs/perf/compile_breakdown/<tag>/compile_breakdown.json`
- `data/outputs/perf/compile_breakdown/<tag>/compile_breakdown.md`

版本库策略：

- `data/outputs/perf/` 下的运行产物默认不纳入版本库。
- 进仓内容以脚本、模板、说明文档为主。
- 若某次测量需要作为里程碑证据保留，优先人工提炼成 `docs/` 或开发任务单中的摘要，而不是直接提交整批原始日志与 CSV。

## 最小报告模板

后续新增稳定入口基线时，报告至少保留这些字段：

1. 入口标识：`key` / `label`
2. 真实 CLI 命令
3. 工作负载假设与输出文件名
4. 机器假设：OS、CPU threads、Julia 版本
5. 样本配置：`samples` / `warmups`
6. 工作负载规模：例如 `len(T)*len(mu)*len(xi)` 或 `len(T)*len(rho)*len(xi)`
7. 汇总指标：总耗时 `min` / `median` / `mean` / `max`
8. 单点平均耗时：至少给出 `mean ms/pt`
9. 原始样本日志位置与输出文件存在性
