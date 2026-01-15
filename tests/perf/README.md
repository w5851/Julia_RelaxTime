# tests/perf

这里存放**性能/基准/剖析(profile)**相关脚本，用于定位热点、比较策略与回归性能；默认不建议作为 CI 的常规单元测试运行。

## 目录约定

- `tests/perf/pnjl/`：PNJL 求解/扫描相关性能脚本
- `tests/perf/relaxtime/`：弛豫时间与散射相关性能脚本
- `tests/perf/results/`：可选的结果落盘目录（JSON/Markdown 等）

## 文件命名约定（建议）

- `benchmark_*.jl`：BenchmarkTools 基准脚本（可重复运行、可落盘）
- `profile_*.jl`：Profile/统计剖析脚本（用于看时间分布/调用栈）
- `timing_*.jl`：手写计时循环（time_ns/time），用于快速定位“哪一步慢”
- `test_*_performance.jl`：带 `@testset` 的性能烟囱测试（通常 `samples=1`），用于快速回归“是否明显变慢/是否能跑通”

如果文件名不足以自解释，务必在文件顶部加入 docstring，至少包含：

- 对应系统流程的哪一步（例如：总截面 → 散射振幅 → 单圈积分/根查找/节点生成）
- 关键入口函数/模块（指向 `src/...`）
- 如何运行（推荐命令）与输出（是否生成 `*.md`/`*.json`）

## 运行方式

在仓库根目录执行（Windows/PowerShell 亦可）：

- `julia --project=. tests/perf/relaxtime/benchmark_*.jl`
- `julia --project=. tests/perf/relaxtime/profile_*.jl`
- `julia --project=. tests/perf/relaxtime/timing_*.jl`

注：这些脚本可能会做预热、运行多次迭代，耗时与机器/线程数有关。
