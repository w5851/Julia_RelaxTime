# benchmark/

本目录保存性能基准与隔离的外部数值 oracle。正确性测试仍位于 `tests/`；`scripts/perf/` 用于聚焦 profiling，不承担回归门禁。

## 环境边界

- 根 `Project.toml` 是 production、稳定 CLI 与常规测试环境。
- `benchmark/Project.toml` 只保存 benchmark 专用依赖；QuadGK 仅在这里声明。
- 首次使用先运行：

```sh
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
```

本仓库是 include-driven，部分基准同时需要根环境依赖。此时应显式叠加环境，而不是把 benchmark-only 包重新加入根项目。

Windows / PowerShell：

```powershell
$env:JULIA_LOAD_PATH = "@;benchmark;@stdlib"
julia --project=. benchmark/relaxtime/benchmark_quadgk_oracle_smoke.jl
```

Linux / macOS：

```sh
JULIA_LOAD_PATH='@:benchmark:@stdlib' julia --project=. benchmark/relaxtime/benchmark_quadgk_oracle_smoke.jl
```

环境变量只应作用于该 benchmark 进程/终端会话；普通开发、测试和 production 不叠加 `benchmark/`。

## 目录与命名

- `bench_<feature>.jl` / `benchmark_<feature>.jl`：可复现性能或 oracle 探针。
- `benchmark/benchmarks.jl`：收集符合 `bench_*.jl` 命名的 BenchmarkTools suite；它不是当前仓库的 PkgBenchmark package entrypoint。
- `benchmark/pnjl/`：PNJL 求解与热力学基准。
- `benchmark/relaxtime/`：散射、传播子、数密度与输运基准。

benchmark 输出不得自动晋升为 production 正确性证据。数值 claim 仍需独立节点/容差收敛、对应 validation/regression 和 provenance。

## charged GBU 相图单点成本

`relaxtime/bench_charged_gbu_contour_point.jl` 使用根环境，不额外引入 oracle。
默认同进程背景重复 5 次、screening 密度重复 3 次、production 密度重复 1 次。`BENCH_CGBU_POINTS` 使用
`T_MeV:muB_MeV` 列表；`BENCH_CGBU_WORKLOAD=screening|production` 分别测当前
无限热核的低分辨率探针和完整 smoke 门禁。production 始终读取并验证原配置，
不降低门槛。失败样本保留原因；成功产点的成本只使用门禁通过的样本。

首调用和预热后的重复样本单独保存，compile/GC/bytes 不从耗时中机械相减；
`@timed` 的编译计数不保证包含调用者预先 inference。进程启动、包加载和 include
在内部样本之外。远程 workflow 另外保存进程 wall time，源码快照、Project/Manifest
及前后 hash，检查是否在运行中变更。纯合成测试不使用速度阈值。

GitHub Actions 的 `Charged GBU contour benchmark` 仅测代表点；初次在
`codex/charged-gbu-contour-benchmark` 的窄路径 push 上验证，随后可手动重复。
二维网格不由此 workflow 启动。网格步长需根据 remote 样本及失败区域另行确定。
