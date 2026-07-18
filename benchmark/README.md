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
