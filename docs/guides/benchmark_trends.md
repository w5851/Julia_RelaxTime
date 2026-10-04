# PNJL benchmark 趋势观测

覆盖现有单点求解与 T–μ、T–ρ 扫描，记录耗时、分配量和扫描收敛率。它用于比较性能；数值正确性仍由 unit、integration、regression 和 validation 负责。

## 运行与查询

[benchmarks workflow](https://github.com/w5851/Julia_RelaxTime/actions/workflows/pnjl-benchmarks.yml) 保留相关 push、PR 和手动触发，并配置为每周一 02:00 UTC（北京时间 10:00）运行。定时触发在该定义进入默认分支后生效。CI 使用 Ubuntu 24.04、Julia 1.12.5、1 个 Julia 线程和 1 个 BLAS 线程。

- [main 的运行记录](https://github.com/w5851/Julia_RelaxTime/actions/workflows/pnjl-benchmarks.yml?query=branch%3Amain) 是最新结果和历史结果的稳定入口。
- 打开运行页即可阅读 Action Summary；下载 `pnjl-benchmark` artifact 获取机器可读记录。成功和失败运行均尝试保留报告，保留期为 **90 天**。
- 首次运行或历史记录已过期时仍记录本次结果；下一次符合条件的运行即可用于比较。

```sh
gh workflow run pnjl-benchmarks.yml --ref main
gh run list --workflow pnjl-benchmarks.yml --branch main --status success --limit 10
gh run download <run-id> --name pnjl-benchmark --dir <destination>
```

## 本地运行

在仓库根目录、已实例化的根 Julia 环境中依次运行。可将输出放在临时目录，避免覆盖已有报告：

```powershell
$env:PNJL_BENCHMARK_OUTPUT_DIR = Join-Path $env:TEMP "pnjl-benchmark"
$env:JULIA_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
julia --project=. benchmark/pnjl/single_point_solver_perf.jl
julia --project=. scripts/perf/pnjl/scan_perf.jl
julia --project=. scripts/dev/check_benchmark_thresholds.jl
# 仅在上面的阈值检查成功后标记 success
python scripts/dev/benchmark_trends.py --threshold-status success
```

任一步失败先处理该失败。报告脚本只读取已有结果，不重跑求解。未设置输出变量时沿用 `tests/perf/results/pnjl/`；这是历史报告位置，不是新增测试脚本目录。Linux/macOS 用 `export PNJL_BENCHMARK_OUTPUT_DIR=...` 等设置相同环境变量。

比较显式保存的本地快照：

```sh
python scripts/dev/benchmark_trends.py --baseline <previous>/benchmark_snapshot.json --threshold-status success
```

可用 `--results-dir` 指向其他结果目录；`--github-baseline --repository w5851/Julia_RelaxTime` 通过已认证的 `gh` 读取远端历史。默认本地命令不联网。

## 比较口径

自动查找最近 **20 次成功的 main 运行**，按新到旧选择首个完整、同配置、未过期的快照。只接受 push、schedule 或 workflow_dispatch；PR、失败阈值、相关源码有未提交改动、来源 commit/run 不匹配的记录不能成为 main 基线。比较年龄上限为 **35 天**，可通过 `--max-baseline-age-days` 显式调整。

同配置包括 Julia/BLAS 线程、Julia 版本、OS/架构、CPU、runner 镜像系列、已解析 Project/Manifest、PNJL/physics 配置及所选 profile、benchmark 脚本以及工作负载参数。runner 镜像补丁版本也会记录，但不单独阻止比较。源代码 commit 可以不同，否则无法观察代码修改后的性能变化。

benchmark 先预热再计时，`evals=1`。报告同时记录请求采样数、BenchmarkTools 时间预算与实际采样数；时间预算可能使实际次数少于请求值。T–ρ 计时包含其现有临时 CSV 写入，扫描收敛率取预热运行结果。

耗时、内存或分配次数增长 **≥25%**，以及扫描收敛率下降 **>2 个百分点**，在 Summary 中标作预警。其他变化也展示具体数值。**相对变化不造成 CI 硬失败**；既有绝对耗时和最低收敛率阈值仍在 `check_benchmark_thresholds.jl` 中执行。

| 比较状态 | 含义 |
| --- | --- |
| `compared` | 已找到同配置基线，显示每项变化 |
| `missing_baseline` | 没有新版快照或尚未建立基线，本次仅记录 |
| `expired_baseline` | 超过比较年龄或 artifact 已过期 |
| `incompatible_baseline` | 环境、配置或指标集合不同，不计算误导性的百分比 |
| `invalid_baseline` | 历史记录缺字段、非有限数值、来源或阈值状态无效 |
| `baseline_unavailable` | GitHub/权限/下载暂不可用，本次报告仍保留 |
| `invalid_current` | 当前报告缺失、环境不一致或数值无效；报告步骤失败 |

## 报告文件

| 文件 | 内容 |
| --- | --- |
| `single_point_benchmark.*`、`scan_benchmark.*` | 原始 JSON 和 Markdown，各自带运行环境及来源 |
| `benchmark_snapshot.json` | 统一指标、单位、工作负载、实际采样数与基线资格 |
| `environment.json` | 可单独查看的环境记录 |
| `comparison.json`、`comparison.md` | 比较状态、选用/拒绝的候选、变化值和预警 |
| `baseline.json` | 本次实际采用的基线副本，便于复算 |

artifact 过期后不承诺永久恢复；需要长期留存的某次观测应在保留期内下载。报告中的性能结论限于所记载的硬件、环境和工作负载。
