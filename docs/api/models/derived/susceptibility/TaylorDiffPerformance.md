# TaylorDiff 性能参考

本页记录 PNJL `chi_B` TaylorDiff backend 的代表性性能探针结果，用作后续优化和回归判断的静态参考。性能数字不是 correctness gate；需要复测时使用 `scripts/perf/pnjl_chi_b_taylordiff_probe.jl`。

## 探针命令

```sh
julia --project=. scripts/perf/pnjl_chi_b_taylordiff_probe.jl --max_order=10 --p_num=8 --t_num=4 --repeats=7
```

## 参考点

- 日期：2026-05-24
- 点位：PNJL xi=0 真 CEP
- `T_CEP_MeV=130.3369140625`
- `muB_CEP_MeV=877.1423865705508`
- 积分规模：`p_num=8`, `t_num=4`
- TaylorDiff 阶数：`max_order=10`
- `linear_solve=:auto`

## 结果快照

| 指标 | 数值 |
| --- | ---: |
| cold wall | 90.138 s |
| first-call | 41.358 s |
| first alloc | 3.603 GB |
| same-process next point | 0.0157 s |
| warm median | 0.0155 s |
| warm min | 0.0153 s |
| warm alloc median | 122 KB |
| peak RSS | 889.2 MB |

一次调用返回 `chi_B` 的 1 到 10 阶：

```text
[1.026550339262278, 18.64931110677032, -326.15207378607215, -652726.4239053279, 9.42551324151821e7, 2.2246863492661136e11, -8.66543098872892e13, -2.0788273590737366e17, 1.604389492911533e20, 3.726625095524057e23]
```

## 解读

- cold wall 和 first-call 主要反映 Julia/TaylorDiff/PNJL AD 专门化编译成本。
- 同一进程扫下一点接近 warm runtime，说明不需要为相邻点重新支付 first-call JIT。
- warm runtime 保持在几十毫秒以下，满足当前 high-order 单方向 `chi_B` 路径的目标。
- `linear_solve=:auto` 当前低阶到 `order=10` 使用 `:refactor_each_order`，`order >= 16` 才切到 `:factorized_each_order`；这是基于旧探针中 order 4/10 未稳定受益、order 16 明显受益的策略。

## Mixed BQS Jet Probe

高阶 mixed BQS 走内部 multivariate Taylor jet，不替代单方向 TaylorDiff fast path。代表性探针命令：

```sh
julia --project=. scripts/perf/pnjl_chi_bqs_mixedjet_probe.jl --p_num=8 --t_num=4 --repeats=5
```

参考点：

- 日期：2026-05-25
- 点位：`T_fm=0.57`, `muB_fm=0.18`, `muQ_fm=0.05`, `muS_fm=0.02`
- 同进程下一点：`T_fm + 0.0025`
- 积分规模：`p_num=8`, `t_num=4`
- `linear_solve=:auto`
- 整体 cold wall：265.702 s

| orders | active axes | first-call | first alloc | next point | warm median | warm alloc median | peak RSS | value |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `(1,1,0)` | 2 | 45.746 s | 1.825 GB | 0.0377 s | 0.0242 s | 74.4 KB | 931.9 MB | 6.083433636767924e-5 |
| `(2,1,0)` | 2 | 32.130 s | 296.6 MB | 0.0721 s | 0.0743 s | 86.4 KB | 931.9 MB | 3.57905081599011e-5 |
| `(1,1,1)` | 3 | 41.218 s | 472.1 MB | 0.0461 s | 0.0353 s | 24.75 MB | 976.2 MB | -4.1028573456403393e-7 |
| `(2,1,1)` | 3 | 23.135 s | 721.9 MB | 0.0985 s | 0.1086 s | 50.57 MB | 1.043 GB | 4.758676296828763e-6 |

解读：

- mixed jet 的 cold/first-call 仍主要是 Julia 专门化成本；同进程下一点接近 warm runtime。
- `D=2` mixed 路径热分配维持在 KB 级；`D=3` 当前仍有 10MB 到 50MB 级热分配，是后续 mixed jet 优化的主要观察点。
- 单方向 `B/Q/S` 不走该通用 mixed jet 路径，因此这些 mixed 常数项不会影响已验证的单变量 TaylorDiff `chi_B` fast path。
