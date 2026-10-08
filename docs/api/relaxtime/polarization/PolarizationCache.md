# 极化函数缓存

`PolarizationCache` 为 `PolarizationAniso.polarization_aniso` 的结果提供进程内缓存，供传播子与输运链复用。普通计算优先使用 `Models` 工作流；直接调用缓存适用于已有同口径夸克参数和预计算 A 函数的调用方。

## 接口与单位

```julia
polarization_aniso_cached(channel::Symbol, k0::Float64, k_norm::Float64,
    m1::Float64, m2::Float64, μ1::Float64, μ2::Float64,
    T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64,
    A1_value::Float64, A2_value::Float64, num_s_quark::Int)
```

| 参数 | 含义与单位 |
| --- | --- |
| `channel` | `:P` 赝标量或 `:S` 标量；其他值抛出 `ArgumentError` |
| `k0`, `k_norm` | 介子能量、三动量模，`fm⁻¹` |
| `m1`, `m2`, `μ1`, `μ2`, `T` | 夸克质量、化学势和温度，`fm⁻¹`；温度为正 |
| `Φ`, `Φbar`, `ξ` | 无量纲；沿用下游极化函数的物理约定与各向异性适用范围 |
| `A1_value`, `A2_value` | 对应两味夸克的 A 函数，`fm⁻²` |
| `num_s_quark` | 下游兼容参数；`1` 启用含奇异味道通道的能量对称化，`0` 不启用 |

返回 `Tuple{Float64,Float64}`，顺序为 `(Π_real, Π_imag)`，两项单位均为 `fm⁻²`。需要复数时由调用方构造 `complex(Π_real, Π_imag)`。该方法不做 MeV 转换，不接受 AD 数作为缓存参数；模型背景、A 函数与极化函数必须使用同一套单位和处方。

其余两个公开函数：

| 函数 | 行为 |
| --- | --- |
| `reset_cache!()` | 清空条目和计数器，返回 `nothing` |
| `get_cache_stats()` | 返回 `total_calls`、`cache_hits`、`cache_misses`、`hit_rate`、`cache_size` 的 `NamedTuple`；无调用时命中率为 `0.0` |

## 最小示例

在仓库根目录运行：

```julia
include("src/models/Models.jl")
using FastGaussQuadrature: gausslegendre
cache = Main.RelaxTime.PolarizationCache
integrals = Main.RelaxTime.OneLoopIntegrals
hbarc = Main.Constants_PNJL.ħc_MeV_fm

m = 300.0 / hbarc
T = 150.0 / hbarc
u, w = gausslegendre(16)
nodes, weights = 10.0 .* (u .+ 1.0), 10.0 .* w
A_value = integrals.A(m, 0.0, T, 0.5, 0.5, nodes, weights)
args = (:P, 1.0, 0.1, m, m, 0.0, 0.0, T, 0.5, 0.5, 0.0,
        A_value, A_value, 0)
cache.reset_cache!()
first_result = cache.polarization_aniso_cached(args...)
second_result = cache.polarization_aniso_cached(args...)
@assert first_result == second_result
@assert cache.get_cache_stats().cache_hits == 1
cache.reset_cache!()
```

示例节点用于验证缓存调用；实际计算的节点与热尾范围沿用所属工作流的收敛设置。

## 缓存语义与边界

- 缓存键包含全部输入参数。Float64 键保留 40 位 mantissa，量化尺度约为 `2⁻⁴⁰`；同桶复用首先写入的计算结果。相近参数也可能跨桶，不能把它理解为任意两值的 `isapprox` 比较。
- 底层计算使用原始参数，键量化不代表极化函数误差或奇点附近误差已有同样的上界。
- 全局字典和计数器没有锁，不能供多个线程并发读写。使用独立进程时，各进程分别加载、统计和清理缓存；当前没有每线程缓存实例接口。
- 条目没有自动淘汰。可在任务结束时清理以控制内存；所有参数都进入键，改变温度本身不要求为正确性清空缓存。
- 缓存仅减少重复求值，收益取决于实际命中率与工作负载；本页不承诺固定耗时或加速倍数。

## 实现与验证

- [缓存实现](../../../../src/relaxtime/PolarizationCache.jl)
- [下游极化函数](../../../../src/relaxtime/PolarizationAniso.jl)
- [极化函数公式](../../../reference/formula/relaxtime/polarization/Polarization_极化函数byB0.md)
- [缓存行为测试](../../../../tests/unit/relaxtime/test_polarization_cache.jl)

```sh
julia --project=. tests/unit/relaxtime/test_polarization_cache.jl
```
