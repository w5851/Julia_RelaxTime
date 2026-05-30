# RelaxationTime 模块 API

弛豫时间计算，基于平均散射率表达式。

## 导出符号
- `REQUIRED_PROCESSES`：需要的散射过程符号列表
- `compute_average_rates(quark_params, thermo_params, K_coeffs; ...)`
- `relaxation_rates(densities, rates)`
- `relaxation_times(quark_params, thermo_params, K_coeffs; densities, ...)`

## 兼容别名与退场节奏

本模块仍保留少量历史兼容入口，便于旧脚本平滑迁移；新代码应统一使用 canonical API。

### 当前兼容别名清单（P1-R4）

- `rate_lookup(rates, key)` -> **替代接口**：`rate_value(rates, key)`
  - 现状：调用时触发 `Base.depwarn`
  - 用途：仅历史调用兼容

- `rate_value` 的历史过程键别名（通过 `RATE_ALIASES` 解析）
  - 示例：`:dubar_to_dubar -> :udbar_to_udbar`
  - 示例：`:subar_to_subar -> :usbar_to_usbar`
  - 示例：`:ubardbar_to_ubardbar -> :ud_to_ud`
  - 示例：`:ubarubar_to_ubarubar -> :uu_to_uu`
  - 示例：`:ubarsbar_to_ubarsbar -> :us_to_us`
  - 示例：`:sbarsbar_to_sbarsbar -> :ss_to_ss`
  - 现状：兼容接受，不主动告警
  - 用途：旧数据/旧脚本过程名兼容

### 退场计划（文档口径）

- **阶段 A（当前）**：保留兼容入口；新文档与新示例仅展示 canonical 名称。
- **阶段 B（不早于 2026-06-30）**：
  - `rate_lookup` 继续保留但强化弃用提示（迁移提示指向 `rate_value`）。
  - 过程键旧别名继续兼容，但在文档中标记为“legacy only”。
- **阶段 C（移除门槛达成后）**：移除 `rate_lookup`。

### 移除门槛

满足以下条件后，才允许进入阶段 C：

1. 兼容测试替换完成（`test_relaxation_time_compat.jl` 不再依赖 `rate_lookup`）。
2. 文档与脚本示例不再出现旧入口。
3. 至少一个完整发布周期内无新的兼容回归反馈。

## 数学定义
$$\tau_i^{-1} = \sum_j \rho_j \; \bar{w}_{ij}$$

## 物理公式

详见 `docs/reference/formula/relaxtime/transport/RelaxationTime_FromAverageScatteringRate.md`。

## 默认行为

当不传入 `p_grid`/`p_w` 时：
- 动量积分范围：`[0, Λ]`（Λ ≈ 3.05 fm⁻¹）
- σ(s) 缓存范围：基于 Λ 截断的 w0cdf 设计
- 数密度积分：半无穷 `[0, ∞)`（物理要求，由调用方提供）

这确保了动量积分范围与 σ(s) 缓存范围的一致性。

`AverageScatteringRate.md` 中的默认截断与缓存口径应与本节保持一致；后者更偏向单过程积分与缓存细节，本页负责 τ 聚合层的默认行为说明。

## 参数说明

- `quark_params`：夸克参数 `(m, μ, A)`
- `thermo_params`：热力学参数 `(T, Φ, Φbar, ξ)`
- `K_coeffs`：有效耦合系数
- `densities`：数密度（**必须使用半无穷积分计算**）
- `p_grid/p_w`：可选，自定义动量积分节点
- `sigma_cutoff`：σ(s) 有效范围的动量截断（默认 Λ）

## 本次提交的变更（e9cb5fbd...，2026-01-29）

- `compute_average_rates(...)` 中对 `cs_cache` 的处理逻辑有调整：
    - 若调用方已在 `cs_caches[process]` 中放置了已填充的 `CrossSectionCache`，`compute_average_rates` 将复用该缓存；否则将把 `cs_cache` 参数以 `nothing` 传递给 `average_scattering_rate`，让后者按需在内部构建缓存并应用阈值减法策略（若请求）。
    - 这通过内部变量 `cs_cache_arg` 实现（若缓存缺失或为空则传入 `nothing`）。

- 推荐调用模式：
    1. 如果需复用截面，调用者应先用 `build_w0cdf_pchip_cache(...)` 或 `precompute_cross_section!` 构建带 `fingerprint` 的 `CrossSectionCache`，放入 `cs_caches[process]`，以避免重复构建开销。
    2. 对于一次性评估，可传 `cs_cache=nothing` 或不传该参数，让内部按需构建（会使用 `threshold_subtraction` 等参数）。
    3. 多个参数点之间只允许在同一生成指纹下复用截面；`quark_params`、`thermo_params`、`K_coeffs`、`n_sigma_points`、`threshold_subtraction` 或 w0cdf 网格策略不一致时，已填充缓存会被拒绝并抛出 `ArgumentError`。

- 兼容策略：
    - 旧的手工 `CrossSectionCache(process, s_vals, sigma_vals)` 和 `load_cross_section_caches_from_dir` 读取的两列 CSV 表没有 `fingerprint`；默认允许复用但会 warning。
    - 若调用方需要严格治理，可传 `require_cache_fingerprint=true`，此时无指纹缓存会被拒绝。

- 关键影响：性能（避免重复构建）和语义（当启用阈值减法时，缓存保存的是 `raw - asym`，返回值会把解析项加回）。


## 典型用法

```julia
using .RelaxationTime
using .GaussLegendre: gauleg
using .Constants_PNJL: Λ_inv_fm

# 方式1：使用默认行为（推荐）
result = relaxation_times(quark_params, thermo_params, K_coeffs;
    densities=densities)

# 方式2：有限截断积分（与 C++ 一致）
p_grid, p_w = gauleg(0.0, 15.0, 64)
result = relaxation_times(quark_params, thermo_params, K_coeffs;
    densities=densities,
    p_grid=p_grid, p_w=p_w,
    sigma_cutoff=Λ_inv_fm)

println(result.tau)      # 各味的 τ
println(result.tau_inv)  # 各味的 τ^-1
```

## 与 C++ 实现对比

| 参数 | Julia 默认 | C++ 实现 | 匹配方式 |
|------|-----------|----------|----------|
| 动量积分范围 | [0, Λ] | [0, 15] fm⁻¹ | 传入 `p_grid=gauleg(0,15,n)` |
| σ(s) 缓存范围 | 基于 Λ 截断 | 基于 Λ 截断 | 已一致 |
| 数密度积分 | [0, ∞) | [0, ∞) | 已一致 |

使用 `p_grid=gauleg(0,15,n)` + `sigma_cutoff=Λ_inv_fm` 后，与 C++ 结果误差 < 0.3%。

若你正在排查单个过程的 σ(s) 缓存、阈值渐近项与平均散射率积分细节，请继续阅读 `AverageScatteringRate.md`，本页不重复展开单过程缓存实现。

## 物理说明

- **数密度积分范围**：必须使用半无穷积分 [0, ∞)，这是物理要求
- **σ(s) 有效范围**：由 PNJL 模型的动量截断 Λ 决定
- **重夸克注意**：当夸克质量 m 接近 Λ 时（如低温下的 s 夸克），相空间被限制
