# ρ 自适应加密模块 `PNJL.AdaptiveRhoRefinement`

`src/pnjl/scans/AdaptiveRhoRefinement.jl` 提供纯数据层的加密建议算法，可在已有 `TrhoScan` 结果上定位 `|Δμ/Δρ|` 接近 0 的区段，并生成需要补样的 ρ 点。该模块不直接调用求解器，而是返回一组新增 ρ，供脚本或上层逻辑再次触发 `run_trho_scan`。

## 导出符号

| 符号 | 说明 |
| --- | --- |
| `AdaptiveRhoConfig` | 配置容器，控制斜率阈值、最小间距、最大加密点数等 |
| `suggest_refinement_points(rho, mu; config)` | 根据现有 `(ρ, μ)` 样本返回新增 ρ 列表 |
| `merge_rho_values(existing, additions; digits)` | 将新增 ρ 与原网格合并并按精度去重 |

## `AdaptiveRhoConfig`

```julia
AdaptiveRhoConfig(; slope_tol=5.0, min_gap=0.002, max_points=64, digits=6)
```

- `slope_tol`：当区段斜率 `|Δμ/Δρ|` 小于该阈值（单位 MeV）即视为“平坦”，需要加密；默认 5 MeV。
- `min_gap`：仅对 ρ 间隔大于此值的区段提出加密建议，避免重复采样已经很密的区域。
- `max_points`：限制每条曲线的最大新增样本数。
- `digits`：返回值的四舍五入精度，同时用于 `merge_rho_values` 去重。

## 建议流程

1. 使用 `load_trho_curves`（或其它 CSV 解析逻辑）获取某温度下的 `(ρ, μ)` 样本；
2. 调用 `suggest_refinement_points` 获得新增 ρ；
3. 将这些 ρ 传入 `PNJL.TrhoScan.run_trho_scan`，并保持 `resume=true` 以避免重复写入；
4. 若需要合并新的 ρ 网格，可通过 `merge_rho_values` 生成统一的、排序的向量。

## 示例

```julia
using PNJL.TrhoScan
using PNJL.AdaptiveRhoRefinement

curve = [(0.10, 325.0), (0.15, 323.0), (0.25, 323.2), (0.40, 335.0)]
rho_vals = [p[1] for p in curve]
mu_vals = [p[2] for p in curve]
config = AdaptiveRhoConfig(; slope_tol=3.0, min_gap=0.01)
extra_rho = suggest_refinement_points(rho_vals, mu_vals; config=config)
# => e.g. [0.125, 0.20]

TrhoScan.run_trho_scan(
    T_values = [131.0],
    rho_values = extra_rho,
    xi_values = [0.0],
    resume = true,
    overwrite = false,
    output_path = "data/outputs/results/pnjl/trho_scan.csv",
)
```

这样可以在检测到 S 形时自动增加低斜率区的数据密度，改善 Maxwell 和 CEP 搜索的稳定性。