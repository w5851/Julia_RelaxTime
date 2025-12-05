# T-μ 扫描模块 `PNJL.TmuScan`

`src/pnjl/scans/TmuScan.jl` 实现了在固定温度/化学势网格上批量调用 `AnisoGapSolver.solve_fixed_mu` 的工具。该模块负责：

- 遍历 (T, μ, ξ) 组合并写入 CSV 结果；
- 自动调用 `SeedCache.find_initial_seed` 以平滑暖启动；
- 支持断点续扫（根据已有 CSV 过滤已完成点）；
- 允许通过关键字参数向底层 NLsolve 传递自定义求解配置。

## 导出符号

| 符号 | 类型 | 说明 |
|------|------|------|
| `run_tmu_scan` | 函数 | 主入口，执行扫描并写入 CSV |
| `DEFAULT_T_VALUES` | `Vector{Float64}` | 默认温度集合（50–200 MeV，步长 10） |
| `DEFAULT_MU_VALUES` | `Vector{Float64}` | 默认化学势集合（0–400 MeV，步长 10） |
| `DEFAULT_OUTPUT_PATH` | `String` | 默认输出 `data/outputs/results/pnjl/tmu_scan.csv` |

## `run_tmu_scan` 接口

```julia
run_tmu_scan(; T_values=DEFAULT_T_VALUES,
                mu_values=DEFAULT_MU_VALUES,
                xi_values=[0.0],
                output_path=DEFAULT_OUTPUT_PATH,
                seed_path=DEFAULT_SEED_PATH,
                overwrite=false,
                resume=true,
                p_num=24,
                t_num=8,
                progress_cb=nothing,
                solver_kwargs...)
```

- `T_values`、`mu_values`、`xi_values`：分别指定温度 (MeV)、平均化学势 (MeV) 与各向异性 ξ 的取值列表。
- `output_path`：结果 CSV 路径；若目录不存在会自动创建。
- `seed_path`：种子表路径，默认使用 `data/raw/pnjl/seeds/sobol_seed_table.csv`。
- `overwrite`：若为 `true` 则重新写入并覆盖旧文件。
- `resume`：若存在旧 CSV，读取其中的 `(T, μ, ξ)` 组合并跳过，支持断点续扫。
- `p_num` / `t_num`：热积分节点数量，直接传递给 `solve_fixed_mu`。
- `progress_cb`：可选回调 `(params, result) -> ...`，每个网格点求解后触发。
- `solver_kwargs...`：透传给 `solve_fixed_mu` 的额外关键字，如 `linesearch = LineSearches.BackTracking()`、`method=:trust_region` 等。

返回值为命名元组 `(total, success, failure, skipped, output)`，用于快速了解扫描统计。

## CSV 字段

`run_tmu_scan` 会写入表头：

```
T_MeV, mu_MeV, xi, pressure_fm4, rho, entropy_fm3, energy_fm4,
phi_u, phi_d, phi_s, Phi1, Phi2, iterations, residual_norm,
converged, message
```

- 热力学量均以 `fm` 自然单位表示；
- `phi_*` 与 `Phi1/Phi2` 对应求解向量；
- `converged` 指示 NLsolve 是否达到公差；
- `message` 记录种子失败或求解异常（为空表示正常）。

## 使用示例

```julia
using PNJL
using PNJL.TmuScan

stats = run_tmu_scan(
    T_values = 80:20:120,
    mu_values = 0:50:200,
    xi_values = [0.0, 0.5],
    output_path = "data/outputs/results/pnjl/tmu_scan_demo.csv",
    p_num = 16,
    t_num = 6,
    solver_kwargs = (; maxiters = 80),
)

@info "scan done" stats
```

运行后可配合 `scripts/pnjl/plot_curves.py` 绘制 `P-μ` 曲线或进一步做 Maxwell 构造。

## 注意事项

1. **种子依赖**：扫描范围应涵盖种子表 `seed_path` 中已有的 `(T, μ, ξ)` 邻域，否则需要扩充种子或提供手动 `seed_state`。
2. **断点续扫**：若启用 `resume`，会逐行读取旧 CSV 并记录完成状态；大文件时建议保持顺序一致以减少额外开销。
3. **性能优化**：
   - 调整 `p_num/t_num` 以在精度与速度间折中；
   - 通过 `progress_cb` 集成实时日志或队列系统；
   - 结合 `Distributed` 或任务队列可并行多个 (T, μ) 区段。
4. **失败处理**：未收敛的点仍会写入 CSV（字段为 `NaN` 或 `false`），便于后续检查或重跑。可用 `message` 字段定位原因。