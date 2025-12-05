# T-ρ 扫描模块 `PNJL.TrhoScan`

`src/pnjl/scans/TrhoScan.jl` 负责在固定温度与归一化密度 (ρ/ρ₀) 网格上批量调用 `AnisoGapSolver.solve_fixed_rho`，并将求解结果写入 CSV，供 ρ-μ 曲线分析、临界点搜索与相图绘制使用。

## 导出符号

| 符号 | 类型 | 说明 |
|------|------|------|
| `run_trho_scan` | 函数 | 主入口，执行 T-ρ 扫描并输出 CSV |
| `DEFAULT_T_VALUES` | `Vector{Float64}` | 默认温度集合（50–200 MeV，步长 10） |
| `DEFAULT_RHO_VALUES` | `Vector{Float64}` | 默认密度集合（0 至 3.0，以 0.05 递增） |
| `DEFAULT_OUTPUT_PATH` | `String` | 默认输出 `data/outputs/results/pnjl/trho_scan.csv` |
| `默认 NLsolve 设置` | `NamedTuple` | `linesearch = LineSearches.BackTracking()`，可用关键字覆盖 |

## `run_trho_scan` 接口

```julia
run_trho_scan(; T_values=DEFAULT_T_VALUES,
                 rho_values=DEFAULT_RHO_VALUES,
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

- `T_values`：温度列表 (MeV)。
- `rho_values`：归一化密度列表（ρ/ρ₀）。默认仅覆盖 0–3 区间，确有需要可显式传入负值测试点。
- `xi_values`：各向异性参数列表。
- `output_path`：结果 CSV 路径；自动创建目录。
- `seed_path`：种子表路径，默认指向 Sobol/LHS 生成的种子。
- `overwrite` / `resume`：与 `TmuScan` 相同，用于重跑或断点续扫。
- `p_num` / `t_num`：高斯–勒让德节点数量，传递给求解器。
- `progress_cb`：可选回调 `(params, result)`。
- `solver_kwargs...`：透传给 `solve_fixed_rho`。默认已启用 `linesearch = LineSearches.BackTracking()`，若需要可显式传入 `linesearch = nothing` 或其他 NLsolve 选项（例如 `method=:trust_region`）。

返回 `(total, success, failure, skipped, output)` 统计信息。

## CSV 字段

输出表头：

```
T_MeV, rho, xi, mu_u_MeV, mu_d_MeV, mu_s_MeV, mu_avg_MeV,
pressure_fm4, entropy_fm3, energy_fm4,
phi_u, phi_d, phi_s, Phi1, Phi2,
iterations, residual_norm, converged, message
```

- `mu_*` 来自 `SolverResult.mu_vec`（已换算为 MeV），`mu_avg_MeV` 为三味平均；
- 其他字段与 T-μ 扫描一致；
- 失败或未收敛的记录会以 `NaN`/`false` 标识，并携带 `message` 描述。

## 使用示例

```julia
using PNJL
using PNJL.TrhoScan
using LineSearches

stats = run_trho_scan(
    T_values = 70:20:110,
    rho_values = 0.2:0.2:1.0,
    xi_values = [0.0],
    output_path = "data/outputs/results/pnjl/trho_scan_demo.csv",
    p_num = 16,
    t_num = 6;
    method = :trust_region,
    # 默认自带 BackTracking，可按需覆盖：
    linesearch = LineSearches.BackTracking(),
)
```

生成的 `ρ-μ` 数据可直接交给 `scripts/pnjl/plot_curves.py` 绘制，也可作为 S 形检测、Maxwell 构造及临界点搜寻的输入。

## 注意事项

1. **多解区域**：在 S 形存在的温度区间，同一 ρ 可能对应多组化学势；为了捕获不同分支，可配合不同初值（修改种子或拆分扫描顺序）。
2. **连续种子**：扫描同一 `(T, ξ)` 时会自动沿 ρ 方向复用上一收敛解，必要时仍可通过 `seed_path`/自定义初值覆盖。
3. **线搜索建议**：固定 ρ 模式更容易出现残差震荡，默认已启用 `LineSearches.BackTracking()`；也可设置 `method=:trust_region` 等其他 NLsolve 选项。
4. **密度范围**：默认密度为 0–3（单位 ρ/ρ₀）。若需探查负密度或更高密度，需显式传入 `rho_values` 并确保种子覆盖。
5. **近似收敛**：若 NLsolve 未报告 `converged=true` 但残差已低于 `1e-4`，扫描会保留该解并在 `message` 中标注 `accepted with residual ...`，以避免舍弃物理可用的点。
6. **自适应加密**：若在 131 MeV 附近需要额外密度点，可先运行基础扫描，再使用 `PNJL.AdaptiveRhoRefinement`（或 `scripts/pnjl/run_adaptive_trho_scan.jl`）读取已有 CSV，自动找出坡度趋零的区段并调用 `run_trho_scan` 补采样。
7. **结果验证**：可结合 `run_solver_smoke.jl` 或局部 `solve_fixed_rho` 调用检查单个点的状态，确保 CSV 数据可信。