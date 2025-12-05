````markdown
# MaxwellRhoMu 模块（ρ–μ Maxwell 等面积）

`src/pnjl/analysis/MaxwellRhoMu.jl` 提供基于 `TrhoScan` 数据的 Maxwell 等面积实现，用来从 ρ(μ) 曲线求解一阶相变的共存化学势与两相密度。

## 导出符号

| 名称 | 类型 | 说明 |
|------|------|------|
| `MaxwellResult` | 结构体 | 记录单个温度的等面积求解结果（共存 μ、ρ 以及残差、迭代信息） |
| `maxwell_rho_mu` | 函数 | 对某条 ρ(μ) 曲线执行 ρ–μ 等面积求解 |
| `build_phase_boundary` | 函数 | 批量地对 `Dict{T => (μ, ρ)}` 数据求解 Maxwell，形成相界结果映射 |
| `phase_boundary_from_rows` | 函数 | 直接接受 `CSV.File` 等行迭代器，自动按温度分组后执行 Maxwell |

## `MaxwellResult`

字段：
- `converged::Bool`：是否找到满足等面积条件的解；
- `mu_coex_MeV::Union{Nothing, Float64}`：共存化学势（MeV）；
- `rho_gas/rho_liquid::Union{Nothing, Float64}`：等面积水平线与 ρ(μ) 曲线左右交点，分别对应气相/液相密度；
- `area_residual::Union{Nothing, Float64}`：最终积分残差 `|∫(μ - μ_coex) dρ|`；
- `iterations::Int`：二分迭代次数（若直接命中解则为 0）；
- `details::Dict{Symbol, Any}`：调试信息（包含 μ 搜索区间、ρ 截断范围、spinodal 初始提示等）。

`MaxwellResult()` 提供一个“未收敛”默认构造。

## `maxwell_rho_mu(mu_vals, rho_vals; kwargs...)`

- 输入为同一温度下的 `(μ, ρ)` 采样（可以是 `Vector{Float64}`、`Vector{Real}` 等）；
- 内部会过滤掉 `NaN`、`Inf`，并按照 ρ 升序构建曲线；
- 通过 `detect_s_shape` 自动获取 spinodal 作为 μ 搜索区间，随后用候选采样 + 二分法寻找 `μ_coex` 使等面积积分为 0；
- 默认关键字：
  - `min_samples=12`：曲线点数不足时直接返回失败；
  - `detect_min_points=6`、`detect_eps=1e-6`：传入 S 形检测的参数；
  - `candidate_steps=64`：初始扫描 μ 值的数量，用于寻找正负残差；
  - `max_iter=60`、`tol_area=1e-4`：二分终止条件；
  - `spinodal_hint`：可选地传入 `SShapeResult`，避免重复检测。

若成功，返回的 `MaxwellResult` 带有气/液相密度与残差；失败时 `details[:reason]` 会指示具体原因（如 `"no_s_shape"`, `"insufficient_points"` 等）。

## `build_phase_boundary(curves; kwargs...)`

- `curves` 需为 `Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}`，可直接由 `CEPFinder.build_curves` 生成；
- 依温度升序遍历，为每条曲线调用 `maxwell_rho_mu`，返回 `Dict{T => MaxwellResult}`；
- 通过 `kwargs...` 将 `min_samples`、`tol_area` 等参数层层传递，可在批量求解时统一调参。

## `phase_boundary_from_rows(rows; xi=0.0, tol=1e-6, kwargs...)`

- 方便从 `CSV.File`、`Tables.rowtable` 等行式数据直接获得相界；
- 内部调用 `group_curves_by_temperature(rows; xi, tol)` 过滤出指定 ξ 的样本，再用 `build_curves` 与 `build_phase_boundary` 求解；
- 典型输入为 `TrhoScan` 生成的 CSV：

```julia
using CSV, DataFrames, PNJL
rows = CSV.File("data/outputs/results/pnjl/trho_scan.csv")
results = PNJL.MaxwellRhoMu.phase_boundary_from_rows(rows; xi=0.0)
for (T, res) in sort(collect(results))
    res.converged || continue
    @info "coexistence" T res.mu_coex_MeV res.rho_gas res.rho_liquid
end
```

## 工作流建议

1. 使用 `TrhoScan.run_trho_scan` 生成 T-ρ 网格数据并写入 CSV；
2. 通过 `CSV.File` 读取后交给 `group_curves_by_temperature` 或直接使用 `phase_boundary_from_rows`；
3. 成功求得的一阶相变点可与 `CEPFinder.find_cep` 的结果结合，用于绘制 T-μ 相图或提取气/液相线。
````