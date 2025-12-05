# CEPFinder 模块（临界点检测）

`src/pnjl/analysis/CEPFinder.jl` 基于 `PhaseTransition.detect_s_shape` 的结果，寻找 PNJL 相图中的临界终点 (Critical End Point, CEP)。

## 导出符号

| 名称 | 类型 | 说明 |
|------|------|------|
| `CEPResult` | 结构体 | 存储 CEP 是否存在及其估计值 |
| `build_curves` | 函数 | 将 `Dict{T=>[(μ,ρ),...]}` 转换为 `(μ,ρ)` 数组形式 |
| `find_cep` | 函数 | 给定温度→曲线的字典，返回 CEP 近似值 |

## `CEPResult`

字段：
- `has_cep::Bool`
- `T_cep_MeV`、`mu_cep_MeV`、`rho_cep::Union{Nothing, Float64}`：若无法估计则为 `nothing`
- `bracket::Union{Nothing, Tuple{Float64, Float64}}`：S 形存在/消失的温度区间 `(T_low, T_high)`
- `details::Dict{Symbol, Any}`：附带调试信息（如自旋odal 值）

`CEPResult()` 返回默认的“未找到 CEP”状态。

## `build_curves(grouped)`

- 输入 `Dict{Float64, Vector{Tuple{Float64, Float64}}}`（多数情况下来自 `PhaseTransition.group_curves_by_temperature`）；
- 过滤掉样本数不足的温度；
- 输出 `Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}`，便于后续计算。

## `find_cep(curves; min_points=6, eps=1e-6, temp_tol=0.01, curve_fetcher=nothing, auto_fetch=:default)`

1. 按温度排序，逐一调用 `detect_s_shape`。
2. 记录最后一个呈现 S 形的温度 `T_low` 及其自旋odal。
3. 寻找之后第一个不再呈现 S 形的温度 `T_high`。
4. 若成功获得 `(T_low, T_high)` 区间，则输出：
   - `T_cep ≈ (T_low + T_high)/2`
   - `μ_cep`、`ρ_cep` 取自 `T_low` 处两自旋odal的平均值

示例：
```julia
grouped = PhaseTransition.group_curves_by_temperature(rows)
curves = build_curves(grouped)
result = find_cep(curves)  # auto_fetch=:default → 自动补样并写回 CSV
if result.has_cep
    @info "CEP" T=result.T_cep_MeV mu=result.mu_cep_MeV
end
```

若需要纯粹的离散温度判定，可显式关闭自动补样：`find_cep(curves; auto_fetch=nothing)`。

### 自动补充缺失温度

- `curve_fetcher`：显式传入 `(T::Float64) -> (mu_vals, rho_vals)` 回调，可与脚本层工具（如 `plot_phase_diagram` 中的 `_fetch_curve_via_trho!`）对接。
- `auto_fetch`：若 `curve_fetcher` 为空，默认值 `:default` 会触发自动补样：
    1. 根据现有 `curves` 在低斜率区（使用 `AdaptiveRhoRefinement.suggest_refinement_points`）生成额外 ρ 网格；
    2. 调用 `TrhoScan.run_trho_scan` 计算缺失温度；
    3. 把新的原始行追加到 `data/outputs/results/pnjl/trho_scan.csv`，并把精简 `(T, μ, ρ, ξ)` 写入 `data/processed/results/pnjl/curves.csv`；
    4. 将该温度的曲线注入 `curves`，供后续二分继续使用。

可通过传入 `NamedTuple`/`Dict` 来覆盖默认参数；若希望彻底关闭自动补样，则设置 `auto_fetch = nothing`。

```julia
result = find_cep(curves; temp_tol=0.01, auto_fetch=(
        xi = 0.0,
        rho_values = 0.0:0.02:3.0,
        seed_path = PNJL.SeedCache.DEFAULT_SEED_PATH,
        processed_path = "data/processed/results/pnjl/curves.csv",
        adaptive_config = (slope_tol = 3.0, min_gap = 0.003),
))
```

| 键 | 默认值 | 说明 |
| --- | --- | --- |
| `xi` | `0.0` | 运行 TrhoScan 时使用的 ξ；应与 `curves` 的来源一致 |
| `rho_values` | `TrhoScan.DEFAULT_RHO_VALUES` | 基础 ρ 网格，会与自适应候选合并 |
| `seed_path` | `SeedCache.DEFAULT_SEED_PATH` | 连续种子文件 |
| `xi_tol` | `1e-6` | Loader 过滤 CSV 时允许的 ξ 误差 |
| `p_num` / `t_num` | `24` / `8` | 透传给 `run_trho_scan` 的积分节点 |
| `solver_kwargs` | 空 NamedTuple | 额外 NLsolve 选项，例如 `(; method=:trust_region)` |
| `output_path` | `data/outputs/results/pnjl/trho_scan.csv` | 若不为 `nothing`，会把新行追加到该原始 CSV |
| `processed_path` | `data/processed/results/pnjl/curves.csv` | 若不为 `nothing`，会把 `(T,μ,ρ,ξ)` 追加到该文件 |
| `adaptive` | `true` | 是否启用自适应 ρ 补点 |
| `adaptive_config` | `AdaptiveRhoConfig()` | 传入 `PNJL.AdaptiveRhoRefinement` 的参数，或其 `NamedTuple`/`Dict` 形式 |
| `adaptive_window` | `5.0` | 仅使用距离目标温度不超过该值 (MeV) 的旧曲线来生成补点；`nothing` 表示全部 |
| `rho_digits` | `6` | 合并 ρ 网格时的舍入精度 |
| `runner` / `loader` / `generator` | `run_trho_scan` / 内置 CSV 解析 / `nothing` | 高级钩子：可替换执行/解析逻辑或直接返回曲线（测试用） |

> **提示**：默认配置已把自动补样写入主要 CSV，可直接复用；若不想持久化，请把 `output_path` 或 `processed_path` 设为 `nothing`。

> **注意**：自动补样同样只给出 CEP 的近似范围；建议在得到最终 bracket 后，再面向该温区运行高分辨率的 TrhoScan/Maxwell 流程以确认结果。