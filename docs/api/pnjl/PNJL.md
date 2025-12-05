# PNJL 模型 API 参考

本文件概述仓库中已实现的 PNJL（多味道 Nambu–Jona-Lasinio + Polyakov）模型核心 API、模块与使用示例，方便在代码与脚本中调用与集成。

**注意**：文档基于当前代码实现（位于 `src/pnjl`），单位说明如下：温度/化学势的外部输入通常以 MeV 为单位（参数名含 `_mev`），内部换算使用 `ħc_MeV_fm` 转换为 `fm⁻¹`。

**概览**
- **模块**：`PNJL`（顶级），包含子模块 `PNJL.AnisoGapSolver`、`PNJL.SeedCache`、`PNJL.SinglePointSolver`。
- **用途**：求解各向异性 PNJL 能隙方程（T-μ 或 T-ρ 定边界条件），计算热力学量（压力、归一化密度、熵、能量），并提供种子（seed）查找与生成工具以暖启动求解器。

**目录**
- 模块概览
- API 参考（逐模块）
- 种子表（CSV）格式与生成脚本
- 使用示例
- 注意事项与实现细节

---

**模块：PNJL**
- 位置：`src/pnjl/PNJL.jl`。
- 作用：汇总并包含子模块：`SeedCache.jl`、`solvers/AnisoGapSolver.jl`、`SinglePointSolver.jl`。

**模块：PNJL.SeedCache**
- 文件：`src/pnjl/SeedCache.jl`
- 导出符号：`SeedPoint`, `load_seed_table`, `find_initial_seed`, `seed_state`, `DEFAULT_SEED_PATH`

- 关键类型与函数：
  - `struct SeedPoint`
    - 字段：`T, mu, rho, xi, phi_u, phi_d, phi_s, Phi1, Phi2, mu_u, mu_d, mu_s, source`
    - 说明：用来表示 CSV 中的一行种子数据（温度以 `fm⁻¹` 存储，CSV 中为 `MeV`，加载时会转换）。

  - `load_seed_table(; path::AbstractString = DEFAULT_SEED_PATH, force::Bool = false)`
    - 作用：加载 CSV 文件并缓存为 `Vector{SeedPoint}`。
    - 返回：`Vector{SeedPoint}`。

  - `seed_state(seed::SeedPoint)`
    - 返回：`(x, mu, phi, Polyakov)`，其中 `x` 是解向量（长度 8），`mu`/`phi`/`Polyakov` 为解的分量。

  - `find_initial_seed(params; weights = DEFAULT_WEIGHTS, path::AbstractString = DEFAULT_SEED_PATH, k_neighbors::Int = 3)`
    - 参数 `params`：可包含 `:T`, `:T_fm`, `:T_mev`（必需），以及 `:mu` / `:mu_fm` / `:mu_mev` 或 `:rho`，和可选的 `:xi`。
    - 返回一个包含字段的 NamedTuple：
      - `seed`：最近的 `SeedPoint`，
      - `distance`：距离值，
      - `used_axes`：匹配的坐标轴数，
      - `state`：`seed_state(seed)`，
      - `neighbors`：邻居列表，
      - `blended_state`：基于 k 个邻居的逆距离加权融合结果（若存在）。
    - 行为：优先使用 KD-tree 空间查找（mu 模式或 rho 模式），若无则回退为线性扫描。

  - `DEFAULT_WEIGHTS`：用于 KD-tree 坐标尺度化的权重元组 `(T, mu, rho, xi)`。

**模块：PNJL.AnisoGapSolver**
- 文件：`src/pnjl/solvers/AnisoGapSolver.jl`
- 导出符号：`SolverResult`, `solve_fixed_rho`, `solve_fixed_mu`, `DEFAULT_RHO_GUESS`, `DEFAULT_MU_GUESS`

- 关键类型：
  - `struct SolverResult`
    - 字段：`mode::Symbol` (`:rho` 或 `:mu`), `converged::Bool`, `solution::Vector{Float64}`, `mu_vec::SVector{3, Float64}`, `pressure::Float64`, `rho::Float64`, `entropy::Float64`, `energy::Float64`, `iterations::Int`, `residual_norm::Float64`, `xi::Float64`。

- 主要函数：
  - `solve_fixed_mu(T_mev::Float64, mu_fm::Float64; xi::Float64=0.0, seed_state::AbstractVector=DEFAULT_MU_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, nlsolve_kwargs...)`
    - 描述：在给定温度 `T_mev`（MeV）与平均化学势 `mu_fm`（`fm⁻¹`）下，求解 5 维场变量（φ_u, φ_d, φ_s, Φ1, Φ2）。
    - 可选参数 `xi`：各向异性参数；`seed_state`：5 元初始猜测；`p_num/t_num`：积分节点数；`nlsolve_kwargs...`：传递到 `NLsolve.nlsolve` 的额外参数（例如 `method`, `xtol` 或 `show_trace` 等）。
    - 返回：`SolverResult`（包含热力学量与收敛信息）。

  - `solve_fixed_rho(T_mev::Float64, rho_target::Float64; xi::Float64=0.0, seed_state::AbstractVector=DEFAULT_RHO_GUESS, p_num::Int=..., t_num::Int=..., nlsolve_kwargs...)`
    - 描述：在给定温度与目标归一化密度 `rho_target` 下求解，解向量包含 8 个元素（5 个场项 + 3 个化学势分量）。
    - 返回：`SolverResult`。

  - `calculate_pressure(...)`, `calculate_core(...)`, `calculate_rho(...)`, `calculate_thermo(...)` 等内部函数用于构造残差、计算梯度与热力学量。`calculate_core` 与 `calculate_rho` 已设计为 ForwardDiff 友好（保留 Dual 元素类型）。

  - `vacuum_integral(mass)`：解析形式的真空项积分实现，用以替代数值截断积分，避免对 Λ 节点的依赖。

  - `cached_nodes(p_num, t_num)`：热项（高斯–勒让德）节点的缓存机制，避免频繁重建积分网格。

**模块：PNJL.SinglePointSolver**
- 文件：`src/pnjl/SinglePointSolver.jl`
- 导出符号：`describe_solver`, `run_single_point`

- 主要函数与行为：
  - `describe_solver()`
    - 返回一个描述字典（ID、说明、参数 schema），可用于构建 API/服务描述。

  - `run_single_point(params::Union{NamedTuple, AbstractDict, Dict}; seed_path::AbstractString = DEFAULT_SEED_PATH, solver_kwargs...)`
    - 参数 `params`：必须包含 `T_mev`（浮点），并包含 `mu_mev` 或 `rho` 中的一个；可选 `xi`。
    - 行为：
      1. 规范化请求参数（单位转换），
      2. 尝试通过 `SeedCache.find_initial_seed` 获取暖启动种子，
      3. 根据模式调用 `solve_fixed_mu` 或 `solve_fixed_rho`，
      4. 返回一个 `Dict`，包含 `status`, `mode`, `converged`, `iterations`, `residual_norm`, `thermo`（pressure, rho, entropy, energy），`solution_vector`，`mu_components`，以及使用的 `seed` 与 `seed_distance`（如果有）。

---

**种子表（CSV）格式与生成**
- 默认路径：`data/raw/pnjl/seeds/sobol_seed_table.csv`（常量 `PNJL.SeedCache.DEFAULT_SEED_PATH` 指向该位置）。
- 列（CSV header）：
  - `T_MeV, mu_MeV, rho, xi, phi_u, phi_d, phi_s, Phi1, Phi2, mu_u_MeV, mu_d_MeV, mu_s_MeV`
- 生成脚本：`scripts/pnjl/generate_seed_table.jl`
  - 使用 LHS/Sobol（实现为 Latin Hypercube）在给定参数区间采样，调用 `AnisoGapSolver.solve_fixed_mu` 逐点求解并输出成功的收敛解到 CSV。
  - 可配置参数（样本数、T/mu/xi 区间、积分节点数等），示例运行：

```pwsh
julia --project scripts/pnjl/generate_seed_table.jl --samples=200 --output=data/raw/pnjl/seeds/sobol_seed_table.csv
```

**种子查找微基准**
- 脚本 `scripts/pnjl/benchmark_seed_lookup.jl` 对 `find_initial_seed` 的热路径（KD-tree 查询、邻居融合）做了 `BenchmarkTools.@btime` 微基准。典型热路径延迟约为几十至百微秒（取决于 k）。

---

**使用示例**

1) 直接调用单点求解（通过 `SinglePointSolver.run_single_point`）：

```julia
using PNJL
using PNJL.SinglePointSolver: run_single_point

params = Dict(:T_mev => 100.0, :mu_mev => 50.0, :xi => 0.0)
out = run_single_point(params)
println(out["status"], ", converged = ", out["converged"]) 
```

2) 直接使用求解器（更细控制）：

```julia
using PNJL.AnisoGapSolver: solve_fixed_mu

T_mev = 100.0
mu_fm = 0.05  # fm^-1
res = solve_fixed_mu(T_mev, mu_fm; xi=0.0, p_num=24, t_num=8)
if res.converged
    println("pressure = ", res.pressure, ", rho = ", res.rho)
end
```

3) 查询种子：

```julia
using PNJL.SeedCache: find_initial_seed
seed_info = find_initial_seed((T=100.0, mu_mev=50.0));
println(seed_info.seed, seed_info.distance)
```

---

**实现细节与注意事项**
- 单位：外部 API 偏好 `MeV`（参数名带 `_mev`），内部使用 `fm⁻¹`（通过 `ħc_MeV_fm` 换算）。
- 自动微分：`calculate_core` / `calculate_rho` 保留输入元素类型并用 ForwardDiff 构造 SVector，因此与 `ForwardDiff.Dual` 兼容，用于 NLsolve 的前向差分。请勿在外部强行将这些输入转换为只含 `Float64` 的结构，否则会失去 AD 兼容性。
- 真空项：已用解析表达式 `vacuum_integral(mass)` 代替数值截断积分，从而移除了对 Λ 节点缓存的必要性。
- 积分节点：热项的 Gauss–Legendre 节点由 `cached_nodes(p_num, t_num)` 缓存，可通过参数 `p_num`/`t_num` 控制精度与速度。默认 `t_num` 已在仓库中被调至较小值（例如 8），以在性能与精度间取得平衡；如需严格精度可增大 `t_num`。
- 求解器：`solve_fixed_mu` 与 `solve_fixed_rho` 调用 `NLsolve.nlsolve(...; autodiff=:forward)`，可通过 `nlsolve_kwargs...` 传递额外参数（例如 `method`, `xtol`, `ftol` 或 `linesearch`），特别是在 fixed-ρ 模式下可传入 `linesearch = LineSearches.BackTracking()` 以提升稳健性。
- 种子融合：`find_initial_seed` 支持 k-NN 逆距离加权融合（若 k>1）。融合开销极小（微基准显示融合成本在 μs 级），但可显著改善暖启动平滑性与求解收敛。

---

若希望把此文档自动生成成网站页面或把函数签名与注释同步为更详细的 API 文档，我可以：
- 将此文件转为 docs 网站页（例如 MkDocs/Documenter.jl 格式）。
- 提取并补充每个导出函数的参数类型注释与更丰富的示例（包含 `nlsolve_kwargs` 常用配置示例）。

请告诉我下一步想要的格式或更详细内容（例如增加 `Trho` 种子生成说明、示例输入/输出 JSON 模板或把文档集成到 README）。
