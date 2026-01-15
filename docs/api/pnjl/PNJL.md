# PNJL 模型 API 参考

本文件概述仓库中已实现的 PNJL（多味道 Nambu–Jona-Lasinio + Polyakov）模型核心 API、模块与使用示例。

**单位约定**：
- 求解器 `PNJL.solve(...)` 使用自然单位 `fm⁻¹`（`T_fm`, `μ_fm`）。
- 扫描脚本/扫描模块通常使用 MeV（`T_MeV`, `μ_MeV`, `μB_MeV`），再通过 `ħc_MeV_fm` 换算到 `fm⁻¹`。

**注意**：仓库已迁移到新架构（位于 `src/pnjl`）。旧版 `SeedCache/AnisoGapSolver/SinglePointSolver` 已移除或废弃；以本文件与 `docs/api/pnjl/*` 为准。

---

## 模块概览

**模块：PNJL**

- 位置：`src/pnjl/PNJL.jl`
- 作用：聚合并导出以下子模块能力：
  - `PNJL.Integrals` / `PNJL.Thermodynamics`：积分与热力学量
  - `PNJL.ConstraintModes` / `PNJL.SeedStrategies` / `PNJL.ImplicitSolver`：求解器主链路
  - `PNJL.ThermoDerivatives`：热力学导数/体粘滞相关导数
  - `PNJL.TmuScan` / `PNJL.TrhoScan` / `PNJL.DualBranchScan`：扫描与相变分析辅助
  - `PNJL.PhaseTransition`：S 形检测、Maxwell 构造、crossover 扫描等

## 核心入口（新架构）

### 平衡求解

```julia
using PNJL

res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=0.0)
res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=0.0, seed_strategy=PNJL.MultiSeed())
```

更多求解器细节见：
- `docs/api/pnjl/ImplicitSolver.md`
- `docs/api/pnjl/SeedStrategies.md`

### 扫描

- `PNJL.run_tmu_scan`：T-μ 网格扫描（见 `docs/api/pnjl/TmuScan.md`）
- `PNJL.run_trho_scan`：T-ρ 网格扫描（见 `docs/api/pnjl/TrhoScan.md`）
- `PNJL.run_dual_branch_scan`：一阶相变区域双分支扫描（见 `docs/api/pnjl/DualBranchScan.md`）

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

## 使用示例

```julia
using PNJL

T_fm = 150.0 / PNJL.ħc_MeV_fm
μ_fm = 0.0
xi = 0.0

res = PNJL.solve(PNJL.FixedMu(), T_fm, μ_fm; xi=xi, seed_strategy=PNJL.MultiSeed())
@show res.converged res.omega res.masses
```

---

## 注意事项

- 若你发现某些点“收敛但落在高 Ω 的非物理/亚稳态分支”，优先用 `MultiSeed()` 或在扫描中使用 `PhaseAwareContinuitySeed(...; bootstrap_multiseed=true)`。
- 一阶相变区域更适合使用 `DualBranchScan` 做物理分析（两分支显式比较 Ω），而不是只靠单分支连续性。

---

若希望把此文档自动生成成网站页面或把函数签名与注释同步为更详细的 API 文档，我可以：
- 将此文件转为 docs 网站页（例如 MkDocs/Documenter.jl 格式）。
- 提取并补充每个导出函数的参数类型注释与更丰富的示例（包含 `nlsolve_kwargs` 常用配置示例）。

请告诉我下一步想要的格式或更详细内容（例如增加 `Trho` 种子生成说明、示例输入/输出 JSON 模板或把文档集成到 README）。
