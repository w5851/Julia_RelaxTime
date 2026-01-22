# MesonMassWorkflow

串联：PNJL 平衡求解 → 介子质量/宽度（MesonMass）→ Mott 阈值与 gap（MottTransition）。

实现位于 [src/pnjl/workflows/MesonMassWorkflow.jl](src/pnjl/workflows/MesonMassWorkflow.jl)。

## 入口

### `solve_gap_and_meson_point`

```julia
solve_gap_and_meson_point(
    T_fm,
    mu_fm;
    xi=0.0,
    mesons=DEFAULT_MESONS,
    k_norm=0.0,
    p_num=..., t_num=...,
    seed_state=..., 
    solver_kwargs=(;),
    mass_kwargs=(;),
)
```

- **输入单位**
  - `T_fm`, `mu_fm`：fm⁻¹
  - `xi`：无量纲各向异性参数
  - `k_norm`：介子外部三动量模（fm⁻¹）；默认 0 表示静止极点
- **流程**
  1. 调用 `PNJL.solve(FixedMu(), T_fm, mu_fm; xi=xi, ...)` 得到平衡解与三味有效质量。
  2. 组装 `quark_params=(m=(u,d,s), μ=(u,d,s))` 与 `thermo_params=(T, Φ, Φbar, ξ)`。
  3. 对 `mesons` 中每个通道调用 `MesonMass.solve_meson_mass` 求解质量与宽度。
  4. 计算 Mott 阈值与 gap：
     - 非混合通道：`threshold = mott_threshold_mass(meson, quark_params)`，`gap = mott_gap(meson, mass, quark_params)`
     - 混合通道（η/η′/σ/σ′）：
       `threshold = mott_threshold_masses(meson, quark_params)` 返回 `(uu, ss, min)`，gap 同理 `mott_gaps`。

- **返回**（NamedTuple）
  - `equilibrium`：PNJL 平衡求解输出（含 `x_state`, `masses` 等）
  - `quark_params`, `thermo_params`
  - `meson_results::Dict{Symbol,NamedTuple}`：每个介子的结果：
    - `mass`, `gamma`, `converged`, `residual`
    - 非混合：`threshold`, `gap`
    - 混合：`threshold=(uu,ss,min)`, `gaps=(uu,ss,min)`

### `build_equilibrium_params`

```julia
build_equilibrium_params(base, T_fm, mu_fm; xi=0.0)
```

将 `PNJL.solve` 的结果 `base` 转成 `(quark_params, thermo_params)`，便于在外层脚本做“先求平衡、后多次复用参数”的组织。

## 默认通道

`DEFAULT_MESONS`：

- `:pi`, `:K`, `:eta`, `:eta_prime`
- `:sigma_pi`, `:sigma_K`, `:sigma`, `:sigma_prime`

## 示例

```julia
include("src/pnjl/workflows/MesonMassWorkflow.jl")
using .MesonMassWorkflow

T = 0.15  # fm^-1
mu = 0.0  # fm^-1

res = solve_gap_and_meson_point(
    T,
    mu;
    xi=0.0,
    mesons=(:pi, :K, :eta, :eta_prime),
    p_num=16,
    t_num=8,
    mass_kwargs=(;),
)

@show res.quark_params
@show res.thermo_params
@show res.meson_results[:pi]
@show res.meson_results[:eta]
```

## 扫描脚本

仓库提供了一个基于该工作流的扫描脚本：

- [scripts/relaxtime/run_gap_meson_mass_scan.jl](scripts/relaxtime/run_gap_meson_mass_scan.jl)

用于扫描 `(T, μ_B, ξ)` 网格并写出 `scan_csv_v1` 格式 CSV（支持续跑与跳过已算点）。

## 数值健壮性

- 单个介子通道的 `MesonMass.solve_meson_mass` 若因数值原因失败（例如求解过程中出现非有限数），工作流会将该通道结果标记为 `converged=false`，并返回 `NaN/Inf` 占位，但不会中断整点/整网格扫描。
- 内部包含“多初值重试”的策略（对质量/宽度的初值做几组尝试），以减少偶发的 `NaN` 与求解失败。

## 与 Fortran 结果对照

该对照依赖外部 Fortran 工作区（非仓库默认依赖）。相关对照脚本与结论已归档：

- [docs/dev/archived/2026-01-19_MesonMass_MottTransition_Fortran_Validation.md](docs/dev/archived/2026-01-19_MesonMass_MottTransition_Fortran_Validation.md)

如需复现，对照脚本代码已归档在：

- [docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md](docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md)

示例（PowerShell）：

```powershell
# 1) 从归档文档中复制脚本 compare_meson_masses_with_fortran.jl
#    docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md
# 2) 保存为本地临时文件，例如：D:/tmp/compare_meson_masses_with_fortran.jl
julia --project=. D:/tmp/compare_meson_masses_with_fortran.jl `
  --fortran-meson-file "D:/Desktop/fortran代码/输运系数/RelaxTime/PNJL-mott-mu_T/PNJL-mu-T/quark_phase/meson.dat" `
  --t-list 150,200,210
```

注意：

- Fortran 的 `meson.dat` 列顺序在不同代码版本中可能存在“标号/通道命名差异”；脚本内置了一套与本仓库参考数据匹配的列映射。
- 混合通道（η/η′、σ/σ′）在不同实现中可能出现“本征态排序互换”的现象。脚本对 η/η′、σ/σ′ 会自动尝试交换匹配并输出提示，但若某一支仍出现较大偏差，通常意味着求解选根/混合角或积分细节差异，需要进一步诊断。
