# P-mu 诊断入口

本页说明相图主题中的 compare-only `P-mu / Omega-mu` 双分支竞争诊断入口，面向想在 CEP 邻域对照 Maxwell 判据的使用者。

## 适用场景

当你满足以下条件时，应优先使用该入口：

- 只想做 fixed-`T` 的 hadron/quark branch competition 诊断
- 需要并排比较 `P-mu` 与 Maxwell 的转变点口径
- 不希望替换现有 production CEP 主结果

当前 v1 主要面向：

- `xi=0.0`
- `solver_backend=:models`
- `p_num=24`
- `t_num=8`
- CEP 邻域温区

## 面向用户入口

首选公开入口：

- `Models.analyze_pm_branch_competition`
- CLI: `scripts/pnjl/diagnose_pm_phase.jl`

最小 Julia 示例：

```julia
using .Models

outdir = mktempdir()

result = Models.analyze_pm_branch_competition(
    T_values=[130.9],
    mu_grid=collect(290.9:0.1:291.1),
    xi=0.0,
    solver_backend=:models,
    p_num=24,
    t_num=8,
    output_dir=outdir,
)
```

最小 CLI 示例：

```bash
julia --project=. scripts/pnjl/diagnose_pm_phase.jl --T_values=130.9 --mu_start=290.9 --mu_stop=291.1 --mu_step=0.1 --output_dir=tmp/pm_diag
```

## 输入与输出约定

- `T_values` 与 `mu_grid` 均采用 `MeV`
- `output_dir` 下当前会写出三类工件：
  - `pm_branch_scan.csv`
  - `pm_phase_summary.json`
  - `pm_vs_maxwell.csv`
- 返回值是一个命名元组，至少包含：
  - `branch_rows`
  - `temperature_summaries`
  - `comparison_rows`
  - `artifacts`

## 职责核心

该诊断入口的职责不是替换 Maxwell，而是补充一个“同配置、同温度、同积分口径”的对照视图：

1. 为每个温度构造 hadron-like / quark-like 初始 seed
2. 在固定 `mu_grid` 上分别追踪两条逻辑分支
3. 对每个 branch point 应用 acceptance 与 continuity 判据
4. 从双分支 overlap 区域提取 `Delta Omega` 零点与 metastable endpoint
5. 用同一组 branch rows 调用 `detect_s_shape` 与 `maxwell_construction` 构造 same-config Maxwell reference
6. 输出 compare-only artifact，而不改变主 pipeline 的 CEP 选择逻辑

这意味着它更接近“诊断台账”，而不是当前生产口径的替代者。

## 稳定性说明

- `analyze_pm_branch_competition` 是新加的 compare-only 入口，当前仍以 CEP 邻域小范围诊断为主
- `_pm_*` helper 虽已导出供测试和维护使用，但不建议把这些内部 helper 当作首选用户入口
- 若需要完整导出清单，请查看自动生成页

## 导出 API 全集

完整导出索引见 [generated/PMPhaseDiagnosticExports.md](docs/api/models/phase/generated/PMPhaseDiagnosticExports.md)。
