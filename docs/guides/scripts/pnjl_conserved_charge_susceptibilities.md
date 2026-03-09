# PNJL 守恒荷广义磁化率脚本入口

本指南对应统一入口脚本：

- [scripts/pnjl/run_conserved_charge_susceptibilities.jl](scripts/pnjl/run_conserved_charge_susceptibilities.jl)

目标是给用户一个稳定入口，统一输出：

- `chi_BQS`
- `cumulant_BQS`
- `baryon_Ssigma`
- `baryon_kappa_sigma2`

---

## 1. 推荐用法

单点计算：

```powershell
julia --project=. scripts/pnjl/run_conserved_charge_susceptibilities.jl --T=150 --muB=0 --muQ=0 --muS=0 --V=64
```

带 CSV 输出：

```powershell
julia --project=. scripts/pnjl/run_conserved_charge_susceptibilities.jl --T=150 --muB=300 --muQ=40 --muS=20 --V=64 --output=data/outputs/results/pnjl/susceptibilities/example_point.csv
```

默认会在终端打印：

- `chi1_B`, `chi2_B`, `chi3_B`, `chi4_B`
- `chi2_Q`, `chi2_S`
- `chi11_BQ`, `chi11_BS`, `chi11_QS`
- `C200`, `C020`, `C002`, `C110`, `C101`, `C011`
- `Ssigma`, `kappa_sigma2`

---

## 2. API 对应关系

脚本内部统一调用：

- `chi_BQS(T, muB, muQ, muS; orders=(i,j,k))`
- `cumulant_BQS(T, muB, muQ, muS, V; orders=(i,j,k))`
- `baryon_Ssigma(T, muB)`
- `baryon_kappa_sigma2(T, muB)`

因此如果后续你需要从脚本迁移到 Julia 代码里，调用对象是一致的。

---

## 3. 关于脚本入口过多

目前仓库里的脚本确实偏多，但不是所有脚本都应该被视为“用户入口”。

更合理的规范是：

- `scripts/pnjl/run_*.jl`
  - 稳定用户入口
  - 应该尽量少，命名明确，帮助文本完整
- `scripts/analysis/`
  - 后处理分析脚本
  - 不承诺稳定 CLI
- `scripts/dev/`
  - 开发期导出、一次性迁移、回归辅助脚本
  - 不应被文档当作正式入口推荐
- `scripts/debug/`
  - 诊断/排障脚本
  - 默认不面向普通用户

结论：

- 是，需要规范化管理。
- 但不建议“删成很少几个巨型脚本”。
- 更好的做法是把稳定入口集中到 `run_*.jl`，其余目录明确降级为内部工具或分析工具。

---

## 4. 建议的后续治理动作

建议后续按以下顺序收敛：

1. 新增功能优先补到已有稳定入口，避免继续散落到新脚本。
2. 对 `scripts/pnjl/` 做一次入口盘点，只保留少量稳定 `run_*.jl` 作为用户文档推荐对象。
3. 给 `scripts/dev/` 和 `scripts/analysis/` 明确加上“非稳定入口”口径。
4. 在 `docs/guides/QUICKSTART.md` 或总脚本指南里维护一个“推荐入口清单”。

当前这条守恒荷广义磁化率线已经按这个方向落地：

- 正式推荐入口是 [scripts/pnjl/run_conserved_charge_susceptibilities.jl](scripts/pnjl/run_conserved_charge_susceptibilities.jl)
- 底层 API 是 `chi_BQS / cumulant_BQS / baryon_Ssigma / baryon_kappa_sigma2`