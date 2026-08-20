# scripts/analysis/

探索性分析与诊断脚本目录（非 CI 路径，仅手动运行）。

## 当前内容

- `convergence/`：数值收敛性分析相关脚本与说明
- `scan_mott_meson_vs_xi_mu0.jl`：在 `mu=0` 条件下扫描 Mott 点与介子质量随 `xi` 的变化
- `mott_reference_mapping.jl`：生成/检查 Mott 参考映射数据
- `relaxtime/audit_bu_meson_density_literature_alignment.jl`：面向 BU 介子数密度文献对齐审查的 fixed-point / phase normalization / charged `μ_K` 规则诊断脚本
- `relaxtime/build_phase_guided_transport_xi001_jump_analysis.py`：基于 xi=0.01 p128 phase-guided transport 正式产物生成 tau-first 突变分析包，并把 `eta_over_s` / `zeta_over_s` 作为 tau 下游响应处理
- `relaxtime/phase_guided_p128_mechanism_scan.jl`：消费 xi001 分析包生成的 mechanism window candidates，对 phase-guided transport 局部窗口执行 denominator-chain / rate-band 机制深拆并写回分析包表格
- `relaxtime/build_phase_guided_pole_sensitive_rendering.py`：验证 v1→v2 tau/rate 机制迁移门槛，生成不改写正式产物的内部极点敏感审计，以及隐藏数值修正痕迹、用星号标示一阶相变点的多曲线论文候选图

## 使用说明

- 这些脚本用于研究与诊断，不作为测试准入门禁。
- 建议从仓库根目录执行并显式指定项目环境。
- xi001 tau-first 突变分析只消费已入库正式产物，不重跑 production。
- xi001 transport 分析包的完整生成顺序是：先运行 Python 生成候选窗口，再运行 Julia 机制深拆，最后重新运行 Python 刷新主报告和 manifest。

```bash
julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl
julia --project=. scripts/analysis/mott_reference_mapping.jl
python scripts/analysis/relaxtime/build_phase_guided_transport_xi001_jump_analysis.py
julia --project=. scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl --case-name first_canonical_v1_p128_xi001_validated_anchored_prod_v1 --candidate-csv docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis/tables/mechanism_window_candidates.csv --out-dir docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis/tables --integration-mode semi_infinite
python scripts/analysis/relaxtime/build_phase_guided_transport_xi001_jump_analysis.py
python scripts/analysis/relaxtime/build_phase_guided_pole_sensitive_rendering.py
```
