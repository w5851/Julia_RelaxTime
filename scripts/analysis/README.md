# scripts/analysis/

探索性分析与诊断脚本目录（非 CI 路径，仅手动运行）。

## 当前内容

- `convergence/`：数值收敛性分析相关脚本与说明
- `scan_mott_meson_vs_xi_mu0.jl`：在 `mu=0` 条件下扫描 Mott 点与介子质量随 `xi` 的变化
- `mott_reference_mapping.jl`：生成/检查 Mott 参考映射数据

## 使用说明

- 这些脚本用于研究与诊断，不作为测试准入门禁。
- 建议从仓库根目录执行并显式指定项目环境。

```bash
julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl
julia --project=. scripts/analysis/mott_reference_mapping.jl
```
