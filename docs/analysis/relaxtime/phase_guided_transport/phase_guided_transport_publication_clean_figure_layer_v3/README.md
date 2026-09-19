# RS `publication_clean_v3` figure-layer review candidate

本目录记录 `publication_clean_v3` 的 byte-preserving publication figure mirror。
实际 PNG 位于 `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v3/`。

v3 只在 `mode_b / T=200 MeV / muB=0 MeV / xi=0.36` 的四个 tau observable 上应用
显式局部线性插值；raw、production registry、solver 输出以及 v2/v2 residual-smoothed
快照均保持不变。audit 图留在 `docs/analysis/...publication_clean_v3/audit/`，不进入公共
publication figure 目录。

`manuscript_eligible=false`。本层仅供作者审查，不能替代 raw/v2 数据作定量或机制结论。

复现：

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v3.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v3.py
```
