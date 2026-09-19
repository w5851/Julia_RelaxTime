# Issue #130 RS `publication_clean_v3` review candidate

## 目的与边界

本包是 `publication_clean_v2_residual_smoothed` 的独立 display-only 派生，供作者审查。
它只对 `mode_b / T=200 MeV / muB=0 MeV / xi=0.36` 的四条 tau 曲线做一次显式局部
线性插值；v2、v2 residual-smoothed、raw CSV、production registry 和 solver 输出均不修改。

这不是传播子有限宽度、分母截断、重新求解或新的物理正则化，也不把该局部结构自动归因于
传播子小分母。

## 变换规则

- 左右锚点：`xi=0.35` 与 `xi=0.37`，锚点值保留父级显示值；
- 目标点：`xi=0.36`；按目标点在锚点区间中的线性位置插值；
- phase gate：区间内必须为 `crossover` 且不含 rendered first-order gap；
- 受影响 observable：`tau_u`, `tau_d`, `tau_ubar`, `tau_dbar`；
- `raw_value`、`v2_clean_value` 和父级 `clean_value` 均在独立审计表中保留。

## Provenance

- parent artifact: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2_residual_smoothed`
- parent manifest SHA256: `b3da0121e21fd94e0734b6aea7c5cd61f3d5a08bed69fac2daabd53be4e0b11e`
- parent plot manifest SHA256: `44c2731fbdb606349aaf33494f0481af8882ab2e5fe169611664c1ea0369c11e`
- adjustment recipe: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_v3_residual_smoothing/tables/tau_xi036_display_adjustments.csv`
- adjustment recipe SHA256: `426c9bd8ce73fbd69ca7cce4111e5f4efc0383583ff06e90f875a47c8808965e`
- derived solver called: `false`
- canonical/raw data modified: `false`
- production write: `false`
- manuscript eligible: `false`

## 局部变换审计

| observable | parent display value | v3 display value | relative change | log residual before → after |
| --- | ---: | ---: | ---: | ---: |
| tau_u | `2.2623236669554796` | `2.2825356240283448` | `0.00893416` | `-0.0088279136091645727` → `6.6570918986985639e-05` |
| tau_d | `2.2623236669554796` | `2.2825356240283448` | `0.00893416` | `-0.0088279136091645727` → `6.6570918986985639e-05` |
| tau_ubar | `2.2623236669566138` | `2.2825356240295021` | `0.00893416` | `-0.0088279136091695687` → `6.6570918987762795e-05` |
| tau_dbar | `2.2623236669566138` | `2.2825356240295021` | `0.00893416` | `-0.0088279136091695687` → `6.6570918987762795e-05` |

当前四条曲线在 `muB=0` 下由对称性几乎重合；本包仍逐 observable 记录变换，避免把对称性
假设隐藏在一个未审计的共享数值中。局部前后对比图见
`audit/tau_xi036_before_after.png`。

## 审查边界

1. v3 只适合图形审查和排版候选，不能替代 raw/v2 数据作精确数值、导数、峰值或机制结论。
2. 若正式论文采用 v3 图，应在内部稿件 provenance 中保留本包及 adjustment map；不得静默覆盖 v2。
3. 一阶 gap、端点和跨分支连接沿用父级语义，未填补任何 gap。
4. `manuscript_eligible=false` 保持到作者明确审查通过为止。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v3.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v3.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v3.py
```

本包生成 publication figures：72 张，完整覆盖父级的 72 张 figure 集。
