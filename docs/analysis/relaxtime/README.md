# Relaxation-Time Analysis Index

本目录保存弛豫时间、散射率、传播子机制和历史比较相关的诊断分析。正式 production CSV/figure root 仍在 `data/outputs`；这里的分析包默认是 diagnostic 或 historical evidence。

## Transport Mechanism Line

### Phase-guided transport

| 阶段 | 当前路径 | 作用 |
| --- | --- | --- |
| v1 tau-first analysis | `phase_guided_transport_p128_xi001_analysis/` | 识别 tau/channel-rate 突变、下游输运响应和 denominator-chain 候选 |
| v2 pole-sensitive rendering | `phase_guided_transport_v2_pole_sensitive_rendering/` | 在 v2 on-shell-kernel production 上迁移机制证据，加入 v2 定点诊断、pole mask 和一阶分支保护 |

v2 的 `tables/window_classification.csv` 明确引用 v1 的 `mechanism_window_summary.csv`，所以两者属于同一逻辑线的连续阶段；v2 不应覆盖或重写 v1。

### T200 tau-u spike package

以下根目录文件共同属于一个 T200 双窗口机制包：

- `t200-dual-window-tauu-spikes-analysis.md`：主分析和证据链。
- `t200-denominator-chain-for-readers.md`：读者导览和最小复现入口。
- `tauu_pos_*.png`：通道、传播子分母、`s/t` 敏感性和零交叉图。

建议目标逻辑路径为 `transport/t200_tauu_spikes/`。当前先保留物理路径，因为多份 T200 分析脚本、文档和图 registry 仍直接指向这些文件。

## Historical Evidence

| 逻辑类别 | 当前路径 | 语义 |
| --- | --- | --- |
| literature comparison | `literature_comparison/` | 文献 target 对照和历史插值比较，不是 strict validation gate |
| meson-density plot review | `meson_density/freezeout_kminus_piminus_mu_pi_100_analysis/` | freeze-out workflow 与 digitized target 的人工 plot review |
| Julia/Fortran validation history | `validation/` | 历史 swap-validated 对照图，不是当前 regression truth |

这些目录有保存价值，但应在总索引中标为 `historical`，不应和当前 phase-guided transport 机制包混合。

## Package Contract

新的 transport 分析包优先采用：

```text
<case>_analysis/
├── README.md
├── manifest.json
├── figures/
└── tables/
```

其中 README 负责范围和结论边界，manifest 负责输入/输出 hash 和生成 provenance，tables 负责可审计数值证据。图像不能单独承载机制结论。
