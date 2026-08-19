# Relaxation-Time Analysis Index

本目录保存弛豫时间、散射率、传播子机制和历史比较相关的诊断分析。正式 production CSV/figure root 仍在 `data/outputs`；这里的分析包默认是 diagnostic 或 historical evidence。

## Transport Mechanism Line

### Phase-guided transport

| 阶段 | 当前路径 | 作用 |
| --- | --- | --- |
| v1 tau-first analysis | `phase_guided_transport/phase_guided_transport_p128_xi001_analysis/` | 识别 tau/channel-rate 突变、下游输运响应和 denominator-chain 候选 |
| v2 pole-sensitive rendering | `phase_guided_transport/phase_guided_transport_v2_pole_sensitive_rendering/` | 在 v2 on-shell-kernel production 上迁移机制证据，加入 v2 定点诊断、pole mask 和一阶分支保护 |

v2 的 `tables/window_classification.csv` 明确引用 v1 的 `mechanism_window_summary.csv`，所以两者属于同一逻辑线的连续阶段；v2 不应覆盖或重写 v1。

该逻辑线的物理分组入口见 [`phase_guided_transport/`](phase_guided_transport/README.md)。目录迁移只改变 namespace；包内生成时 manifest、图 manifest、旧路径快照和 provenance 保持原样。

### T200 tau-u spike package

以下文件共同属于 `transport/t200_tauu_spikes/` T200 双窗口机制包：

- `transport/t200_tauu_spikes/t200-dual-window-tauu-spikes-analysis.md`：主分析和证据链。
- `transport/t200_tauu_spikes/t200-denominator-chain-for-readers.md`：读者导览和最小复现入口。
- `transport/t200_tauu_spikes/tauu_pos_*.png`：通道、传播子分母、`s/t` 敏感性和零交叉图。

该目录的物理迁移已完成；历史 figure registry 仍作为迁移前快照保留，不作为当前路径索引。

## Historical Evidence

| 逻辑类别 | 当前路径 | 语义 |
| --- | --- | --- |
| literature comparison | `historical/literature_comparison/` | 文献 target 对照和历史插值比较，不是 strict validation gate |
| meson-density plot review | `historical/meson_density/freezeout_kminus_piminus_mu_pi_100_analysis/` | freeze-out workflow 与 digitized target 的人工 plot review |
| Julia/Fortran validation history | `historical/validation/` | 历史 swap-validated 对照图，不是当前 regression truth |

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
