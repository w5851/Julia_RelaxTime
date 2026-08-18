# T200 Tau-u Spike Analysis

这是 T200 `tau_u` 双窗口异常的诊断 evidence package。它沿着

```text
denominator -> propagator -> sigma/rate -> tauinv -> tau -> transport
```

组织读图和机制证据，不是 production result、regression truth 或 reference promotion。

## Contents

- `t200-dual-window-tauu-spikes-analysis.md`：主分析、通道归因和分母机制证据链。
- `t200-denominator-chain-for-readers.md`：面向读者的最小读图和复现入口。
- `tauu_pos_*.png`：分母复轨迹、传播子强度、`s/t` 敏感性、项切换和 `Re(M08)` 零交叉图。
- `manifest.json`：迁移后的文件 hash、脚本入口和 provenance 边界。

## Scope

- case：`T=200`、`muB=0` 的 `tau_u` 双窗口及相关通道。
- status：`diagnostic_only`。
- 主要机制口径：不同窗口可能由 simple `1-4KPi` 过零或 mixed `detM` 近零/项间重分配主导，不能将两个窗口压缩为单一机制。
- 图像只呈现已有诊断结果；本次迁移没有重跑数值、重绘图片或修改正式 result/figure root。

## Reproduction Entry Points

数据和机制诊断脚本位于 `scripts/analysis/relaxtime/t200_*.jl`；绘图脚本的默认输出目录已统一为本目录。代表性入口包括：

```powershell
julia --project=. scripts/analysis/relaxtime/t200_dual_window_full_chain_decomposition.jl
julia --project=. scripts/analysis/relaxtime/t200_imag_path_evidence.jl
python scripts/analysis/relaxtime/plot_t200_denominator_chain.py
```

完整参数、外部临时输入和逐层证据仍以主分析文档及 `manifest.json` 为准。
