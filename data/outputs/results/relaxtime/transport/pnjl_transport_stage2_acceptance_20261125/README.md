# Stage 2 作者集中签收（2026-11-25）

作者一次性审阅并签收 8 个图件、16 张 PNG（彩色 + 灰度）。本目录存放逐图哈希绑定的签收文件，
Stage 3 已用它们零改码导出全部 8 个矢量 PDF。

## 作者原话

> 我已审阅八张新图（C2 相图、两张 T–μB 带、π/K 谱、弛豫时间、输运系数、固定温度扫描、Mott 机制图），可以正式接受并完成矢量交付。

## 签收范围与矢量产物

| # | 资产 | PNG 复核 case | 矢量交付 case |
|---|---|---|---|
| 1 | `pnjl.transport.phase.native_issue130_c2` | `native_issue130_c2_20261124_stage1_png_review` | `native_issue130_c2_20261125_vector` |
| 2 | `meson_freezeout_andronic_Tmu_envelope_20261018` | `andronic_envelope_20261124_Tmu_stage1_png_review` | `andronic_envelope_20261125_Tmu_vector` |
| 3 | `meson_freezeout_Tmu_envelope_20261018` | `freezeout_envelope_20261124_Tmu_stage1_png_review` | `freezeout_envelope_20261125_Tmu_vector` |
| 4–8 | `ward_moller_*`（5 图） | `ward_moller_manuscript_20261125_stage1r3_png_review` | `ward_moller_manuscript_20261125_vector` |

签收以 manifest 与 PNG 的 SHA-256 逐字节绑定；`pdf_export_authorized = true`，
`manuscript_eligible` 保持 `false`（留待 Stage 4 由作者单独授权）。

## 关于 r2 → r3

ward_møller 的 PNG 复核 case 初版（`..._stage1r2_png_review`）渲染时漏传 `--review-results-dir`，
缺少交付合约要求的 `production figure SOP snapshot` 输入，被 `add_vector_pdf()` 在交付时拦下。
重出 `..._stage1r3_png_review` 补齐该输入后，**10 张 PNG 的 SHA-256 与已签收版本逐字节一致**
（已程序核对：不一致条目 = 0），故仅把 manifest 绑定从 r2 改挂到 r3，作者原话与 PNG 哈希均未改动。
r2 目录随后删除，避免出现两个等价复核 case。

## 复核依据

- 视觉终审：8 个图件分 3 批送审（C2 + 两张 T–μB、ward_møller 彩版、ward_møller 灰度版），
  机制图因灰度堆叠柱不可辨返工一版，复审全通过。
- 采样密度质询：作者质疑机制图点位偏少，复核确认 4 点即原生 Δξ=0.01 网格、2 根柱即穿越前后断点，
  u+d 扇区携带 97–98% 电流故其余物种贴零，**不需要补点重算**。
- 数字口径：正文 `+1.40% / −1.35% / +0.05%` 已与冻结 `mott_window_summary.csv` 对齐。
- 矢量校验：4 个 manifest 全部 `validate_plot_artifact.py` OK；8 个 PDF 均单页、0 个 Type 3 字体、0 个未嵌入字体。
