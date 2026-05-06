---
title: strict BW Stage1 扫描结果说明
archived: true
original: docs/dev/active/2026-05-02_strict BW Stage1扫描结果说明.md
archived_date: 2026-05-06
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# strict BW Stage1 扫描结果说明

更新日期：2026-05-02

目的：记录 reduced strict BW（Stage 1）正式 workflow 扫描在 `T = 208:2:220 MeV`、`muB = 0`、`xi = 0` 窗口下的最小结果，并明确它与 stable / Phase-E3 的当前关系。

相关输出：

- strict BW 扫描：
  - `data/outputs/results/relaxtime/scan/strict_bw_meson_density_scan_208_220_step2.csv`
- stable 扫描：
  - `data/outputs/results/relaxtime/scan/meson_density_scan_208_220_step2.csv`
- Phase-E3 扫描：
  - `data/outputs/results/relaxtime/scan/phase_shift_meson_density_scan_208_220_step2_v2.csv`

## 1. 当前 strict BW Stage1 口径

本页对应的 strict BW 仍是 **reduced Stage 1**，不是 full strict BW：

1. 上游 `mass/gamma` 来自 meson workflow 当前点输出；
2. 数密度积分采用
   - `E(q) = sqrt(q^2 + M^2)`
   - `Gamma(q) = Gamma(q=0)`
   - `omega = E(q) + Delta omega`
3. 当前还没有进入 `q` 依赖复极点 `z_p(q)` 的 full strict BW。

因此，这里记录的是：

- 一个正式挂在 workflow 上的 strict BW Stage 1 对照层；
- 不是最终 full strict BW 生产口径。

## 2. strict BW 扫描直接结果

`strict_bw_meson_density_scan_208_220_step2.csv` 当前给出：

- `gamma_pi`：约 `0.079 -> 0.457`
- `gamma_K`：约 `0.223 -> 0.651`
- `kpi_ratio`：约 `0.510 -> 0.588`

即：

1. `pi/K` 两个通道宽度在该窗口内单调增大；
2. `K/pi` 比值也随温度单调增大；
3. 当前 strict BW 输出字段、workflow 契约与 scan CSV 契约均已稳定落地。

## 3. 与 stable / Phase-E3 的并排比较

按同一温区重新跑 stable / strict BW / Phase-E3 后，当前最小比较结果为：

1. `stable`
   - `kpi_ratio`：约 `0.567 -> 0.638`
2. `strict BW Stage1`
   - `kpi_ratio`：约 `0.510 -> 0.588`
   - 相对 stable 的偏移：
     - `delta_kpi_bw`：约 `-0.057 -> -0.050`
3. `Phase-E3`
   - `kpi_ratio`：约 `1.369 -> 1.399`
   - 相对 stable 的偏移：
     - `delta_kpi_phase`：约 `+0.801 -> +0.761`

## 4. 当前最重要的观察

### 4.1 strict BW Stage1 没有贴近 stable

当前 reduced strict BW 的 `n_pi` / `n_K` 绝对值都明显低于 stable，因此 `K/pi` 也整体低于 stable。

这说明：

1. 当前 Stage 1 strict BW 不是一个“对 stable 的小修正”；
2. 它与此前 `BW-proxy` 给出的“略高于 stable”的结果并不等价；
3. `BW-proxy` 与 reduced strict BW 必须继续严格区分，不应混称为同一个 BW 结果。

### 4.2 strict BW Stage1 与 Phase-E3 的差异仍然很大

当前 `Phase-E3` 的 `K/pi` 仍在 `1.37-1.40`，而 strict BW Stage1 只有 `0.51-0.59`。

这意味着：

1. strict BW Stage1 不能替代完整相移 / BU 主线；
2. 它的职责仍然是“极点近似对照层”，不是主生产口径；
3. 现阶段仍不应把 strict BW Stage1 视为对 `Phase-E3` 的数值收口。

## 5. 当前解释边界

对这组结果，当前只宜作工程与口径层解释：

1. reduced strict BW 已经是一个正式、可回归、可扫描的 workflow 层；
2. 但它仍包含 `Gamma(q)=Gamma(0)` 与 `omega = E(q) + Delta omega` 这两个 Stage 1 降阶假设；
3. 因此它更适合作为：
   - strict BW 的第一版正式实现
   - full strict BW（`q` 依赖复极点）之前的过渡对照层

当前不宜把这组结果直接上升为：

1. 最终 strict BW 物理解读；
2. 或对 BU 主线的最终数值判据。

## 6. 下一步建议

基于本页结果，当前更合理的后续顺序是：

1. 先提交 Stage 1 reduced strict BW 这一轮改动；
2. 更新 PR draft，明确：
   - 新增 reduced strict BW workflow
   - `BW-proxy` 仍只是旧诊断层
   - full strict BW 尚未实现
3. 后续再决定是否进入 Stage 2：
   - 在 `q` 网格上直接求 `z_p(q)=E(q)-iGamma(q)/2`
   - 形成真正的 full strict BW 对照层
