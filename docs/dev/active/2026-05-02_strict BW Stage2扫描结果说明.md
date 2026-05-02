# strict BW Stage2 扫描结果说明

更新日期：2026-05-02

目的：记录 `q` 依赖复极点版本 strict BW（Stage 2）在 `T = 208:2:220 MeV`、`muB = 0`、`xi = 0` 窗口下的首轮正式 workflow 扫描结果，并明确当前可接受的解释边界。

相关输出：

- Stage2 strict BW 扫描：
  - `data/outputs/results/relaxtime/scan/strict_bw_meson_density_scan_stage2_208_220_step2_converged.csv`
- Stage1 strict BW 扫描：
  - `data/outputs/results/relaxtime/scan/strict_bw_meson_density_scan_208_220_step2.csv`
- stable 扫描：
  - `data/outputs/results/relaxtime/scan/meson_density_scan_208_220_step2.csv`
- Phase-E3 扫描：
  - `data/outputs/results/relaxtime/scan/phase_shift_meson_density_scan_208_220_step2_v2.csv`

## 1. 本页对应的 Stage2 口径

本页对应的 strict BW 已进入 `q` 依赖复极点首轮实现：

1. workflow 先给出当前温点的 `pi/K` 介子结果；
2. strict BW 后处理在 `q` 网格上逐点重解复极点；
3. 将得到的 `E_M(q)`、`Gamma_M(q)` 回填到 BW 双积分。

本轮正式扫描采用：

1. `qmax = 12`
2. `q_nodes = 8`
3. `omega_max = 10`
4. `omega_nodes = 16`
5. `solver_iterations = 12`
6. `pole_residual_norm_max = 1e-4`
7. `pole_require_converged = true`

## 2. 直接结果

在 `208:2:220 MeV` 窗口内，当前 Stage2 扫描给出：

1. `n_pi`：约 `0.0732 -> 0.0634`
2. `n_K`：约 `0.0216 -> 0.0287`
3. `K/pi`：约 `0.296 -> 0.452`
4. `gamma_pi`：约 `0.073 -> 0.457`
5. `gamma_K`：约 `0.212 -> 0.650`

当前最直接的数值特征是：

1. `n_pi(T)` 随温度下降；
2. `n_K(T)` 随温度缓慢上升；
3. `K/pi(T)` 单调上升。

## 3. 与 Stage1 的最小对比

相比 Stage1 reduced strict BW，本轮 Stage2 首轮结果更偏向压低 `K/pi`：

1. Stage1 `K/pi`：约 `0.510 -> 0.588`
2. Stage2 `K/pi`：约 `0.296 -> 0.452`

这说明：

1. `q` 依赖复极点信息不是一个小修正；
2. reduced Stage1 与 q-pole Stage2 不能混写成同一个 strict BW 数值结果；
3. 后续若论文或正式图表需要引用 strict BW，必须显式标注是 `Stage1 reduced` 还是 `Stage2 q-pole`。

## 4. 当前工程结论

本轮更重要的是工程收口，而不是立即做物理解读：

1. Stage2 已正式挂在 workflow 主链上，而不是脚本层自行拼装；
2. 扫描脚本、测试、API 文档、CSV 契约已随 Stage2 一并打通；
3. `pole_require_converged=true` 的正式口径明显比 `best-effort` 更平滑，应作为当前默认扫描口径；
4. Stage2 现在已经可以作为 strict BW 的正式对照层，而不再只是草稿算法。

## 5. 当前解释边界

这组结果当前只应作以下层面的解释：

1. 它是当前传播子口径下的 `q-pole strict BW` 首轮实现结果；
2. 它已经比 Stage1 更接近“由传播子复极点定义出发的 strict BW”；
3. 但它仍不是 generalized BU，也不能替代 full phase-shift / BU 主线。

当前不宜据此直接下结论：

1. Stage2 strict BW 必须接近 stable；
2. Stage2 strict BW 必须接近 Phase-E3；
3. Stage2 strict BW 已足够作为最终物理产出标准。

## 6. 下一步建议

基于当前状态，更合理的后续顺序是：

1. 将 Stage2 这一轮 workflow/scan/doc 改动独立提交；
2. 在 PR draft 中明确：
   - Stage1 reduced strict BW 已稳定
   - Stage2 q-pole strict BW 已进入首轮正式实现
   - `BW-proxy` 继续只保留为旧诊断层
3. 继续把主精力放回更严格 BU / generalized BU 收口，而不是要求 strict BW 与主生产口径数值贴合。
