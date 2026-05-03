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

本页对应的 strict BW 已进入 `q` 依赖复极点首轮实现，并采用严格的内层 `\omega` 积分口径：

1. workflow 先给出当前温点的 `pi/K` 介子结果；
2. strict BW 后处理在 `q` 网格上逐点重解复极点；
3. 将得到的 `E_M(q)`、`Gamma_M(q)` 回填到 BW 双积分；
4. 内层按 `\omega \in [0,\omega_{\max}]` 直接积分，而不是只取 `\omega \ge E_M(q)` 的右半边。

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

1. `n_pi`：约 `0.210 -> 0.604`
2. `n_K`：约 `0.583 -> 1.035`
3. `K/pi`：约 `2.78 -> 1.71`
4. `gamma_pi`：约 `0.073 -> 0.457`
5. `gamma_K`：约 `0.212 -> 0.650`

当前最直接的数值特征是：

1. `n_pi(T)` 随温度单调上升；
2. `n_K(T)` 处于更高量级，并在窗口内有起伏但保持正值；
3. `K/pi(T)` 整体大于 `1`，并在窗口内总体下降。

## 3. 与 Stage1 的最小对比

相比 Stage1 reduced strict BW，本轮严格 `\omega` 口径的 Stage2 结果明显抬高了 `K/pi`：

1. Stage1 `K/pi`：约 `0.510 -> 0.588`
2. Stage2 `K/pi`：约 `2.78 -> 1.71`

这说明：

1. “是否按严格 `\omega \in [0,\omega_{\max}]` 积分”不是小修正，而是 strict BW 是否闭环的关键差别；
2. reduced Stage1 与 q-pole Stage2 不能混写成同一个 strict BW 数值结果；
3. 后续若论文或正式图表需要引用 strict BW，必须显式标注是 `Stage1 reduced` 还是 `Stage2 q-pole`。

## 4. 当前工程结论

本轮更重要的是工程收口，而不是立即做物理解读：

1. Stage2 已正式挂在 workflow 主链上，而不是脚本层自行拼装；
2. 扫描脚本、测试、API 文档、CSV 契约已随 Stage2 一并打通；
3. `pole_require_converged=true` 的正式口径明显比 `best-effort` 更平滑，应作为当前默认扫描口径；
4. Stage2 现在已经可以作为 strict BW 的正式对照层，而不再只是草稿算法。
5. 当前这版已经补上了“严格 `\omega` 积分口径”这一闭环条件，strict BW 这一层可视为阶段性闭环完成。

## 5. 当前解释边界

这组结果当前只应作以下层面的解释：

1. 它是当前传播子口径下、已按严格 `\omega` 口径闭环的 `q-pole strict BW` 结果；
2. 它已经比 Stage1 更接近“由传播子复极点定义出发的 strict BW”；
3. 但它仍不是 generalized BU，也不能替代 full phase-shift / BU 主线。

当前不宜据此直接下结论：

1. Stage2 strict BW 必须接近 stable；
2. Stage2 strict BW 必须接近 Phase-E3；
3. Stage2 strict BW 已足够作为最终物理产出标准。

## 6. 下一步建议

基于当前状态，更合理的后续顺序是：

1. 将 strict BW 视为已阶段性闭环，后续论文若需要 stable / strict BW / BU 并排对照，当前实现已经可直接复用；
2. `BW-proxy` 继续只保留为旧诊断层；
3. 主精力转回更严格 BU / generalized BU 收口，而不是继续返工 strict BW 同一层积分定义。
