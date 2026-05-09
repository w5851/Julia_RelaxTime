# T200 双异常窗口（`tau_u`）证据链分析

## 1) 现象与范围

- 目标：针对 `T=200, muB=0` 的 `tau_u` 双异常窗口做与 T190 同构深拆，并回答“机制是否同构、虚部是否为 0”。
- 观测图：`data/outputs/figures/relaxtime/fixed_temperature_xi_scan_muB0/T200/multi_y_tau_u_tau_ubar_tau_s_tau_sbar_vs_xi.png`。
- 本文聚焦两个窗口：
  1. `xi∈[-0.4,-0.2]`
  2. `xi∈[0.2,0.4]`
- 主数据口径：`D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_main.csv`（本轮按固定参数重跑得到）。

## 2) 观测层：`tau_u` 与 `sigma_over_T` 的分段异常

- `[-0.4,-0.2]`：出现“先缓降 -> 深谷 -> 回升”
  - `tau_u: 0.6346 (-0.4) -> 0.3320 (-0.32) -> 0.2393 (-0.30) -> 0.5489 (-0.20)`。
  - `sigma_over_T` 同步呈深谷后回升：`0.01327 -> 0.00678 -> 0.00514 -> 0.00770`。
- `[0.2,0.4]`：出现“平台 + 局部凹陷 + 反弹”
  - `tau_u` 在 `0.34->0.36` 下跌（`1.8304 -> 1.4669`），随后 `0.36->0.38` 反弹（`1.4669 -> 2.0483`）。
  - `sigma_over_T` 同步 `0.34->0.36` 下行（`0.01264 -> 0.01026`），`0.36->0.38` 上行（`0.01026 -> 0.01346`）。

## 3) 串行链路证据（propagator -> sigma -> w_ij -> tauinv -> tau）

- 产物：
  - `D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_full_chain_detail.csv`
  - `D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_full_chain_summary.csv`
  - `D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_adjacent_transition_summary.csv`
- 关键相邻段（`uubar_to_ddbar`，面积比右/左）：
  1. 负窗口主下跌段 `-0.34 -> -0.32`：
     - `tau_u ratio ≈ 0.636`
     - `sigma_over_T ratio ≈ 0.667`
     - `ratio_area_abs_Dmixed_sq ≈ 1.017`
     - `ratio_area_abs_detM_sq ≈ 0.982`
  2. 负窗口主回升段 `-0.22 -> -0.20`：
     - `tau_u ratio ≈ 1.765`
     - `sigma_over_T ratio ≈ 1.633`
     - `ratio_area_abs_Dmixed_sq ≈ 1.047`
     - `ratio_area_abs_detM_sq ≈ 0.948`
  3. 正窗口主下跌段 `0.34 -> 0.36`：
     - `tau_u ratio ≈ 0.898`
     - `sigma_over_T ratio ≈ 0.894`
     - `ratio_area_abs_Dmixed_sq ≈ 1.165`
     - `ratio_area_abs_detM_sq ≈ 1.049`
  4. 正窗口主回升段 `0.36 -> 0.38`：
     - `tau_u ratio ≈ 1.247`
     - `sigma_over_T ratio ≈ 1.191`
     - `ratio_area_abs_Dmixed_sq ≈ 1.332`
     - `ratio_area_abs_detM_sq ≈ 1.052`

## 4) 与 T190 的机制对照：同链路，不同主导子机制

- T190 主异常段（`-0.10->-0.08`）的核心是 `|detM|` 抬升导致 `|D_mixed|` 塌缩（`<1`）。
- T200 双窗口里，`|D_mixed|` 在关键相邻段基本都 `>1`，且 `|detM|` 并未出现 T190 那种同段突升+放大量级。
- 因此：
  - **同构层面**：仍是同一串行链路（上游传播子变化传到 `tau`/输运）。
  - **差异层面**：T200 不是 T190 那种“单一 detM 分母增强触发的 Dmixed 塌缩型毛刺”，而是更偏向“多段重分配 + 通道权重切换”。

## 5) 通道贡献层（`tauinv_u`）主导项：两窗口都由 `uubar_to_ddbar / uubar_to_uubar / udbar_to_udbar` 主导切换

- 产物：`D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_channel_diag.csv`
- 负窗口主下跌段 `-0.34 -> -0.32`（`tauinv_u: 1.850 -> 3.012`）：
  - `udbar_to_udbar` 贡献增量 `+0.863`
  - `uubar_to_ddbar` 贡献增量 `+0.171`
  - `uubar_to_uubar` 贡献增量 `+0.130`
- 负窗口主回升段 `-0.22 -> -0.20`（`tauinv_u: 3.781 -> 1.822`）：
  - `udbar_to_udbar` 贡献减量 `-1.539`
  - `uubar_to_ddbar` 贡献减量 `-0.250`
  - `uubar_to_uubar` 贡献减量 `-0.209`
  - `usbar_to_usbar` 仅小幅补偿 `+0.043`
- 正窗口凹陷与反弹：
  - `0.34 -> 0.36`（`tauinv_u: 0.546 -> 0.682`）主要由
    - `uubar_to_ddbar +0.0736`
    - `uubar_to_uubar +0.0665`
  - `0.36 -> 0.38`（`tauinv_u: 0.682 -> 0.488`）主要由
    - `uubar_to_ddbar -0.1036`
    - `uubar_to_uubar -0.0854`

## 6) 虚部路径核验（公式条件逐条落位）

- 产物：
  - `D:\Desktop\Temp\relaxtime_t200_window\t200_imag_path_evidence_detail.csv`
  - `D:\Desktop\Temp\relaxtime_t200_window\t200_imag_path_evidence_summary.csv`
- 脚本：`scripts/analysis/relaxtime/t200_imag_path_evidence.jl`
- 结论：
  1. 采样点全部是 `s` 道口径，`k_s=0`（`k_min=k_max=0`），故本轮走 `k=0` 虚部分支。
  2. `simple` 分支（`usbar_to_usbar`）在各采样点 `pole_fraction=0.5`，`den_im` 非零（`~1e-6` 到 `~2e-5`）。
  3. `mixed_uu` 分支 `pole_fraction=0.5`，`pi_im_abs_max` 非零；`detM_im` 非零（`~1e-9` 到 `~1e-8`）。
  4. `mixed_ss` 分支 `pole_fraction=0` 且 `pi_im_abs_max=0`（严格为 0）。
- 口径：
  - “虚部完全为 0 的可能性”在 T200 同样存在（`Π_ss^P` 已见到严格 0）；
  - 但最终分母虚部（`den_simple_im` / `detM_im`）在关键异常段仍非零，不支持“异常段分母虚部全为 0”。

## 7) 可复现命令（本轮）

1. 双窗口主数据+通道贡献：
```bash
julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --output "D:/Desktop/Temp/relaxtime_t200_window/t200_dual_window_main.csv" --channel-diagnostics-output "D:/Desktop/Temp/relaxtime_t200_window/t200_dual_window_channel_diag.csv" --overwrite --no-resume --tmin 200 --tmax 200 --tstep 10 --mubmin 0 --mubmax 0 --mubstep 60 --xi-list -0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4 --mode finite_15 --tau-n-sigma 6 --sigma-grid-n 128
```
2. 全链路分解（双窗口）：
```bash
julia --project=. scripts/analysis/relaxtime/t200_dual_window_full_chain_decomposition.jl
```
3. 虚部路径核验：
```bash
julia --project=. scripts/analysis/relaxtime/t200_imag_path_evidence.jl
```

## 8) 当前结论与后续

- 当前可确认：
  1. T200 的 `tau_u` 双异常区并非 T190 单异常区的简单平移复刻；
  2. 同一串行链路成立，但主导结构是多段通道重分配（尤其 `udbar_to_udbar` 与 `uubar` 两主通道）而非单一 `detM` 增强导致 `Dmixed` 塌缩；
  3. 虚部在异常区不是“全为 0”，且“可 0 与可非 0”两类状态并存（取决于分支与极点条件）。
- 建议下一步（若继续加深）：
  1. 对 `udbar_to_udbar` 复制 T190 的 `sigma(s)` 阈值邻域分项反事实（K/blocking/t 区间）专项拆解；
  2. 在两窗口各选一个“下跌段+回升段”做 `M_s/M_t/干涉` 全分解对照表，以形成与 T190 结构完全对齐的章法。

## 9) 回答关键追问：三个主通道是否都由“分母实部先近 0 再远离 0”主导

- 追问对象：`udbar_to_udbar` 与 `uubar_to_ddbar/uubar_to_uubar` 在 T200 两异常窗口。
- 新增追踪脚本：`scripts/analysis/relaxtime/t200_target_channels_denominator_trace.jl`
- 产物：`D:\Desktop\Temp\relaxtime_t200_window\t200_target_channels_denominator_trace.csv`

- 结论先行：
  1. **`udbar_to_udbar`**：在负窗口主异常段，其主导变化主要由 `s` 道 simple 分母实部 `Re(1-4KΠ)` 过零触发；
  2. **`uubar_to_ddbar/uubar_to_uubar`**：在负窗口部分段与上条同向（simple 分母过零参与），但在正窗口主异常段更关键的是 mixed 分支 `detM` 变小（接近零）导致 `|D_mixed|` 放大；
  3. 因此“三通道两窗口”不是单一同号机制，而是“simple 过零 + mixed 近零”在不同子区间权重切换。

- 逐通道证据：
  1. `udbar_to_udbar`
     - `s` 道 simple 分母实部：`Re(den_simple)` 从 `+0.00699 (xi=-0.34)` 到 `+0.00273 (-0.32)` 再到 `-0.00204 (-0.30)`，在负窗口中段附近过零；
     - 同期 `|D_simple|^2` 从 `1.83e4` 跃升到 `1.20e5` 再到 `2.15e5`，与 `tau_u` 谷值形成同向；
     - `t` 道 mixed 的 `Re(detM)` 约 `0.00366~0.00378`，未见近零，不是主导。
  2. `uubar_to_ddbar`
     - `s` 道 simple 分母与上类似在负窗口中段过零（`+0.00699 -> +0.00273 -> -0.00204`），对负窗口异常段有贡献；
     - 但正窗口 `s` 道 mixed 的 `Re(detM)` 从 `7.94e-4 (xi=0.2)` 下降到 `3.98e-5 (xi=0.4)`，明显逼近零；
     - 对应 `|D_mixed|^2` 在 `0.36->0.38`、`0.38->0.4` 继续放大，解释正窗口后段反弹。
  3. `uubar_to_uubar`
     - 与 `uubar_to_ddbar` 同型：负窗口含 simple 分母过零，正窗口 mixed `detM` 近零增强更显著；
     - 该通道在 `tau_u` 贡献里与 `uubar_to_ddbar` 同为主导切换项（见第 5 节）。

- 小结口径（直接回答问题）：
  - “是否也都是分母实部先近0再远离0”——**部分成立但不完全**：
    - 在负窗口，中段确实可见 simple 分母过零链路；
    - 在正窗口，更强的是 mixed `detM` 接近零并持续增强，不是同一 simple 过零模式的镜像复现。

## 10) T200 的 `tau_s` 异常段（`xi∈[-0.3,-0.1]`）同构研究

- 数据与诊断：
  - 主数据：`D:\Desktop\Temp\relaxtime_t200_window\t200_taus_window_main.csv`
  - 通道贡献：`D:\Desktop\Temp\relaxtime_t200_window\t200_taus_window_channel_diag.csv`
  - 分母追踪：`scripts/analysis/relaxtime/t200_taus_usbar_denominator_trace.jl`
  - 产物：`D:\Desktop\Temp\relaxtime_t200_window\t200_taus_usbar_denominator_trace.csv`

- 观测层：
  - `tau_s` 在 `[-0.3,-0.1]` 呈显著“先降后升”：
    - `1.23 (-0.30) -> 0.448 (-0.22) -> 0.417 (-0.20) -> 1.73 (-0.10)`；
  - 对应 `tauinv_s`：`0.81 -> 2.23 -> 2.40 -> 0.58`，在 `-0.22/-0.20` 附近形成峰值。

- 通道层：
  - `tauinv_s` 变化几乎由 `usbar_to_usbar` 主导：
    - `-0.24 -> -0.22`: `+0.939`
    - `-0.22 -> -0.20`: `+0.169`
    - `-0.20 -> -0.18`: `-0.840`
  - 其他 `s` 相关通道变化量都仅 `O(1e-3)`。

- 分母同构性（与 T190/T200 其他段统一）：
  - `usbar_to_usbar` 的 `s` 道 simple 分母
    `den_simple = 1 - 4K_{4567}^+ Π_{us}^P` 的实部在该窗口发生过零：
    - `Re(den_simple) = +0.0193 (xi=-0.22)`
    - `Re(den_simple) = -0.00112 (xi=-0.20)`
  - 同点 `|D_simple|^2` 从 `1.81e3` 激增到 `5.49e5`，随后在 `xi=-0.18` 回落到 `1.29e3`；
  - `t` 道 mixed 的 `detM`（`~3.68e-3`）稳定且远离零，贡献次级。

- 结论：
  1. T200 的 `tau_s` 异常段与此前 simple 分支异常**同构**：核心是 `1-4KΠ` 的实部过零导致传播子增强；
  2. 该段不是 mixed `detM` 近零主导，而是 `usbar_to_usbar` 的 `s` 道 simple 主导；
  3. 因而可与 T190 的 `tau_s` 机制在数学结构上并轨（小分母非线性放大），只是窗口中心位置与幅值不同。

## 11) T190 vs T200 机制对照（一页表）

- 对照口径：统一按
  `propagator/amplitude -> sigma(s) -> w_ij -> tauinv -> tau -> transport`。
- 数据来源：
  - T190：`D:\Desktop\Temp\relaxtime_t190_window\t190_xi_window_adjacent_transition_summary.csv`、`t190_imag_path_evidence_summary.csv`
  - T200：`D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_adjacent_transition_summary.csv`、`t200_imag_path_evidence_summary.csv`

| 对象 | 异常窗口 | 主导通道（`tauinv` 贡献） | 主导分母对象 | 过零/近零证据 | 虚部路径（本口径） | 机制判定 |
|---|---|---|---|---|---|---|
| T190 `tau_u` 主异常 | `xi≈[-0.12,-0.08]`（核心 `-0.10->-0.08`） | `uubar_to_ddbar`、`uubar_to_uubar` | mixed `detM^P` | `ratio_area_abs_detM_sq≈1.126` 同段 `ratio_area_abs_Dmixed_sq≈0.349`（分母抬升、传播子塌缩） | `k=0`；`mixed_uu` 有极点项（`pole_fraction=0.5`），`mixed_ss` 可为 0；`detM_im` 非零小量 | **detM 分母增强主导** |
| T190 `tau_s` 异常 | `xi≈[-0.46,-0.40]` | `usbar_to_usbar` | simple `1-4K_{4567}^+\Pi_{us}^P` | `Re(den_simple)` 在窗口内过零；`|D_simple|^2` 阶跃放大 | `k=0`；simple 有极点项（`pole_fraction=0.5`），`den_im` 非零小量 | **simple 小分母放大主导** |
| T200 `tau_u` 负窗口 | `xi∈[-0.4,-0.2]`（核心 `-0.34->-0.30`） | `udbar_to_udbar` 为首，`uubar_to_ddbar/uubar_to_uubar` 次之 | simple 为主，mixed 为辅 | `udbar` 的 `Re(den_simple)`：`+0.00699(-0.34) -> +0.00273(-0.32) -> -0.00204(-0.30)`；对应 `|D_simple|^2` 放大；同期 `t` 道 mixed `Re(detM)~3.7e-3` 远离 0 | `k=0`；simple 与 `mixed_uu` 均有极点项，`mixed_ss` 可为 0 | **simple 过零主导（局部 mixed 参与）** |
| T200 `tau_u` 正窗口 | `xi∈[0.2,0.4]`（核心 `0.34->0.40`） | `uubar_to_ddbar` 与 `uubar_to_uubar` 双主导 | mixed `detM` | `Re(detM)` 持续逼近 0：`~7.94e-4(0.2) -> 3.98e-5(0.4)`；`ratio_area_abs_Dmixed_sq` 在 `0.34->0.36/0.36->0.38/0.38->0.4` 依次 `1.165/1.332/3.024` | `k=0`；`mixed_uu` 有极点项且 `detM_im` 非零小量；`mixed_ss` 可为 0 | **detM 近零增强主导** |
| T200 `tau_s` 异常 | `xi∈[-0.3,-0.1]`（核心 `-0.22/-0.20`） | `usbar_to_usbar` 绝对主导 | simple `1-4K_{4567}^+\Pi_{us}^P` | `Re(den_simple): +0.0193(-0.22) -> -0.00112(-0.20)`；`|D_simple|^2: 1.81e3 -> 5.49e5 -> 1.29e3` | `k=0`；simple 极点项在，`den_im` 非零小量 | **simple 过零主导** |

- 统一结论：
  1. T190 与 T200 **同链路**，但 `tau_u` 的窗口主导对象不同：T190 以 mixed `detM` 分母增强/塌缩为主，T200 呈“负窗口 simple 过零 + 正窗口 mixed 近零”的分段切换；
  2. `tau_s` 在两温度点均可落到 simple 小分母同构机制（`1-4KΠ` 过零触发）；
  3. 虚部在异常段不是“全为 0”：本口径均为 `k=0` 分支，且呈现“部分子分支严格 0 + 最终分母虚部非零小量”并存。

## 12) 严格版总命题 + 正窗口反例条款（按本轮复核）

- 严格版总命题（可被当前证据支持）：
  - `T190/T200` 已识别的 `u/s` 异常窗口可统一归入：
    **某些散射通道的某类传播子分母在局部成为小量（过零或近零）且伴随小虚部**，从而在 `1/|den|^2` 结构中触发非线性放大/塌缩。

- 不应再使用的强命题：
  - “所有窗口都由 `Re(den)` 线性穿零 + `Im(den)` 非零小量触发”。

- 反例条款（T200 `tau_u` 正窗口）：
  1. 在主异常段（约 `xi=0.34->0.36` 下跌、`0.36->0.38` 回升），主导通道 `uubar_to_ddbar` 的 `s` 道 mixed 分母 `detM` 满足：
     - `Re(detM): 2.3707e-4 (0.34) -> 1.6891e-4 (0.36) -> 1.0319e-4 (0.38)`，同号单调减小；
     - `Im(detM): -8.03e-8 -> -8.17e-8 -> -8.30e-8`，为非零小量；
     - `|detM|` 对应 `2.37e-4 -> 1.69e-4 -> 1.03e-4`，即“更近零”而非穿零。
  2. 对应传播子增强：
     - `|D_mixed|^2: 1.3177e2 (0.34) -> 2.5733e2 (0.36) -> 6.8377e2 (0.38)`，连续放大；
     - 相邻段面积比（`uubar_to_ddbar`）：`ratio_area_abs_Dmixed_sq(0.34->0.36)=1.165`，`(0.36->0.38)=1.332`。
  3. 与 `tau_u` 层联动（同段）：
     - `tau_u: 1.8304 (0.34) -> 1.4669 (0.36) -> 2.0483 (0.38)`；
     - `tauinv_u: 0.5463 -> 0.6817 -> 0.4882`；
     - `uubar_to_ddbar` 贡献：`0.1306 -> 0.2042 -> 0.1005`。

- 对“近0后远离0”这句的判定（仅针对你指定主段）：
  - 在 `0.34->0.38` 这一主段里，证据显示的是“持续近0”，**尚未出现远离0**；
  - “远离0”发生在离开该主段之后（未在你指定的两步主波动中出现）。

## 13) 为什么主通道贡献会“先升后降”：rate 核分解

- 你指出的问题成立：
  - `uubar_to_ddbar` 的 `s` 道 mixed 分母在 `0.34->0.36->0.38` 持续近零；
  - 但该通道贡献在诊断表中呈现“先升后降”。

- 这不矛盾，原因是：
  - 通道贡献定义是 `contribution = multiplicity * density * rate`；
  - 其中 `rate` 是 5D 热平均核的结果（含 `f_i f_j`, `v_rel`, 全 `s` 区间插值 `σ(s)` 与截断），不是“单个近阈值点的传播子值”。

- 新增脚本与产物：
  - 脚本：`scripts/analysis/relaxtime/t200_rate_kernel_decomposition_034_038.jl`
  - 产物：`D:\Desktop\Temp\relaxtime_t200_window\t200_rate_kernel_decomposition_034_038.csv`

- 核分解结果（`xi=0.34/0.36/0.38`）：
  1. `rate_func`（同参数口径重算）为 `0.6606 -> 0.6719 -> 0.6166`，在 `0.36->0.38` 明确回落；
  2. 该回落主要来自热平均核中的 `FV_mean_sigma` 下行：
     - `0.4881 -> 0.4980 -> 0.4586`；
  3. 近阈值极窄区（`Δs<=0.02`）在总核中的份额仅 `~7e-4`，`Δs<=0.10` 份额也仅 `~5e-3`，说明最终 rate 不是由“最靠近阈值的局部点”单独控制；
  4. 因而即使某局部传播子持续增强，热平均后 `rate` 仍可能因为权重重排而下降。

- 结论口径修正：
  - “分母近零”是该通道局部增强的**必要条件**；
  - 但“该通道贡献单调随之增强”不是充分条件；
  - 充分决定量是热平均核 `rate`（再乘密度与重数）。

## 14) 在线埋点复核：穿零窗口是否足以主导 + 正窗口“下一层”数值链条

- 你要求“同路径、可审计”是正确的。本节改为在 `AverageScatteringRate._omega_integral_5d` 内加可开关 band 埋点（不是离线重算另一路），并由分析脚本读取同路径输出。

- 代码与脚本：
  1. `src/relaxtime/AverageScatteringRate.jl`
     - 新增可选参数：`band_edges`, `band_omega_out`, `band_omega_sigma_out`；
     - 在 5D 主积分循环内按 `Δs=s-s_th` 分段累积：
       - `omega_bin = Σ (w * f_i * f_j * v_rel)`
       - `omega_sigma_bin = Σ (w * f_i * f_j * v_rel * σ)`
  2. `scripts/analysis/relaxtime/t200_rate_band_window_checks.jl`
     - 三组窗口同路径导出并比较：
       - `tau_u` 负窗：`udbar_to_udbar`, `xi=-0.34/-0.32/-0.30`
       - `tau_s` 异常窗：`usbar_to_usbar`, `xi=-0.22/-0.20/-0.18`
       - `tau_u` 正窗主段：`uubar_to_ddbar`, `xi=0.34/0.36/0.38`
  3. 输出：
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_rate_band_window_checks_detail.csv`
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_rate_band_window_checks_pair_delta.csv`

- 问题 1：穿零窗口是否“足以主导”
  1. 先给审计修正：`scripts/analysis/relaxtime/t200_rate_band_window_checks.jl` 已与扫描主路径完全对齐（读取主 CSV 的 `tau_*` 元数据、复用 workflow 风格 `A` 构造并保持 `K_coeffs` 口径一致），`rate_func` 与 `channel_diag` 的 `rate` 现已基本一致：
     - `tau_s/usbar` 最大绝对差约 `1.08e-4`；
     - `tau_u/udbar` 最大绝对差约 `1.23e-2`；
     - `tau_u/uubar_to_ddbar` 最大绝对差约 `4.60e-3`。
  2. 在这一对齐口径下，三组窗口都显示“主导增量主要来自中高 `Δs` 带（尤其 `2~20`）”，不是极近阈值窄带单独主导。
  3. 因此第一问的审计结论更新为：
     - **“分母近零/穿零”是局部增强触发条件，但在这批 T200 窗口中不是单独充分条件；最终方向由全带宽热平均核分配决定。**

- 问题 2：你要的“下一层”不是黑话——给出可计算分解
  1. rate 的同路径分解是：
     - `rate = prefactor * Σ_b omega_sigma_bin(b)`
     - `omega_sigma_bin(b) = omega_bin(b) * sigma_eff_bin(b)`
     - 其中
       - `omega_bin(b)=Σ(w*f_i*f_j*v_rel)`（纯权重核）
       - `sigma_eff_bin(b)=omega_sigma_bin/omega_bin`（该带有效截面）
  2. 所谓“更宽 s 带重排”，在本节里有明确数值定义：
     - 不是一句话，而是指 `Δomega_sigma_bin` 在中高 `Δs` 带（如 `Δs>0.2`）出现系统负增量，且绝对值超过近阈值带正增量。
  3. 以正窗主段 `0.36->0.38`（`uubar_to_ddbar`）为例（同路径在线埋点）：
     - 极近阈值带转负且幅值有限：
       - `Δs∈(0,0.02]`: `Δomega_sigma ≈ -8.23e-1`
       - `Δs∈(0.02,0.05]`: `≈ -3.36`
     - 中高带同样持续负增量并累计主导：
       - `Δs∈(2,5]`: `≈ -6.53e-1`
       - `Δs∈(5,10]`: `≈ -6.92e-1`
       - `Δs∈(10,20]`: `≈ -7.18e-1`
       - `Δs∈(20,∞]`: `≈ -4.17e-1`
     - 这给出更稳健的链条：不是“某一窄带决定一切”，而是多带协同下行（近阈值与中高带同向）。

- 本轮口径更新：
  1. 不再把“穿零/近零”当作 rate 单调性的充分条件；
  2. 对每个异常窗口坚持分层判断：
     - 层A：传播子分母小量机制（局部增强触发）；
     - 层B：`ωσ` 分段重排（决定 rate 与贡献最终方向）；
  3. 当前三窗口里，层B主导带均落在中高 `Δs`，不存在“仅靠极近阈值窄带即可解释全部方向”的证据。

- 对“下一层解释是否严格对应”的补充（回答复核问答）：
  1. 在 `tau_u` 正窗 `0.36->0.38`，中高 `Δs` 带负增量来自两层并行：
     - `ω_bin` 下降（热权重核变小）；
     - `σ_eff_bin` 下降（该带有效截面变小）。
  2. 典型带（`uubar_to_ddbar`）：
     - `Δs∈(2,5]`：`omega_bin 3.4165 -> 3.3078`，`sigma_eff 0.6316 -> 0.4550`，`ωσ` 负增量；
     - `Δs∈(5,10]`：`omega_bin 5.2238 -> 5.0623`，`sigma_eff 0.4392 -> 0.3166`，`ωσ` 负增量；
     - `Δs∈(10,20]`：`omega_bin 7.6026 -> 7.3898`，`sigma_eff 0.3155 -> 0.2274`，`ωσ` 负增量。
  3. 因此更准确的说法是：
     - 不是“单点传播子变强/变弱”一句话；
     - 而是“多 `Δs` 带上权重核与有效截面的耦合演化共同决定总符号”。

- 普遍性判定（避免过度泛化）：
  1. 该机制不是单一模板可覆盖，需窗口逐个审计；
  2. 当前最稳妥统一模板是“先判 A（分母小量），再判 B（分段 `ωσ` 主导带），并以同路径 rate 一致性作为前置验收”。

## 15) “tau 先减后增”的统一解释（含传播子层、全带宽）

- 先写成可判定形式：
  1. `tau = 1/tauinv`；
  2. 所以“tau 先减后增”等价于“tauinv 先增后减”；
  3. 而 `tauinv_i = Σ_j contribution_{i<-j}`，每一项由对应通道 `rate` 决定。

- 对应到你关心的 `tau_u` 正窗（`xi=0.34->0.36->0.38`）：
  1. `0.34->0.36`：`tauinv_u` 增加（约 `+0.1354`），所以 `tau_u` 从 `1.830` 降到 `1.467`；
  2. `0.36->0.38`：`tauinv_u` 降低（约 `-0.1935`），所以 `tau_u` 从 `1.467` 升到 `2.048`；
  3. 触发该“符号翻转”的主通道对是同一组：
     - `uubar_to_ddbar`：`+0.0736` 后转 `-0.1036`；
     - `uubar_to_uubar`：`+0.0665` 后转 `-0.0854`。

- 传播子层的“全带宽”解释（不局限近阈值）：
  1. 每个通道可写成
     - `σ(s,xi) ~ |N(s,xi)/D(s,xi)|^2 × K(s)`，
     - 其中 `D` 对应 mixed 分母（如 `detM`），`N` 对应投影/耦合组合（`JMJ` 及相关结构）；
  2. A 层（`D` 变小）只提供“可放大性”，但真正控制通道增量符号的是 `σ(s)` 在整个 `Δs` 带上的谱形变化（峰值位置、峰宽、尾部斜率、带间再分配）；
  3. 同路径在线分段结果显示，`0.34->0.36` 与 `0.36->0.38` 的符号翻转来自“多带协同翻转”：
     - `0.34->0.36`：主要带（含近阈值与中高带）`Δomega_sigma_bin` 以正为主；
     - `0.36->0.38`：同一批主要带转为负，且累计负增量更大。

- 这个“多带翻转”在数值上的直接证据（`uubar_to_ddbar`）：
  1. 近阈值带：
     - `sigma_eff(0.02,0.05]`：`8.81 -> 640.77 -> 6.89`（先暴涨后回落）；
  2. 中高带同样出现先升后降：
     - `sigma_eff(2,5]`：`0.5719 -> 0.6316 -> 0.4550`；
     - `sigma_eff(5,10]`：`0.3979 -> 0.4392 -> 0.3166`；
     - `sigma_eff(10,20]`：`0.2858 -> 0.3155 -> 0.2274`；
  3. 同时这些中高带 `omega_bin` 也在 `0.36->0.38` 下降（例如 `(2,5]`: `3.416->3.308`, `(5,10]`: `5.224->5.062`, `(10,20]`: `7.603->7.390`），
     最终使 `omega_sigma_bin = omega_bin * sigma_eff_bin` 同向下降。

- 因而“主导通道增量符号翻转”的传播子层原因可归纳为：
  1. 不是仅由“分母是否近零”单变量决定；
  2. 而是 `N/D` 共同决定的 `σ(s)` 谱形在多 `Δs` 带上的再分配，叠加热权重核后发生净符号翻转；
  3. 这就是异常窗口里常见“tau 先减后增”的统一根因链。

## 16) “分母层 -> 传播子层 -> 截面层 -> rate层”的同带宽最终汇总

- 你要求的最终口径已改为“**每个 band 一行**、四层同带宽并表”，并且新增了“**不积分**的 `Δs` 离散采样”用于看分母取 `abs` 前后的趋势。

- 脚本与输出（同路径、可审计）：
  1. 脚本：`scripts/analysis/relaxtime/t200_denominator_fullband_scan.jl`
  2. 同带宽四层总表（最终汇总图表口径）：
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_fullchain_band_table.csv`
  3. 带宽积分层（保留）：
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_fullband_detail.csv`
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_fullband_pair_delta.csv`
  4. `Δs` 离散采样（不做积分）：
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_ds_samples.csv`

- 最终汇总表列定义（`t200_fullchain_band_table.csv`）：
  1. 分母层：`area_invabs_den_simple`, `area_invabs_detM`
  2. 传播子层：`area_abs_D_total_sq`
  3. 截面层：`area_sigma`, `sigma_eff_bin`
  4. rate 层：`omega_bin`, `omega_sigma_bin`, `rate_bin`（并保留 `rate_func/rate_direct/contrib_direct` 作为通道总量锚点）

- 正窗口主翻转段（`uubar_to_ddbar`, `0.36->0.38`）按 band 的四层同向性：
  1. `Δs∈(2,5]`：
     - 分母层 `delta_area_invabs_detM = -2.52e4`；
     - 传播子层 `delta_area_abs_D_total_sq = -7.77e-2`；
     - 截面层 `delta_area_sigma = -9.18e-3`；
     - rate层 `delta_omega_sigma_bin = -6.53e-1`，`delta_rate_bin = -3.50e-2`。
  2. `Δs∈(5,10]`：
     - `delta_area_invabs_detM = -3.22e4`；
     - `delta_area_abs_D_total_sq = +1.19e-3`（近零小回弹）；
     - `delta_area_sigma = +8.89e-4`（近零）；
     - 但 rate层仍为 `delta_omega_sigma_bin = -6.92e-1`，`delta_rate_bin = -3.71e-2`。
  3. `Δs∈(10,20]`：
     - `delta_area_invabs_detM = -3.19e4`；
     - `delta_area_abs_D_total_sq = -3.46e-3`；
     - `delta_area_sigma = +1.14e-3`（小量）；
     - rate层 `delta_omega_sigma_bin = -7.18e-1`，`delta_rate_bin = -3.85e-2`。

- `Δs` 离散采样（你关心的“分母取 abs 前趋势”）
  1. `uubar_to_ddbar, xi=0.36`：`detM_re(Δs)` 在采样点中发生符号翻转
     - `Δs=0.001: +1.19e-5`
     - `Δs=0.02: +2.84e-6`
     - `Δs=0.05: -1.14e-5`
     - `Δs=0.2: -8.11e-5`
     - `Δs=5: -1.25e-3`
     这说明在同一参数下，`detM` 的“取 abs 前核”随 `Δs` 从近零正值跨到负值并逐步远离零。
  2. `uubar_to_ddbar, xi=0.38`：`detM_re(Δs)` 全采样点同号为负且幅值增大
     - `Δs=0.001: -5.80e-5`
     - `Δs=0.1: -1.04e-4`
     - `Δs=1: -4.75e-4`
     - `Δs=5: -1.29e-3`
     即该参数下“取 abs 前”不存在再过零，而是随 `Δs` 继续向负侧展开。
  3. simple 分母 `den_simple` 也可直接看取 abs 前趋势：
     - 例如 `udbar_to_udbar, xi=-0.34`：`den_simple_re` 在 `Δs≈0.5` 前后由正转负（`0.5: +1.12e-2`, `1.0: -3.30e-2`）；
     - 对应 `invabs_den_simple` 在近零区峰化（`Δs=0.5` 约 `2.41e3`）后随远离阈值回落（`Δs=10` 约 `1.18`）。

- 对“被积函数分母核在取绝对值前随 `Δs` 的变化趋势”给出统一解释：
  1. 先看复分母本体：`den_re(Δs), den_im(Δs)` 或 `detM_re(Δs), detM_im(Δs)`；
  2. 再看 `1/|den(Δs)|^2` 的结果只是在该复轨迹上的模长映射；
  3. 因而“`rate` 翻转”不需要要求 `den_re` 在所有 `Δs` 都过零，而是由多个 band 上复分母轨迹变化共同决定。

- 本节最终结论（按你要求的根因层）：
  1. 近阈值的过零/近零只是局部触发条件；
  2. 真正决定 `tau` 异常段方向的是同带宽四层量在中高 `Δs` 带的合成符号（尤其 `omega_sigma_bin`/`rate_bin`）；
  3. `Δs` 离散采样已证明：在同一参数点，分母“取 abs 前核”的趋势可以与 `xi` 方向直觉不同，必须按 `Δs` 逐点看其复轨迹。

## 17) 两张图的读图结论（审核意见）

- 图文件：
  1. `docs/analysis/relaxtime/tauu_pos_uubarddbar_uubar_to_ddbar_denominator_chain.png`
  2. `docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_denominator_chain.png`
- 审核结论：**同意你的判断**。这两张图已经能够把异常链路讲清楚：
  1. 在“分母层”（前两幅子图），`detM_re(Δs)` 与 `detM_im(Δs)` 的轨迹变化直接决定 `|detM|` 的收缩/放大节奏；
  2. 在“传播子层”（第三幅子图），该轨迹映射为 `inv|detM|^2` 与 `|D_total|^2` 的数量级变化；
  3. 在“截面/rate层”（第四幅子图），再通过 band 加权汇总到 `omega_sigma_bin` 与 `rate_bin`，从而决定 `tauinv` 增减与 `tau` 的“先降后升”。
- 对“异常可由分母实部/虚部随 `Δs` 变化趋势解释”这句话，建议定稿为：
  - “分母复轨迹（`Re/Im` 随 `Δs`）给出局部增强触发机制，最终异常方向由全带宽加权后的 `omega_sigma_bin` 共同决定。”
  - 这样既保留你强调的核心，也避免把结论误写成“只看分母就足够”。
- 与你后续研究问题直接相关的观察：
  1. 在 `xi=0.36/0.38` 的采样中，`detM_im(Δs)` 明显出现“先更负、后回升并转正”的段落；
  2. 该“先减后增”对应图中高 `Δs` 端的拐点，是下一步值得单独建模的对象（例如拆 `Π` 的实虚部来源与耦合系数项）。
- 为便于面向非专业读者串联叙事，新增了一份“零基础读图版”说明：
  - `docs/analysis/relaxtime/t200-denominator-chain-for-readers.md`

## 18) 按你指定口径复核：正窗口关键是小 `Δs` 区 `detM_re` 的“接近 0 程度”

- 口径范围：`tauu_pos_uubaruubar / uubar_to_uubar`，`xi=0.34/0.36/0.38`，重点看 `Δs<=0.2`。
- 数据源：`D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_ds_samples.csv`。

- 先给结论（与你的判断一致）：
  1. 在小 `Δs` 区，三条曲线的 `detM_im` 都是同号非零小量且量级接近，不能解释主要差异；
  2. 主要差异来自 `detM_re` 的“离 0 的距离”：`xi=0.36` 在近阈值区最贴近 0；
  3. 这会直接压低 `|detM|^2 = detM_re^2 + detM_im^2`，从而放大 `inv|detM|^2`，在 `Δs` 小区间形成更强局部增强。

- 离散点证据（`Δs=0.001/0.005/0.01/0.02`）：
  1. `detM_im` 的跨 `xi` 差异很小（例如 `Δs=0.01`：`-7.86e-5 / -7.99e-5 / -8.10e-5`），
     说明虚部在这三点上“同量级、同趋势”；
  2. `detM_re` 的差异更关键：
     - `Δs=0.01`：`|detM_re|` 为 `8.01e-5 (xi=0.34)`、`7.60e-6 (xi=0.36)`、`6.22e-5 (xi=0.38)`；
     - `Δs=0.02`：`|detM_re|` 为 `7.53e-5`、`2.84e-6`、`6.69e-5`；
     其中 `xi=0.36` 明显最接近 0。

- 对放大量的直接影响（同一离散点）：
  1. `Δs=0.001`：`inv|detM|^2` 比值
     - `xi=0.36 / 0.34 ≈ 9.97`
     - `xi=0.36 / 0.38 ≈ 5.16`
  2. `Δs=0.005`：约 `3.01` 与 `2.09`；
  3. `Δs=0.01`：约 `1.96` 与 `1.62`；
  4. `Δs=0.02`：约 `1.41` 与 `1.38`。

- 物理解释口径（可直接用于正文）：
  - 在近阈值区，当 `detM_im` 已固定在“非零但相近”的小量带时，
    `detM_re` 是否更贴近 0 决定了 `|detM|` 的主导项结构；
  - `xi=0.36` 对应更小的 `|detM_re|`，因此在小 `Δs` 上得到更强的分母放大，
    这正是你强调的“关键不在过零先后，而在接近 0 程度”的定量证据。

## 19) 启动 `detM_im`“先减后增”机制研究（阶段结论）

- 在相同对象（`tauu_pos_uubaruubar / uubar_to_uubar`，`xi=0.34/0.36/0.38`）下，
  `detM_im(Δs)` 呈一致分段：
  1. `Δs: 0.001 -> 5.0`：单调变得更负（下降）；
  2. `Δs: 5.0 -> 10.0`：由负转正（回升并过零）；
  3. `Δs: 10.0 -> 20.0`：继续增大为正。

- 数值锚点（`xi=0.36` 为例）：
  - `detM_im(0.2)=-3.65e-4`
  - `detM_im(2.0)=-1.31e-3`
  - `detM_im(5.0)=-1.94e-3`（最负）
  - `detM_im(10.0)=+2.40e-4`（过零后）
  - `detM_im(15.0)=+1.18e-2`
  - `detM_im(20.0)=+3.09e-2`

- 跨 `xi` 的一致性：
  - 三个 `xi` 的最小值都在 `Δs=5.0` 附近，且数值非常接近（约 `-1.93e-3`）；
  - 过零点都出现在离散采样的 `Δs=10.0` 档位（说明该形态更像“共享结构驱动”，而非单点噪声）。

- 当前阶段判断：
  - `detM_im` 的“先减后增”是一个对 `xi=0.34/0.36/0.38` 近乎共形的 `Δs` 结构信号；
  - 因而下一步应优先拆 `detM` 虚部的构成项（`Π` 的虚部来源与系数组合），
    而不是先从 `xi` 差异入手。

## 20) 上游机制图（按你要求：只看 `detM_im` 上游量）

- 图文件：
  - `docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_detM_im_upstream_mechanism.png`
- 生成链路：
  1. 上游追踪脚本：`scripts/analysis/relaxtime/t200_detm_im_upstream_trace.jl`
  2. 追踪数据：`D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_trace.csv`
  3. 机制绘图脚本：`scripts/analysis/relaxtime/plot_t200_detm_im_upstream_mechanism.py`

- 三幅子图读图结论：
  1. `Im(Π_uu^P)` 随 `Δs` 单调抬升，且三条 `xi` 曲线几乎重合；
  2. `Im(Π_ss^P)` 在大部分区间近似 0，仅在高 `Δs` 端（`Δs=20` 采样点）出现明显抬升；
  3. `detM_im` 的符号翻转由两项分解竞争决定：
     - 低/中 `Δs`：`Im(M00*M88)<0` 且 `Im(M08^2)>0`，故
       `detM_im = Im(M00*M88) - Im(M08^2)` 更负（先减）；
     - `Δs≈10`：`Im(M00*M88)` 由负转正并接近 `Im(M08^2)`，出现过零；
     - 高 `Δs`：`Im(M08^2)` 转负，减去负项等效“加法”，与 `Im(M00*M88)>0` 同向叠加，
       使 `detM_im` 快速转正并增大（再增）。

- 关键锚点（`xi=0.36` 示例）：
  1. `Δs=5.0`：`Im(M00*M88)=-3.61e-4`，`Im(M08^2)=+1.65e-3`，`detM_im=-2.01e-3`（谷底附近）；
  2. `Δs=10.0`：`Im(M00*M88)=+1.77e-3`，`Im(M08^2)=+1.69e-3`，`detM_im=+8.72e-5`（过零）；
  3. `Δs=20.0`：`Im(M00*M88)=+1.84e-2`，`Im(M08^2)=-1.19e-2`，`detM_im=+3.03e-2`（正向放大）。

- 阶段结论（上游层）：
  - `detM_im` 的“先减后增”不是单一 `Π` 项的单调映射，而是
    `Im(M00*M88)` 与 `Im(M08^2)` 的符号/幅值重排结果；
  - 且该重排在 `xi=0.34/0.36/0.38` 上高度一致，支持“共享结构机制”判断。

## 21) 图线可见性与 `t` 口径澄清（你最新两点问题）

- 问题 A：为什么看起来有些线段“不完整”
  1. 该图中三条 `xi` 曲线在多处高度重合（尤其 `Im(Π_uu^P)` 与分解项中段），
     在像素级会表现为“上层曲线遮住下层曲线”，视觉上像断线；
  2. `Im(Π_ss^P)` 在多数采样点非常接近 0，若 y 轴范围偏宽，也会出现“贴轴看不清”的效果；
  3. 为减少该假象，绘图脚本已做三项改进：
     - 增大画布并外置图例；
     - 给不同 `xi` 使用不同 `linestyle/marker/markevery`；
     - 对 `Im(Π_ss^P)` 子图使用更聚焦的 y 轴范围。

- 问题 B：这些图里的 `t` 到底怎么处理
  1. 本文“分母/上游机制诊断图”（含本节机制图）是**定点切片诊断**：
     - 对每个采样 `s`，取
       `t_mid = 0.5 * (t_min(s) + t_max(s))`，
       然后在 `(s, t_mid)` 上评估 `Π`、`M`、`detM` 等上游量；
  2. 这类图的用途是解释“上游结构如何随 `Δs` 变化”，不是直接替代物理观测量；
  3. 真实的截面/输运率计算口径是：
     - 在 `total_cross_section` 中对 `t∈[t_min,t_max]` 做数值积分；
     - 在 `average_scattering_rate` 中再对动量与角度做热平均（5D 核）；
  4. 因此你可把机制图理解为“用于机理定位的 `t_mid` 代表切面”，
     最终 `tau` 结论仍以“全 `t` 积分 + 全相空间平均”的主路径结果为准。

- 口径结论（防误读）
  1. 机制图回答的是“上游复分母结构为何会产生 `detM_im` 的先减后增”；
  2. 物理量方向（`tau` 先降后升）仍需并且已经通过同路径积分量（`omega_sigma_bin`/`rate_bin`）闭环验证；
  3. 两者不是冲突关系，而是“机理切片 + 物理积分”的分层证据链。

- 若继续加深（下一步可选）
  1. 增加 `t`-敏感性附录：在固定 `s` 比较 `t=t_min/t_mid/t_max` 三点，
     对 `detM_im` 分解项给出偏差带；
  2. 再与 `t`-积分后的带宽结果并排，量化“切片解释与积分结论”的一致区间。

## 22) `t`-敏感性附录（`t_min/t_mid/t_max` 三点复核）

- 你要求的附录已按同路径补齐：
  1. 数据脚本：`scripts/analysis/relaxtime/t200_detm_im_t_sensitivity_trace.jl`
  2. 图脚本：`scripts/analysis/relaxtime/plot_t200_detm_im_t_sensitivity.py`
  3. 输出图：`docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_detM_im_t_sensitivity.png`
  4. 汇总表：`D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_t_sensitivity_summary.csv`

- 结果（本对象：`uubar_to_uubar`，`s` 道 mixed 上游）非常明确：
  1. 在全部采样点（`xi=0.34/0.36/0.38`，`Δs=1e-3..20`）中，
     `detM_im(t_min) = detM_im(t_mid) = detM_im(t_max)`；
  2. 汇总表中 `range_abs` 与 `range_rel_to_mid` 全部为 `0`；
  3. 即本附录下 `detM_im` 对 `t` 三点切片无可见敏感性。

- 解释（为何会出现“零敏感性”）：
  1. 当前追踪对象使用 `:s` 道口径，极化函数输入由 `k0,k` 决定；
  2. 在此实现路径下，`s` 已固定时该口径对 `t` 不显式展开，
     因而 `t_min/t_mid/t_max` 给出相同上游 `detM_im`；
  3. 这与“物理截面最终仍需全 `t` 积分”不冲突：
     - 上游 `detM` 诊断切片可对 `t` 不敏感；
     - 最终可观测量依然通过 `t`-积分与热平均闭环。

- 因而本问题可定稿为：
  - 对你当前关心的上游 `detM_im` 机制图（`uubar_to_uubar`, `:s` 通道），
    使用 `t_mid` 不会引入额外偏差；
  - 机制图中的“先减后增”不是 `t` 切片选择造成，而是分解项本身的 `Δs` 结构。

## 23) 你这个理论问题的代码级回答：`t` 是否会“间接影响”当前上游图

- 先给结论（针对当前这套图的对象与通道）：
  1. 你说的理论路径是成立的：一般来说 `t` 可以通过影响 `(k0, k)` 进而影响极化函数与传播子；
  2. 但在**当前分析对象**（`uubar_to_uubar` 的 `:s` 道 mixed 上游）里，这条间接路径在实现上退化为常量：
     - `k_s = 0`；
     - `k0_s = sqrt(s)`，与 `t` 无关；
  3. 所以该对象下 `Π`、`detM_im`、`inv|detM|^2`、`|D_mixed|^2` 对 `t` 扫描全部平坦。

- 代码依据（对应实现位置）：
  1. `src/relaxtime/TotalPropagator.jl` 的 `calculate_cms_momentum(...)` 中，
     `channel == :s` 分支显式写为 `k0 = abs(E1 + E2)`、`k = 0.0`；
  2. 同一函数中 `E1=(s+m1^2-m2^2)/(2*sqrt(s))`、`E2=(s-m1^2+m2^2)/(2*sqrt(s))`，
     因此 `E1+E2 = sqrt(s)`，不含 `t`；
  3. 当前上游追踪脚本调用的就是该 `:s` 分支（`calculate_cms_momentum(..., :s, ...)`），
     所以即便传入不同 `t`，输出 `(k0_s, k_s)` 仍不变。

- 你要求的“随 `t` 出图（与 `Δs` 图同构）”已完成：
  1. 数据脚本：`scripts/analysis/relaxtime/t200_detm_im_vs_t_scan.jl`
  2. 绘图脚本：`scripts/analysis/relaxtime/plot_t200_detm_im_vs_t_scan.py`
  3. 输出图：`docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_detM_im_vs_t_scan.png`
  4. 汇总表：`D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_vs_t_scan_summary.csv`

- 该图与汇总表给出的数值结果：
  1. `k0_min = k0_max`、`k_abs_max = 0`（每个 `xi,Δs` 组合都如此）；
  2. `detM_im_range_abs = 0`；
  3. `invabs_detM_range_rel_to_mid = 0`、`abs_Dmixed_sq_range_rel_to_mid = 0`；
  4. 即“按 `t` 扫描的同构图”是严格平线，和理论推导及实现一致。

- 口径补充（防止过度外推）：
  1. 上述“无 `t` 依赖”是**当前 `:s` 道对象**结论，不是全局结论；
  2. 对 `:t`/`:u` 道，`calculate_cms_momentum` 中 `k0` 与 `k` 含 `t/u`，通常会有显著 `t` 依赖；
  3. 因而你这个问题本质上很关键：要先分清“看的是哪个散射道”。

## 24) 按你要求完成 `:t` 道对照图（与 `:s` 道并排）

- 你问的第一句先直接回答：
  1. **不是 `t=0`**。当前 `:s` 道分析里 `t` 取的是运动学区间内的值（可扫全区间），
     只是 `:s` 道定义使 `k0,k` 不随 `t` 变；
  2. 也就是说“`k0/k` 固定”来自散射道运动学关系，不是把 `t` 人为固定到 0。

- 你要的对照产物已生成：
  1. 数据脚本：`scripts/analysis/relaxtime/t200_detm_im_s_vs_t_channel_tscan.jl`
  2. 绘图脚本：`scripts/analysis/relaxtime/plot_t200_detm_im_s_vs_t_channel_tscan.py`
  3. 对照图：`docs/analysis/relaxtime/tauu_pos_uubaruubar_detM_im_s_vs_t_channel_tscan.png`
  4. 汇总表：`D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_s_vs_t_channel_tscan_summary.csv`

- 对照图读数（核心结论）：
  1. `:s` 道：
     - `k0` 常数、`k=0`，`detM_im`/`inv|detM|^2`/`|D_mixed|^2` 对 `t` 扫描均平坦；
  2. `:t` 道：
     - 在本过程下 `k0≈0`（等质量对称性），但 `k=\sqrt{-t}` 随 `t` 明显变化；
     - `detM_im` 仍基本不变（相对变化≈0）；
     - 但 `inv|detM|^2` 与 `|D_mixed|^2` 出现显著 `t` 依赖，且随 `Δs` 增大更强
       （如 `Δs=10` 时相对幅度约 `9.25~9.32` 与 `3.02~3.03`）。

- 这给出你问题的精确定稿：
  1. `t` 对“虚部链条”（本对象下的 `detM_im`）影响很弱；
  2. `t` 对“传播子模长/分母模长映射”（`inv|detM|^2`, `|D_mixed|^2`）影响显著；
  3. 所以“上游虚部机制图”和“传播子强度图”对 `t` 的灵敏度可以不同，这并不矛盾。

## 25) 继续上游：`detM_im` 的代数贡献切换（你要求的“更上游解释”）

- 你这次追问是对的：仅看 `Im(M00*M88)` 与 `Im(M08^2)` 还不够。
  本轮已再上推一层到恒等分解：
  - `Im(detM) = Re(M00)Im(M88) + Im(M00)Re(M88) - 2Re(M08)Im(M08)`；
  - 记作 `c1 + c2 + c3`。

- 新增产物（同路径、可复现）：
  1. 数据脚本：`scripts/analysis/relaxtime/t200_detm_im_upstream_term_switch_trace.jl`
  2. 符号翻转表脚本：`scripts/analysis/relaxtime/t200_detm_im_upstream_signflip_table.py`
  3. 绘图脚本：`scripts/analysis/relaxtime/plot_t200_detm_im_upstream_term_switch.py`
  4. 图：`docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_detM_im_term_switch.png`
  5. 表：
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_term_switch_trace.csv`
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_term_switch_summary.csv`
     - `D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_signflip_table.csv`

- 关键数值结论（`xi=0.34/0.36/0.38` 一致）：
  1. `detM_im` 的符号翻转都发生在 `Δs=5 -> 10`；
  2. `c1=Re(M00)Im(M88)` 的符号翻转更早，发生在 `Δs=1 -> 2`；
  3. `c3=-2Re(M08)Im(M08)` 的符号翻转发生在 `Δs=10 -> 15`；
  4. `Re(M08)` 也在 `Δs=10 -> 15` 翻符号，和 `c3` 翻转同位点；
  5. `c2=Im(M00)Re(M88)` 在全段保持负号，是稳定“下拉项”。

- 因而“先减后增”的机制可写成清晰三阶段：
  1. `decrease`（`Δs<=5`）：`c2<0` 与 `c3<0` 共振下拉，`detM_im` 下降到负谷；
  2. `recover`（`Δs≈10`）：`c1` 已转正且放大，与负的 `c2,c3` 竞争，`detM_im` 过零回正；
  3. `grow`（`Δs>=15`）：`Re(M08)` 翻号使 `c3` 也转正，`c1` 与 `c3` 同向上推，`detM_im` 继续增长。

- 对你问题的直接回答：
  - 是的，**有必要**往更上游推；而且这一层已经给出“项级别主导权切换”证据。
  - 下一层若再推进，最值得追的是：`Re(M08)` 为何在 `Δs≈10-15` 过零（即 `M08` 的耦合项与 `Π` 项竞争平衡点），它正是 `c3` 翻号与后段上冲的拐点控制量。

## 26) 再上推一层：`Re(M08)` 过零的根因（耦合常数项 vs 极化差项）

- 本轮按你“继续上推”的要求，直接定位 `Re(M08)` 的控制方程：
  - `M08 = K08_plus + c_offdiag*(Π_uu - Π_ss)`，
  - 因而 `Re(M08)=0` 的条件是
    `Re(Π_uu - Π_ss) = -K08_plus/c_offdiag`。

- 新增产物：
  1. 数据脚本：`scripts/analysis/relaxtime/t200_m08_re_zero_cross_trace.jl`
  2. 绘图脚本：`scripts/analysis/relaxtime/plot_t200_m08_re_zero_cross.py`
  3. 图：`docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross.png`
  4. 汇总：`D:\Desktop\Temp\relaxtime_t200_window\t200_m08_re_zero_cross_summary.csv`
  5. 细表：`D:\Desktop\Temp\relaxtime_t200_window\t200_m08_re_zero_cross_trace.csv`

- 关键结果（`xi=0.34/0.36/0.38`）：
  1. `Re(M08)` 零点线性插值都落在 `Δs≈12.48~12.51`；
  2. 对应目标阈值 `-K08_plus/c_offdiag` 分别为
     `0.16738 / 0.16602 / 0.16472`；
  3. 扫描中 `Re(Π_uu-Π_ss)` 随 `Δs` 单调下行并在上述区间穿过目标阈值；
  4. `Re(M08)` 由正转负与该阈值穿越一一对应；
  5. 公式残差 `Re(M08) - [K08_plus + c_offdiag*(ReΠ_uu-ReΠ_ss)]` 为机器精度量级（`~1e-18`），说明分解闭合。

- 与第 25 节对齐后的机制闭环：
  1. 第 25 节已证实 `c3=-2Re(M08)Im(M08)` 在 `Δs=10->15` 翻号；
  2. 本节给出其根因不是偶然数值波动，而是 `Re(M08)` 控制方程中的阈值穿越；
  3. 因而“`detM_im` 后段上冲”可以精确归因到
     `Re(Π_uu-Π_ss)` 穿越 `-K08_plus/c_offdiag` 后触发的 `c3` 翻转+增益。

- 你提到“横坐标改线性是否更好”这个建议成立：
  1. 已补充线性横轴版本图：
     `docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross_linear_x.png`；
  2. 同时保留原对数横轴图用于全尺度可见性；
  3. 在线性图中加入 `Δs∈[8,16]` 的局部线性拟合（每条 `xi` 单独拟合），
     用于直观看 `Re(M08)` 在过零邻域的近线性段与交点位置。

## 27) 图像“损坏”问题复盘与线性聚焦图修复说明

- 现象：
  1. `docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross.png`
     在部分查看器中无法打开，表现为“疑似损坏”；
  2. 线性图 `..._linear_x.png` 可正常打开。

- 诊断结论（文件本体未损坏）：
  1. PNG 签名、chunk CRC、IEND 结构均通过；
  2. 问题根因是导出画布高度异常（曾出现数十万像素级），超出常见查看器/库可处理上限。

- 根因（绘图导出配置）：
  1. 旧脚本全局启用了 `savefig.bbox="tight"`；
  2. log 分支中部分文字注释使用“数据坐标”，而第 3 子图 `residual` 量级接近机器精度；
  3. tight-bbox 在该组合下把画布包围盒拉到极端尺寸，导致导出图像高度异常。

- 已实施修复（仅改出图稳健性，不改物理数据）：
  1. 去除全局 `savefig.bbox="tight"`；
  2. log 分支注释统一改为轴坐标 `transAxes`；
  3. 新增线性聚焦参数并默认用于过零窗口：
     - `--focus-ds-min`（默认 `10.0`）
     - `--focus-ds-max`（默认 `15.0`）
     - `--focus-y-pad`（默认 `0.20`）
  4. 线性模式下改为 x/y 轴均线性，并自动聚焦 `Re(M08)` 过零邻域。

- 修复后产物：
  1. `docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross_linear_x.png`
     （x/y 线性 + 过零邻域聚焦）；
  2. `docs/analysis/relaxtime/tauu_pos_uubaruubar_uubar_to_uubar_reM08_zero_cross.png`
     （log 版同步修复后可正常打开）。

- 对图中 `M08 formula residual` 的定义与意义：
  1. 定义：
     `residual_formula = Re(M08)_total - [K08_plus + c_offdiag*(Re(Π_uu)-Re(Π_ss))]`；
  2. 这是“代数闭合误差”而非新的物理量，目的是检查
     `M08 = K08_plus + c_offdiag*(Π_uu-Π_ss)` 在实部层面的数值一致性；
  3. 若实现和数据链路一致，residual 应接近机器精度（本轮约 `1e-18` 量级），
     这说明我们用来解释过零点的控制方程是数值自洽的。

## 28) 术语与读图说明（Re(M08) 过零图专用）

- 为避免把“验证量”误读成“物理驱动量”，这里统一术语口径：
  1. **`Re(M08)`**：本节主解释对象，直接决定 `c3=-2Re(M08)Im(M08)` 的符号；
  2. **`Re(Piuu-Piss)`**：驱动变量，控制 `Re(M08)` 是否达到过零阈值；
  3. **`target_delta = -K08_plus/c_offdiag`**：阈值线（常数），与 `Re(Piuu-Piss)` 的交点对应 `Re(M08)=0`；
  4. **`M08 formula residual`**：代数闭合误差，仅用于数值自洽核验，不参与物理驱动判定。

- 三幅子图的推荐阅读顺序：
  1. 先看第 2 子图：`Re(Piuu-Piss)` 是否穿过 `target_delta`（判定是否会过零）；
  2. 再看第 1 子图：对应位置 `Re(M08)` 是否从正到负翻号（验证阈值映射）；
  3. 最后看第 3 子图：`residual` 是否维持机器精度量级（确认公式闭合与实现一致）。

- 误读防护（写作/汇报建议）：
  1. 不要把 residual 的波动当成“新机制”；
  2. residual 若显著偏离机器精度，应优先排查实现或数据链路，而非先做物理解读；
  3. 机制陈述应优先使用“阈值穿越 -> Re(M08) 翻号 -> c3 翻号 -> detM_im 后段上冲”这条主链。

- 与本节数据一致的口径锚点：
  1. 本轮 `residual_formula` 维持在 `~1e-18` 量级；
  2. 该量级支持“控制方程闭合”结论，可将解释重心放在阈值穿越与符号切换，而非数值噪声。
