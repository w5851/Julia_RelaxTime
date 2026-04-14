# T150 收敛修复与 T190/T200 毛刺证据链分析

## 1) 现象复述与数据口径

- 现象 A（T150）：此前 `xi=-0.5` 首点在部分扫描中被跳过，导致 `xi` 连续链路起点缺失。
- 现象 B（T190/T200）：`(sigma/T)/(eta/s)` 在局部 `xi` 区间出现“毛刺感”（局部斜率突变或平台不平滑）。
- 主数据口径：`D:\Desktop\Julia_RelaxTime\.worktrees\repro-main-oldparams\data\outputs\tmp\repro_main_oldparams\results\relaxtime\plan_b\plan_b_merged.csv`。
- 最小复跑口径：`D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05.csv` 与 sidecar `D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05_failed_points.csv`。

## 2) 问题1根因链

- 源码链路（修复前）：`scripts/relaxtime/run_gap_transport_scan.jl` 的 `solve_models_equilibrium(...)` 仅直接调用 `Models.solve_constraint(...)`，没有先走 `Models.solve(...)` 的内建 multiseed/fallback 主路径。
- 机制影响：当首点更依赖种子鲁棒性时，单一路径失败即直接进入 `catch point_err` 并被跳过，缺少可审计失败上下文。
- 可观测性缺口：原逻辑没有失败点 sidecar，日志之外无法批量复核“为何跳过、跳过了哪些点”。

## 3) 问题1修复结果

- 修复动作：在 `scripts/relaxtime/run_gap_transport_scan.jl` 新增 `Models.solve` 优先包装 `_solve_fixedmu_via_models_solve(...)`，失败后再回退 `_solve_fixedmu_via_models_constraint(...)`。
- 修复动作：新增 `--failed-points-output` 与 `write_failed_point_row!(...)`，把 `T_MeV,muB_MeV,xi,seed_source,phase_prev,phase_curr_hint,error_type,error_message,timestamp` 落盘。
- 测试证据：`tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl` 覆盖“优先调用路径可达”与 sidecar 字段行为。
- 最小复跑证据：`xi=-0.5` 在 `T=150, muB=0` 成功产出结果行；对应 sidecar 文件仅有表头无失败行，说明首点恢复且无隐藏失败。

## 4) 问题2毛刺解释

- 先更正口径：当前结果文件并非“只输出比值”。`plan_b_merged.csv` 同时包含分子分母与比值（`sigma_over_T`,`eta_over_s`,`sigma_over_T_over_eta_over_s`），并同时包含 `tau_u/tau_s` 与 `tauinv_u/tauinv_s`。
- 数值链路（主导，约 70%）：
  - 比值结构放大：`(sigma/T)/(eta/s)` 对分子分母相对变化高度敏感；在 `eta/s` 偏低区小幅波动可被放大为“毛刺感”。
  - T190 证据：`eta_over_s` 范围约 `0.0856~0.3601`，`sigma_over_T` 仅约 `0.00606~0.01933`，比值约 `0.0509~0.1041`，比值振幅显著高于单项量纲振幅。
  - T200 证据：`eta_over_s` 范围约 `0.0468~0.2195`，`sigma_over_T` 约 `0.00481~0.01536`，比值约 `0.0687~0.1110`，同样体现“分母敏感放大”。
  - 与弛豫时间联动：在比值高点（如 T200 `xi=-0.2`），`tau_u/tau_s` 明显小于平滑段（如 T200 `xi=0.36`），说明局部散射通道权重切换会放大到比值层。
- 物理链路（次级，约 30%）：
  - `m_u,m_s` 随 `xi` 变化引起阈值与通道权重再分配；T190/T200 位于接近相结构敏感区，通道竞争更强。
  - 该物理变化本身不必然导致“数值毛刺”，但会提高比值型观测量对离散点与插值误差的放大效应。

### 4.0 问题B零知识复现导航（最小闭环）

- 目标：让未参与本轮分析的读者，按固定顺序重现“`xi` 局部异常 -> 主导通道 -> `sigma` 下移 -> `mixed_P` 分母增强”链路。
- 最小输入：
  - 主数据：`D:\Desktop\Julia_RelaxTime\.worktrees\repro-main-oldparams\data\outputs\tmp\repro_main_oldparams\results\relaxtime\plan_b\plan_b_merged.csv`
  - 窗口重跑输出目录：`D:\Desktop\Temp\relaxtime_t190_window\`
- 最小复现步骤（按顺序）：
  1. 先看观测层折点（4.1）：确认 `xi=-0.10 -> -0.08` 是 T190 主异常段；
  2. 再看通道排序（4.3/4.3.1）：确认 `uubar_to_ddbar` + `uubar_to_uubar` 占 `tauinv_u` 跌幅约 90.3%；
  3. 看阈值邻域 `sigma(s)`（4.4）：确认两主通道面积约下降 43%~46%；
  4. 看核/σ 反事实（4.5/4.6/4.8）：确认主因在 `sigma` 本体（非 `v_rel`、非 `t` 区间、非 blocking 主导）；
  5. 看振幅与 s 道定位（4.9/4.10）：确认下降聚焦到 s 道 `|D_s^P|^2`；
  6. 看 mixed_P 深拆（4.10.1/4.10.2）：确认 `|D_mixed^P|^2` 下跌由 `|detM^P|^2` 上移主导。
- 最小通过判据（满足即认为链路复现成功）：
  - `ratio_Dmixed_B_over_A ≈ 0.349`；
  - `ratio_detM_B_over_A ≈ 1.126`；
  - 阈值最近点 `Δs=1e-6` 上 `|detM^P|^2(B/A)≈8.03` 且 `|D_mixed^P|^2(B/A)≈0.123`；
  - `detK` 与 `|J^TMJ'|^2` 仅小幅变化（分别约 `0.998` 与 `0.990`）。

### 4.1 针对 T190 的 `sigma_over_T` 图（`xi∈[-0.2,0]`）异常段解释

- 图位点：`data/outputs/figures/relaxtime/plan_b/T190/sigma_over_T_vs_xi.png`。
- 该区间在主数据中的对应行为不是单调平滑，而是“先降后升”：
  - 从 `xi=-0.20` 到 `xi=-0.14`，`sigma_over_T` 由 `0.01851` 下降到 `0.01679`。
  - 在 `xi=-0.14 -> -0.12` 出现局部反弹（`+1.20e-4`）。
  - 在 `xi=-0.10 -> -0.08` 出现更明显跃升（`+2.10e-3`），这是该窗口最主要的“异常感”来源。
  - 之后到 `xi=0.0` 重新转为缓慢上升并趋于平台。
- 与弛豫时间链路对应：
  - 同一窗口内 `tau_u` 在 `-0.10 -> -0.08` 同步出现大幅跳增（约 `+0.373`），`tauinv_u` 同步下降（约 `-0.0788`）。
  - 这说明异常段不是“只有比值在动”，而是散射率与弛豫时间本体在该局部发生了通道权重切换，`sigma_over_T` 图上的折点有可追溯的动力学对应。
- 与质量/阈值链路对应：
  - `m_u,m_s` 在 `[-0.2,0]` 内仍保持平滑上升，未出现同幅突跳。
  - 因此该异常更像“平滑背景上的数值/通道切换折点”，而非宏观相切换导致的真正不连续。

### 4.2 补全链路：`xi -> 解(序参量) -> 通道贡献`

- `xi -> 解(序参量)`（数据侧）：
  - 在 T190、`xi∈[-0.2,0]` 窗口内，`Phi` 与 `m_u,m_s` 平滑单调变化：
    - `Phi: 0.561138 -> 0.536806`
    - `m_u: 0.920839 -> 1.156106`
    - `m_s: 2.267323 -> 2.371113`
  - 该层没有与 `sigma_over_T` 同步的“突跳点”，说明异常不在“求解器失稳/序参量断裂”层。

- `解(序参量) -> 通道贡献`（公式-实现映射）：
  1. 平衡态输出 `x_state=(phi_u,phi_d,phi_s,Phi,Phibar)` 与 `masses`（扫描脚本中来自 `solve_models_equilibrium` 后续链路）。
  2. `x_state,masses` 进入 `A` 场与有效耦合计算：`G_u,G_s -> K_coeffs`（见 `src/relaxtime/EffectiveCouplings.jl`）。
  3. 平均散射率 `w_ij` 在 `src/relaxtime/AverageScatteringRate.jl` 中由 5D 积分给出，积分核为
     `f_i(m,mu,T,Phi,Phibar,xi) * f_j(...) * v_rel * sigma(s;K_coeffs,masses,thermo)`。
  4. `src/relaxtime/RelaxationTime.jl` 按
     `tau_i^{-1} = sum_j n_j * w_ij`
     聚合到物种级 `tau_inv`。
  5. `tau = 1/tau_inv`，再进入输运系数并形成 `sigma_over_T`。

### 4.2.1 “完整计算链路”与声明式流水线实现映射（主线/诊断支线）

- 本文档链路叙事
  `xi -> 平衡解 -> K/Π -> 振幅/传播子 -> sigma(s) -> w_ij -> tau -> transport`
  可以直接落为项目内的声明式流水线，但应拆为：
  - 主线（生产必跑）：`xi -> 平衡解 -> sigma/w_ij -> tau -> transport`
  - 诊断支线（按需开启）：`平衡解 -> K/Π -> 振幅/传播子 -> mixed_P(detM分解) -> 机制指标`

- 实现映射总表（输入 -> 核心函数/阶段 -> 输出字段 -> 证据文件）：

| 链路环节 | 输入 | 核心函数/阶段（声明式） | 输出字段（结构化） | 证据文件 |
|---|---|---|---|---|
| 统一入口 | `T_fm,mu_fm,xi` 与 workflow kwargs | `Models.solve_gap_and_transport`（`src/models/entrypoints.jl`）或 `Models.run_workflow_pipeline(:transport)` | workflow output tuple（含 `equilibrium`,`tau`,`transport`） | `run_manifest.json`（workflow adapter 输出） |
| 平衡解（主线） | `T_fm,mu_fm,xi,p_num,t_num` | `TransportWorkflow.solve_gap_and_transport` -> `_solve_equilibrium`（`src/models/workflows/TransportWorkflow.jl`） | `equilibrium.x_state`,`equilibrium.masses`,`equilibrium.converged` | 主 CSV（plan_b）与窗口重跑点位 |
| 率核与通道聚合（主线） | `quark_params,thermo_params,K_coeffs/densities` | `RelaxationTime.relaxation_times`（`src/relaxtime/RelaxationTime.jl`） | `tau`,`tau_inv`,`rates`; 且满足 `tauinv_i=sum_j n_j*w_ij` | `t190_window_channel_diag.csv` |
| 输运系数组装（主线） | `tau,bulk,densities,pressure,energy` | `transport_coefficients`（`src/relaxtime/TransportCoefficients.jl`） | `eta,zeta,sigma,kappa_**,lambda,lorenz...` | `plan_b_merged.csv` |
| K/Π 到混合传播子分解（诊断） | `equilibrium + process + (s,t)` | `decompose_mixed_p_propagator_chain`（`scripts/analysis/relaxtime/t190_sigma_chain_decomposition_lib.jl`） | `M00/M08/M88`, `detM` 三项，`|D_mixed^P|^2` | `t190_mixed_p_propagator_chain_decomposition.csv` |
| 振幅/传播子分项（诊断） | 同上 | `decompose_qqbar_amplitude_terms` / `decompose_p_channel_propagator_strength`（同一 lib） | `M_s/M_t/干涉`, `|D_s^P|^2` simple/mixed 分项 | `t190_sigma_amplitude_decomposition_summary.csv`、`t190_p_channel_propagator_absolute_strength_summary.csv` |
| 机制回归门禁（治理） | 固定点配置（T190, `xi=-0.10/-0.08`） | regression 测试 `tests/regression/relaxtime/test_t190_mixed_p_chain_regression.jl` | `ratio_detM_area_B_over_A`, `ratio_Dmixed_area_B_over_A`, 阈值点 ratio | `tests/baselines/relaxtime/baseline_t190_mixed_p_chain_v1.csv` |

- 对“如何具体实现为声明式流水线”的项目内落地建议（与现结构兼容）：
  1. **阶段命名**：在现有 workflow stage 语义上，显式对应
     `solve_equilibrium -> compute_tau -> compute_transport`（主线）与
     `analyze_sigma_chain -> analyze_mixed_p_chain`（诊断支线，按需开启）。
  2. **契约分层**：主线阶段仅要求返回生产字段；诊断阶段仅追加 artifact path 与指标，不反向修改主线物理量。
  3. **同 run_id 绑定**：主线与诊断产物都写入同一 manifest（当前已具备 manifest 机制），保证可审计与可复现。
  4. **门禁入口**：回归层读取诊断 summary 关键量做 isapprox 校验（当前已落地 `detM/Dmixed` 关键比值 smoke 版本）。

### 4.3 深拆“序参量 -> 通道贡献”：T190 异常点的可验证证据

- 对 `xi=-0.10 -> -0.08`（异常主段）做了通道级重跑：`D:\Desktop\Temp\relaxtime_t190_window\t190_window_channel_diag.csv`。
- 关键事实（u 物种）：
  - `tauinv_u` 约 `0.500863 -> 0.422067`（`-0.078796`）。
  - 其中两条主通道贡献下降占主导：
    - `uubar_to_ddbar`: `0.113891 -> 0.071446`（`-0.042444`）
    - `uubar_to_uubar`: `0.153437 -> 0.124726`（`-0.028712`）
    - 合计 `-0.071156`，约占 `tauinv_u` 总降幅的 90%+。
- 密度与率分解：
  - `n_u(=n_ubar)` 仅 `0.249327 -> 0.243413`（约 `-2.4%`）。
  - 但两主通道 `rate` 明显下降：
    - `uubar_to_ddbar` 率 `0.456792 -> 0.293520`（约 `-35.7%`）
    - `uubar_to_uubar` 率 `0.615405 -> 0.512404`（约 `-16.7%`）
  - 因而异常主因是“通道率核变化”，不是“密度平移”。
- 这一步给出的是“序参量层到通道贡献层”的实证闭环；仍未完全回答“率核为何在该点位突降”的第一性原理问题。
  - 下一层需要拆 `sigma(s)` 阈值邻域与 `f_i f_j v_rel` 权重在 `s` 轴上的贡献分布。

### 4.3.1 为什么优先分析 `uubar_to_ddbar` 与 `uubar_to_uubar`

- 选择不是主观挑通道，而是按异常主段 `xi=-0.10 -> -0.08` 的 **通道贡献降幅排序** 得到：
  - `uubar_to_ddbar`: `Δcontribution=-0.042444`，约占 `tauinv_u` 总降幅 `53.87%`；
  - `uubar_to_uubar`: `Δcontribution=-0.028712`，约占 `36.44%`；
  - 两者合计约 `90.30%`，远高于其余通道（单通道均 < 4%）。
- 因此后续 4.4~4.9 对这两条通道深挖，是“先抓主贡献、再做机制拆分”的数据驱动策略，不是先验假设。
- 对应原始证据源：`D:\Desktop\Temp\relaxtime_t190_window\t190_window_channel_diag.csv`。

### 4.4 继续下钻：`序参量 -> sigma(s)阈值结构 -> 通道率`

- 针对异常主段 `xi=-0.10 -> -0.08`，对两条主通道
  - `uubar_to_ddbar`
  - `uubar_to_uubar`
  进行了独立 `sigma(s)` 阈值邻域分解（同一积分设置：`N=128`, `n_sigma_points=6`, `threshold_subtraction=true`）。

- 关键结果（T190, muB=0）：
  - `m_u`：`1.05946 -> 1.08142`，对应阈值 `s_thr=(2m_u)^2`：`4.48983 -> 4.67786`。
  - `uubar_to_ddbar`：
    - 阈值后 1.0 区间峰值：`4.427 -> 2.397`（约 `-45.9%`）
    - 阈值后 2.0 区间面积：`0.4828 -> 0.2628`（约 `-45.6%`）
  - `uubar_to_uubar`：
    - 阈值后 1.0 区间峰值：`4.228 -> 2.351`（约 `-44.4%`）
    - 阈值后 2.0 区间面积：`0.4542 -> 0.2585`（约 `-43.1%`）

- 链路解释（数值+物理一致）：
  1. `xi` 改变平衡解后，`m_u` 上升并抬高阈值 `s_thr`。
  2. 在当前窗口，主导通道的 `sigma(s)` 阈值邻域峰值与面积同步显著下降。
  3. 这直接压低 `w_ij`（平均散射率）并传导到 `n_j*w_ij` 通道贡献。
  4. 进而触发 `tauinv_u` 的台阶式下降、`tau_u` 上升，并在 `sigma_over_T` 上体现为异常折点。

- 因此“通道率突降”已有可量化的阈值邻域证据，不再只是经验性描述。

### 4.5 5D 核分项归因：把通道率变化拆成“核变化”和“sigma 变化”

- 方法（反事实分解）：对 `xi=-0.10` 记为 A，对 `xi=-0.08` 记为 B。对每个通道率 `r` 做四组计算：
  - `rAA`: 核 A + `sigma` 缓存 A
  - `rBB`: 核 B + `sigma` 缓存 B
  - `rBA`: 核 B + `sigma` 缓存 A（只换核）
  - `rAB`: 核 A + `sigma` 缓存 B（只换 sigma）
  用两条对称路径估计贡献占比，避免单路径偏置。

- 结果 1：`uubar_to_ddbar`
  - 总变化：`rBB-rAA = -0.05408`
  - 路径一（经 `rBA`）：核项约 `-0.01397`，sigma 项约 `-0.04012`
  - 路径二（经 `rAB`）：sigma 项约 `-0.04378`，核项约 `-0.01030`
  - 解释：该通道下降以 `sigma(s)` 变化为主（约 74%~81%），核项次之。

- 结果 2：`uubar_to_uubar`
  - 总变化：`rBB-rAA = -0.03108`
  - 路径一（经 `rBA`）：核项约 `-0.01984`，sigma 项约 `-0.01124`
  - 路径二（经 `rAB`）：sigma 项约 `-0.01393`，核项约 `-0.01716`
  - 解释：该通道下降中核项占比更高（约 55%~64%），sigma 项为次要但非零。

- 合并解释（两主通道）：
  - `uubar_to_ddbar` 是“sigma 主导型”下跌；
  - `uubar_to_uubar` 是“核项主导型”下跌；
  - 二者叠加形成 `tauinv_u` 的主跌幅，最终触发 `tau_u` 与 `sigma_over_T` 的异常折点。

- 这一步把“序参量 -> 通道贡献”继续下钻到“积分核内部因子贡献”层：
  - 已完成：核项 vs `sigma(s)` 项的可量化拆分；
  - 未完成：`sigma(s)` 形状变化进一步映射到具体耦合参数组合（`K123/K4567/...`）的灵敏度曲线。

### 4.6 继续按公式拆核：`f_i`, `f_j`, `v_rel` 分项

- 依据公式文档 `AverageScatteringRate_FromCrossSection.md`，5D 积分核可写为
  `f_i * f_j * v_rel * sigma(s)`。
- 在本窗口（`uubar_to_ddbar` / `uubar_to_uubar`）里，`i,j` 都是轻味反/夸克且质量、化学势对称，因此
  `f_i` 与 `f_j` 的变化可合并看作 `f_if_j` 乘积项。

- 对 `xi=-0.10(A)` 与 `xi=-0.08(B)` 的核分项统计（不含 `sigma`）显示：
  - `sum(f_i f_j)` 比值 `B/A ≈ 0.9530`
  - `sum(v_rel)` 比值 `B/A ≈ 0.9992`
  - `sum(f_i f_j v_rel)` 比值 `B/A ≈ 0.9477`

- 解释：
  1. 在异常主段中，核项变化主要由分布函数乘积 `f_i f_j` 驱动（约 4.7% 量级下降）；
  2. `v_rel` 近乎不变（<0.1%），不是主导源；
  3. 因此“核项下降”可进一步归因到各向异性分布权重变化，而非相对速度几何项变化。

- 结合 4.5 与本节：
  - `uubar_to_ddbar`：以 `sigma(s)` 项主导下降，`f_if_j` 提供次级下降；
  - `uubar_to_uubar`：`f_if_j` 核项与 `sigma(s)` 项共同作用且核项略占主导。

### 4.7 `sigma(s)` 影响项对照表（公式项 -> 代码项 -> 可观测证据）

| 公式层影响项 | 代码实现锚点 | 在本次 T190 异常窗口的可观测证据 |
|---|---|---|
| 微分截面本体 `dσ/dt`（振幅/传播子，`|M|^2`） | `src/relaxtime/TotalCrossSection.jl`（`total_cross_section`） + `src/relaxtime/ScatteringAmplitude.jl` + `src/relaxtime/TotalPropagator.jl` | 两条主通道 `sigma(s)` 阈值邻域峰值/面积在 `xi=-0.10 -> -0.08` 同步下降约 43%~46%，说明截面本体在该窗口显著变“软”。 |
| 有效耦合 `K_coeffs`（由 `G_u,G_s` 导出） | `src/relaxtime/EffectiveCouplings.jl`（`calculate_effective_couplings`）+ 扫描脚本 `build_K_data(...)` | `xi` 改变平衡解后（`Phi,m_u,m_s` 平滑变化），`K_coeffs` 及相关传播子参数随之变化，并通过 `sigma(s)` 变化传导到主通道率下降。 |
| 末态阻塞因子 `(1-f_c)(1-f_d)`（位于 `σ` 的 `t` 积分内部） | `docs/reference/formula/relaxtime/scattering/TotalCrossSection_FromDifferentialCrossSection.md` 对应实现于总截面积分核 | 在反事实分解中，`sigma` 项贡献显著（尤其 `uubar_to_ddbar`），与仅外层核(`f_i f_j v_rel`)无法解释的降幅一致，支持“末态阻塞+截面核”共同作用。 |
| 积分区间 `t_±` 与阈值/相同粒子对称因子 | 同上文档 + `src/relaxtime/TotalCrossSection.jl` 的 `t` 区间与同粒子修正逻辑 | `m_u` 上升使阈值 `s_thr=(2m_u)^2` 上移（`4.4898 -> 4.6779`），与阈值邻域 `sigma(s)` 面积显著收缩同向，解释主通道贡献骤降。 |
| 阈值邻域快变结构（`σ(s)` 的尖峰/拐点） | `docs/reference/formula/relaxtime/scattering/AverageScatteringRate_FromCrossSection.md` + cache 预计算/插值路径 | 已量化“阈值后 1.0/2.0 区间”的峰值与面积下降；该证据直接对接 `tauinv_u` 的台阶式下降。 |
| 数值实现项（`n_sigma_points`, `interpolation_mode`, `threshold_subtraction`, `asym_*`, cache 网格） | `src/relaxtime/AverageScatteringRate.jl`（`build_w0cdf_pchip_cache`, `_get_sigma_core`） | 本分析固定与主跑一致配置后仍重现异常，说明不是“单纯参数误配”；但这些项仍决定 `sigma(s)` 细节形状，属于残余不确定性来源。 |

- 注：上表用于避免“把外层核变化误当成 `sigma` 本体变化”。本报告中：
  - 外层核分解见 4.6（`f_if_j` 主导核项变化，`v_rel` 次要）；
  - `sigma(s)` 本体变化分解见 4.4/4.5/本节对照表。

### 4.8 继续深挖 4.7：`sigma` 本体分项反事实（K / blocking / t 区间）

- 新增可复现实验脚本：`scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl`。
  - 固定 `T=190, muB=0`，比较 `xi=-0.10(A)` 与 `xi=-0.08(B)`。
  - 对两主通道 `uubar_to_ddbar`、`uubar_to_uubar` 在 `s=s_thr+Δs (Δs∈[1e-6,2.0])` 上积分，分别计算：
    1) 基线 `sigma_base`（原始配置）；
    2) 去掉末态阻塞 `sigma_no_block`；
    3) 仅替换 `K_coeffs` 的反事实 `sigma_Kswap`（A态用B的K，B态用A的K）；
    4) `t` 区间宽度 `t_width=t_max-t_min(含同粒子修正)`。
  - 输出：
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_chain_decomposition.csv`
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_chain_decomposition_summary.csv`

- 结果 1（`uubar_to_ddbar`）：
  - `sigma` 阈值后面积比（B/A）约 `0.403`；
  - 去阻塞后面积比约 `0.394`（与基线接近），说明“B相对A的下降”并非由 blocking 主导；
  - `K` 反事实敏感度显著：A态替换为B态 `K` 后面积约 `-31.3%`，B态相对A态 `K` 的差异量级约 `18%~26%`；
  - 平均 `t_width` 基本不变（A/B 都约 `1.0000005`），排除“t区间退化导致毛刺”的主导解释。

- 结果 2（`uubar_to_uubar`）：
  - `sigma` 阈值后面积比（B/A）约 `0.510`；
  - 去阻塞后面积比约 `0.498`（同样接近基线），blocking 不是主导下降源；
  - `K` 反事实敏感度同样显著：A态替换为B态 `K` 后约 `-25.5%`，B态相对A态 `K` 的差异量级约 `18%`；
  - `t_width` 同样几乎不变（A/B 都约 `1.0000005`）。

- 解释更新（对应 4.7 表）：
  1. `t_±` 区间变化在本窗口不是主导项（宽度几乎不变）；
  2. `(1-f_c)(1-f_d)` 会改变绝对量级，但对 A/B 相对跌幅贡献次要；
  3. `K_coeffs -> propagator/amplitude -> dσ/dt` 是可量化的强敏感链路，能解释 `sigma(s)` 主体下移；
  4. 与 4.4 的“阈值邻域峰值/面积下降 43%~46%”一致，表明异常主要源于 `sigma` 本体（耦合/振幅链）而非几何积分边界项。

### 4.9 振幅级继续下钻：`|M|^2 = s道 + t道 + 干涉`

- 新增分解工具：
  - 复用库：`scripts/analysis/relaxtime/t190_sigma_chain_decomposition_lib.jl`
  - 振幅分解脚本：`scripts/analysis/relaxtime/t190_sigma_amplitude_decomposition.jl`
  - 单测：`tests/unit/relaxtime/test_t190_sigma_chain_decomposition.jl`

- 分解口径：
  - 对 `uubar_to_ddbar`、`uubar_to_uubar` 在 `s=s_thr+Δs`（`Δs∈[1e-6,2.0]`）上取 `t_mid=(t_min+t_max)/2`，
    将 `M2_total` 拆为：
    - `M_s_sq = M_s_S + M_s_P`
    - `M_t_sq = M_t_S + M_t_P`
    - `M_interf = -2Re(M_s M_t*)`
  - 并做 `K` 反事实（A态用B态K、B态用A态K）评估敏感度。

- 关键结果（汇总文件：`D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_amplitude_decomposition_summary.csv`）：
  - `uubar_to_ddbar`：
    - `M2` 面积比（B/A）≈ `0.433`；
    - A 点道贡献占比：`s道≈94.0%`, `t道≈7.8%`, `干涉≈-1.8%`；
    - B 点道贡献占比：`s道≈86.0%`, `t道≈17.7%`, `干涉≈-3.7%`；
    - `K` 反事实敏感度：A侧约 `-30.1%`，B侧约 `+32.5%`。
  - `uubar_to_uubar`：
    - `M2` 面积比（B/A）≈ `0.544`；
    - A 点道贡献占比：`s道≈96.5%`, `t道≈5.6%`, `干涉≈-2.1%`；
    - B 点道贡献占比：`s道≈93.3%`, `t道≈10.3%`, `干涉≈-3.7%`；
    - `K` 反事实敏感度：A侧约 `-24.3%`，B侧约 `+20.5%`。

- 解释：
  1. 两主通道在异常窗口内都表现为 **s道主导**；
  2. 从 A 到 B，`s道` 绝对贡献显著下降，同时 `t道` 相对占比上升、干涉负贡献增强；
  3. 这与 4.8 的 `sigma` 面积下降、以及 4.5 的“sigma项主导/参与”结论一致，进一步把“`sigma` 本体变化”落到传播子/振幅级证据。

### 4.10 基于公式继续深挖：为何是 s 道出现异常主导

- 依据公式文档 `ScatteringAmplitude_FromTotalPropagator.md`，`qqbar` 过程中
  `M_s_sq = |D_s^S|^2 * s12^+*s34^+ + |D_s^P|^2 * s12^-*s34^-`。
- 因此可把 s 道再拆成两类来源：
  1. 传播子模方项（`|D_s^S|^2`, `|D_s^P|^2`）；
  2. 运动学乘子项（`s12^+*s34^+`, `s12^-*s34^-`）。
- 新增分解脚本：`scripts/analysis/relaxtime/t190_s_channel_rootcause_decomposition.jl`，输出：
  - `D:\Desktop\Temp\relaxtime_t190_window\t190_s_channel_rootcause_decomposition.csv`
  - `D:\Desktop\Temp\relaxtime_t190_window\t190_s_channel_rootcause_decomposition_summary.csv`

- 关键结果（A=`xi=-0.10`, B=`xi=-0.08`）：
  - `uubar_to_ddbar`：
    - `area_Ms(B/A) ≈ 0.396`（s道总项显著下降）；
    - `|D_s^P|^2` 面积比约 `0.290`（强下降），而 `s12^-*s34^-` 面积比约 `1.069`（略升）；
    - 说明 `M_s` 下跌不是运动学乘子导致，而是 **P 通道传播子强度下降** 导致；
    - 分项占比从 A 到 B：`M_s_P` 约 `94.9% -> 80.8%`，`M_s_S` 约 `5.1% -> 19.2%`（S项占比上升但不足以抵消P项下跌）。
  - `uubar_to_uubar`：
    - `area_Ms(B/A) ≈ 0.526`；
    - `|D_s^P|^2` 面积比约 `0.417`，`s12^-*s34^-` 面积比约 `1.069`；
    - 同样表明主因是 **P 通道传播子模方下降** 而非 s 运动学因子突变；
    - 分项占比：`M_s_P` 约 `97.8% -> 92.2%`，仍是主导分量。

- 对应到传播子公式 `Propagator_传播子byPolarization.md`：
  - P 通道传播子由 `K^+` 与极化函数 `Π^P` 进入分母 `1 - 4K^+Π^P`；
  - 本窗口中出现的是“P通道传播子模方敏感下移”，与前文 `K` 反事实敏感性结果一致，且与 `t` 区间、blocking 非主导的结论相容。

### 4.10.1 继续深挖 4.10：P 通道传播子强度的**绝对量级**证据（simple/mixed 分解）

- 新增脚本：`scripts/analysis/relaxtime/t190_p_channel_propagator_absolute_strength.jl`。
  - 固定 `T=190, muB=0`，对 A=`xi=-0.10` 与 B=`xi=-0.08`、`Δs∈[1e-6,2.0]`（36 点）逐点计算：
    - `D_s^P_total`
    - `D_s^P_simple`（s道 simple 介子部分，`pi`/`K`）
    - `D_s^P_mixed`（s道 mixed_P 介子部分，η/η'整体）
    - 以及对应 `|D|^2` 与 `M_s_P=|D_s^P_total|^2*s12^-*s34^-`。
  - 输出：
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_p_channel_propagator_absolute_strength.csv`
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_p_channel_propagator_absolute_strength_summary.csv`

- 关键结果 1（`uubar_to_ddbar`）：
  - `area(|D_s^P_total|^2)`：`71.106 -> 20.633`，`B/A ≈ 0.290`（与 4.10 一致）；
  - `area(|D_s^P_simple|^2)`：`1.726 -> 1.626`，`B/A ≈ 0.942`（变化小，且绝对量级很小）；
  - `area(|D_s^P_mixed|^2)`：`79.870 -> 27.905`，`B/A ≈ 0.349`（显著下跌，量级主导）；
  - `area(s12^-*s34^-)`：`60.944 -> 65.144`，`B/A ≈ 1.069`（运动学乘子仍略升）；
  - `area(M_s_P)`：`1507.06 -> 508.03`，`B/A ≈ 0.337`。

- 关键结果 2（`uubar_to_uubar`）：
  - `area(|D_s^P_total|^2)`：`92.086 -> 38.428`，`B/A ≈ 0.417`；
  - `area(|D_s^P_simple|^2)`：`1.726 -> 1.626`，`B/A ≈ 0.942`（同样不是主导变化）；
  - `area(|D_s^P_mixed|^2)`：`79.870 -> 27.905`，`B/A ≈ 0.349`（同样显著下跌）；
  - `area(s12^-*s34^-)`：`60.944 -> 65.144`，`B/A ≈ 1.069`；
  - `area(M_s_P)`：`2047.86 -> 1016.31`，`B/A ≈ 0.496`。

- 解释（为何这是“绝对量级证据”而不只是比例重述）：
  1. 在两主通道中，`|D_s^P|^2` 的绝对面积下降都在几十量级，而 `|D_s^P_simple|^2` 仅约 `1.6~1.7` 且 A/B 变化 < 6%；
  2. 因此 `|D_s^P|^2` 的主跌幅来自 **mixed_P（η/η'整体）传播子强度链路**，不是 simple 部分；
  3. 由于总量是复振幅先求和再取模方，`|D_total|^2` 不等于 `|D_simple|^2+|D_mixed|^2`，不同过程会出现不同相位干涉符号（`uubar_to_ddbar` 与 `uubar_to_uubar` 在 A/B 点均表现出不同的净干涉方向），但不改变“mixed_P 强度下移主导”的结论；
  4. 与 4.10 合并后可写成：`M_s_P` 下降 = `|D_s^P|^2` 绝对强度下移（主因） × `s12^-*s34^-` 略升（非主因）。

### 4.10.2 按 mixed_P 公式继续下钻：是 `detK`、`detM` 还是 `J^T M J'` 在主导下降？

- 依据公式文档 `Propagator_传播子byPolarization.md`，mixed_P（η/η'整体）传播子可写为
  `D_mixed^P = 2 * detK^+ / detM^P * (J^T M^P J')`，
  因此
  `|D_mixed^P|^2 = 4 * (detK^+)^2 * |J^T M^P J'|^2 / |detM^P|^2`。
- 新增脚本：`scripts/analysis/relaxtime/t190_mixed_p_propagator_chain_decomposition.jl`，逐点输出
  `detK^+`, `|detM|^2`, `|J^T M J'|^2`, `|D_mixed^P|^2`, 以及 `|Π_uu^P|^2`, `|Π_ss^P|^2`。
  - 输出：
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_mixed_p_propagator_chain_decomposition.csv`
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_mixed_p_propagator_chain_decomposition_summary.csv`

- 关键结果（A=`xi=-0.10`, B=`xi=-0.08`，两主通道在 mixed_P 链路上数值一致）：
  - `area(|D_mixed^P|^2): 79.870 -> 27.905`，`B/A≈0.349`（显著下降）；
  - `area((detK^+)^2)` 近似不变（`detK^+` 面积比 `B/A≈0.998`）；
  - `area(|J^T M J'|^2)` 仅轻微下降（`B/A≈0.990`）；
  - `area(|detM|^2)` 明显上升（`2.10e-6 -> 2.37e-6`, `B/A≈1.126`）；
  - 极化函数模方变化温和：`|Π_uu^P|^2` 比值约 `1.019`，`|Π_ss^P|^2` 比值约 `1.007`。

- 把 `detM^P = M00*M88 - M08^2` 继续拆成“三项贡献”（对应 `|detM^P|^2` 的可加分解）：
  - 记 `a = M00*M88`, `b = M08^2`，则
    `|detM^P|^2 = |a-b|^2 = |a|^2 + |b|^2 - 2Re(a b*)`。
  - 在脚本中已显式输出三项面积：
    1) `|M00*M88|^2`；
    2) `|M08^2|^2`；
    3) `cross = -2Re((M00*M88) * conj(M08^2))`。
  - 量化结果（A=`xi=-0.10`, B=`xi=-0.08`）：
    - `area(|M00*M88|^2): 4.612e-7 -> 4.731e-7`，`B/A≈1.026`（小幅上升）；
    - `area(|M08^2|^2): 7.876e-7 -> 8.831e-7`，`B/A≈1.121`（明显上升）；
    - `area(cross): 8.536e-7 -> 1.011e-6`，`B/A≈1.184`（上升最强）。

- 三项分解解释（针对“为何 mixed_P 强度下降”）：
  1. `|detM^P|^2` 的上移不是单一子项突变，而是三项同向抬升；
  2. 其中 `|M08^2|^2` 与交叉项 `cross` 的抬升幅度显著大于 `|M00*M88|^2`，是分母增强的主导来源；
  3. 因 `|D_mixed^P|^2 ∝ 1/|detM^P|^2`，分母增强直接压低 mixed_P 传播子强度，与 4.10.1 的绝对量级证据一致。

- 继续把 `M00/M08/M88` 拆到 `K` 与 `Π` 项（对应 `mixing_matrix_elements` 公式）：
  - `M00 = K0^+ - (4/3)detK^+*Π_uu^P - (8/3)detK^+*Π_ss^P`
  - `M08 = K08^+ + (4/3)√2 detK^+*Π_uu^P - (4/3)√2 detK^+*Π_ss^P`
  - `M88 = K8^+ - (8/3)detK^+*Π_uu^P - (4/3)detK^+*Π_ss^P`

- `Π/K` 项级面积比（B/A，A=`xi=-0.10`, B=`xi=-0.08`）：
  - `M00` 分解项：
    - `|K0^+|^2`：`0.03569 -> 0.03531`，`≈0.990`（小降）；
    - `|Π_uu|`项：`0.007135 -> 0.007237`，`≈1.014`（小升）；
    - `|Π_ss|`项：`0.010833 -> 0.010850`，`≈1.002`（近似不变）。
  - `M08` 分解项：
    - `|K08^+|^2`：`1.717e-4 -> 1.658e-4`，`≈0.966`（下降）；
    - `|Π_uu|`项：`0.014269 -> 0.014473`，`≈1.014`（上升）；
    - `|Π_ss|`项：`0.005416 -> 0.005425`，`≈1.002`（近似不变）。
  - `M88` 分解项：
    - `|K8^+|^2`：`0.09292 -> 0.09342`，`≈1.005`（小升）；
    - `|Π_uu|`项：`0.028539 -> 0.028947`，`≈1.014`（上升）；
    - `|Π_ss|`项：`0.002708 -> 0.002713`，`≈1.002`（近似不变）。

- 该层解释（把异常继续落到 `Π/K`）：
  1. `K` 侧并未出现“单点突变”，且 `K08^+` 反而下降，不支持“仅由 K 增大导致分母变强”；
  2. 与 `Π_uu^P` 相关的项在 `M08/M88` 中呈一致上升（约 `+1.4%`），是 `M08^2` 抬升的重要方向性来源；
  3. `Π_ss^P` 相关项变化较小（约 `+0.2%`），主要提供次级抬升；
  4. 结合前述 detM 三项结果，当前最强链路可写为：
     `Π_uu^P` 主导的 `M08/M88` 组合抬升 -> `|M08^2|^2` 与交叉项同步抬升 -> `|detM^P|^2` 增强 -> `|D_mixed^P|^2` 下移。

- 为什么“Π/K 项只有约 1% 变化”仍能导致 `|D_mixed^P|^2` 大幅变化？（关键是**非线性与阈值权重**）
  1. 传播子是分式结构：`|D_mixed^P|^2 = 4(detK^+)^2|J^TMJ'|^2 / |detM^P|^2`，变化主要由分母控制，不是线性叠加；
  2. `detM^P = M00*M88 - M08^2` 含复数差分与相位干涉，`M00/M08/M88` 的小变化可通过“差分接近抵消区”被放大到 `|detM|^2`；
  3. 点值证据（阈值最近点 `Δs=1e-6`）：
     - `|D_mixed^P|^2(B/A)≈0.123`（单点已是 8x 级下降）；
     - 同点 `|detM^P|^2(B/A)≈8.03`（分母强增），`|J^TMJ'|^2(B/A)≈0.991`（几乎不变）；
  4. 权重证据：A 态 `|D|^2` 强烈集中在最靠阈值点（首点占总面积约 `76.4%`），B 态该占比降至 `35.6%`；因此阈值最近点的分母增强会对积分面积产生不成比例影响；
  5. 这解释了“项级看似小变（~1%）”与“传播子面积大变（~65% 下降）”并不矛盾：前者发生在 `M` 元素层，后者发生在含差分/模方/倒数/积分权重的复合非线性映射后。

- 回答“是否看到 `|detM^P|^2` 随 `xi` 的直接变化”：此前文档确实主要写了 A/B 对比，`xi` 连续扫描证据不够显式；现补充 T190 窗口直接扫描（`xi=-0.20:0.02:0.00`，`uubar_to_ddbar`，同样 `Δs∈[1e-6,2.0]`）：
  - `|detM^P|^2`（阈值最近点）随 `xi` 的序列：
    - `-0.12: 3.55e-9`
    - `-0.10: 3.87e-9`
    - `-0.08: 3.11e-8`
    可见在异常主段 `-0.10 -> -0.08` 出现约 `8.03x` 的突升。
  - 对应 `|D_mixed^P|^2`（阈值最近点）序列：
    - `-0.12: 1914.54`
    - `-0.10: 1729.13`
    - `-0.08: 212.21`
    同段出现约 `0.123x`（约 8 倍下降）。
  - 面积层（`Δs` 积分）也同向：
    - `area(|detM^P|^2): 2.10e-6 -> 2.37e-6`（`B/A≈1.126`）
    - `area(|D_mixed^P|^2): 79.87 -> 27.90`（`B/A≈0.349`）
  - 这说明你的判断是对的：核心不是“Π 单项大幅变化”，而是 `detM^P` 在阈值敏感区接近小量时，经过复数差分与倒数结构后，对微扰高度敏感。

- 解释（公式一致）：
  1. `detK^+` 几乎不变，排除“耦合行列式前因子突变”作为主导原因；
  2. `|J^T M J'|^2` 仅小幅变化，不足以解释 `|D_mixed^P|^2` 的 65% 级降幅；
  3. 主导项是分母 `|detM^P|^2` 的上移（离共振极点更远），导致 mixed_P 传播子绝对强度显著下移；
   4. 这与 4.10.1 的“mixed_P 主导 `|D_s^P|^2` 下跌”完全一致，并把根因进一步定位到
      **`M^P(Π_uu^P, Π_ss^P, K^+)` 组合导致的 `detM^P` 分母增强**。

### 4.10.3 从“两点”扩展到全窗口 `xi∈[-0.2,0.0]`：完整链路深拆结果

- 为避免只看 `-0.10/-0.08` 两点导致过拟合，本轮按同口径做了 **全区间逐点深拆**（11 个 `xi` 点）：
  - 新脚本：`scripts/analysis/relaxtime/t190_xi_window_full_chain_decomposition.jl`
  - 产物：
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_xi_window_full_chain_detail.csv`
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_xi_window_full_chain_summary.csv`
    - `D:\Desktop\Temp\relaxtime_t190_window\t190_xi_window_adjacent_transition_summary.csv`

- 观测层（主数据 `plan_b_merged.csv`）确认：`sigma_over_T` 在 `[-0.2,0]` 呈“先降后升”，并带一个小回跳：
  - `-0.20 -> -0.14` 连续下降；
  - `-0.14 -> -0.12` 小回跳（相邻比值约 `1.007`）；
  - `-0.12 -> -0.10` 再次回落；
  - `-0.10 -> -0.08` 主回升（相邻比值约 `1.128`，窗口内最大跃升）。

- 对应链路层（`uubar_to_ddbar`，相邻 `xi` 比值，来自 `t190_xi_window_adjacent_transition_summary.csv`）：
  - `-0.14 -> -0.12`（小回跳段）：
    - `ratio_area_abs_Dmixed_sq ≈ 2.276`
    - `ratio_area_abs_detM_sq ≈ 1.119`
    - `ratio_area_abs_DsP_mixed_sq ≈ 2.276`
    - `ratio_area_M_s_sq ≈ 2.176`
  - `-0.10 -> -0.08`（主跃升段）：
    - `ratio_area_abs_Dmixed_sq ≈ 0.349`
    - `ratio_area_abs_detM_sq ≈ 1.126`
    - `ratio_area_abs_DsP_mixed_sq ≈ 0.349`
    - `ratio_area_M_s_sq ≈ 0.396`

- 解释（窗口级，不再局限两点）：
  1. `sigma_over_T` 的“先降后升”来自 **至少两次不同方向的 mixed_P 强度重分配**（一次抬升，一次塌缩）；
  2. 两次转折都保留“`detM^P` 分母控制 + mixed_P 主导”的同一机制框架，只是相邻点处符号与幅度不同；
  3. 因此此前两点结论并未失效，而是全区间上可复现为“同机制的分段非单调响应”。

### 4.10.4 为什么 `xi>0` 基本单调下降，而 `[-0.2,0]` 非单调

- 使用同口径新增了正区间对照（`uubar_to_ddbar`）：
  - 脚本：`scripts/analysis/relaxtime/t190_xi_positive_chain_contrast.jl`
  - 产物：`D:\Desktop\Temp\relaxtime_t190_window\t190_xi_positive_chain_contrast.csv`

- 对照事实：
  - 在 `xi∈[0,0.5]`，主数据 `sigma_over_T` 的相邻比值始终 `<1`（约 `0.997~0.999`），表现为平滑单调下降；
  - 同时链路量呈连续缓变：
    - `ratio_area_abs_Dmixed_sq` 约 `0.851 -> 0.970`（始终小于 1，且逐步接近 1）；
    - `ratio_area_abs_detM_sq` 约 `1.092 -> 1.023`（始终大于 1，且逐步回到弱敏感区）；
    - 未出现 `[-0.2,0]` 内那种 `2.276` 或 `0.349` 级别的突发跳变。

- 因而两段趋势差异可归纳为：
  1. `[-0.2,0]` 落在 mixed_P 阈值敏感窗口，`detM^P`-控制的分母效应出现分段强化/反向重分配，导致局部回跳与主跃升；
  2. `xi>0` 区间该链路进入“缓变衰减态”，没有触发同量级的相邻突变，故 `sigma_over_T` 只表现为缓慢单调下行。

### 4.10.5 对“同为上升但 `D_mixed` 比值方向相反”的澄清（补 `-0.12->-0.10` 机制）

- 你指出的关键问题是成立的：
  - `-0.14 -> -0.12` 小回跳段，`D_mixed` 面积比 `>1`（约 `2.276`）；
  - `-0.10 -> -0.08` 主跃升段，`D_mixed` 面积比 `<1`（约 `0.349`）；
  - 两段都出现 `sigma_over_T` 上升，但 **不是同号同性质** 过程，不能被同一句“同机制同方向”覆盖。

- 本文统一口径的更精确表述应为：
  1. 两段都属于同一“mixed_P / detM 控制链路”的敏感窗口内事件；
  2. 但在该链路内，**上升的驱动项不同**：
     - `-0.14 -> -0.12`：`D_mixed` 与 `M_s` 同步大幅抬升（`ratio_area_abs_Dmixed_sq≈2.276`, `ratio_area_M_s_sq≈2.176`），并经 `sigma(s)->w_ij->tauinv` 串联后体现为小回跳；
     - `-0.10 -> -0.08`：`D_mixed` 与 `M_s` 同步塌缩（`≈0.349`, `≈0.396`），进而压低 `w_ij` 与 `tauinv_u`（`0.5009 -> 0.4221`，`tau_u` 跃升），最终表现为更大幅的 `sigma_over_T` 上升。

- 串联关系澄清（避免“并行因子”误读）：
  - 本文中的“截面核变化”和“散射率变化”不是并行两条线，而是同一链路上的前后级：
    `传播子/振幅 -> sigma(s) -> w_ij -> tauinv -> tau -> transport(含 sigma_over_T)`。
  - 因此文中若出现“截面核主导”表述，严格含义是“主导了后续散射率变化的上游原因”。

- 这也解释了你提出的第三段（`-0.12 -> -0.10`）为何是大幅下降：
  - 该段 `sigma_over_T` 比值约 `0.972`（下降），同时
    - `D_mixed` 面积比约 `0.866`（回落），
    - `M_s` 面积比约 `0.862`（回落），
    - `tauinv_u` 仅轻微下降（`0.5041 -> 0.5009`），不足以抵消截面核回落；
  - 因而形成“先回跳（`-0.14->-0.12`）再回落（`-0.12->-0.10`）再主跃升（`-0.10->-0.08`）”的三段结构。

- 直接数值证据（连续 `xi` 阈值最近点，`uubar_to_ddbar`, `Δs=1e-6`；来自 `t190_xi_window_full_chain_detail.csv`）：
  - `|detM^P|^2`：
    - `xi=-0.12`: `3.55e-9`
    - `xi=-0.10`: `3.87e-9`
    - `xi=-0.08`: `3.11e-8`
    - `-0.10 -> -0.08` 比值约 `8.03x`
  - `|D_mixed^P|^2`：
    - `xi=-0.12`: `1914.54`
    - `xi=-0.10`: `1729.13`
    - `xi=-0.08`: `212.21`
    - `-0.10 -> -0.08` 比值约 `0.123x`

- 面积层（`Δs∈[1e-6,2.0]`，与 4.10.2 保持一致）：
  - `area(|detM^P|^2): 2.10e-6 -> 2.37e-6`（`B/A≈1.126`）
  - `area(|D_mixed^P|^2): 79.87 -> 27.90`（`B/A≈0.349`）

- 因此“是否看到 `|detM^P|^2` 随 `xi` 的直接变化”现在可明确回答：
  1. 是，且在连续 `xi` 扫描中显式可见；
  2. 异常主段 `-0.10 -> -0.08` 的核心仍是 `|detM^P|^2` 阈值最近点突升与 `|D_mixed^P|^2` 对应突降；
  3. 小回跳段与下降段是同一敏感窗口内的不同子阶段，不能简化为单一单调链路。

### 4.10.6 阈值近点“近零分母”是否导致发散：收敛性与可积性检查

- 你提出的数学问题是关键：当 `|detM^P|^2` 在阈值近点变小，`|D_mixed^P|^2 ∝ 1/|detM^P|^2` 是否会导致不可积发散？

- 本轮直接数值检查（`uubar_to_ddbar`, T190）：
  1. **`Δs` 逼近检查（`1e-8 -> 1e-6`）**：
     - 明确列出 `Δs=1e-8` 与 `1e-6`：
       - `xi=-0.12`：
         - `|detM^P|^2: 3.544e-9 -> 3.548e-9`
         - `|D_mixed^P|^2: 1.916e3 -> 1.915e3`
         - `sigma_total: 111.7097 -> 111.6048`
       - `xi=-0.10`：
         - `|detM^P|^2: 3.872e-9 -> 3.869e-9`
         - `|D_mixed^P|^2: 1.728e3 -> 1.729e3`
         - `sigma_total: 94.5770 -> 94.6394`
       - `xi=-0.08`：
         - `|detM^P|^2: 3.109e-8 -> 3.108e-8`
         - `|D_mixed^P|^2: 2.122e2 -> 2.122e2`
         - `sigma_total: 10.9661 -> 10.9693`
     - 这些量在 `1e-8~1e-6` 内变化平滑、无“随 `Δs->0` 爆炸增长”的迹象。
  2. **积分离散收敛检查（`n_points=12..96`）**：
     - 以阈值近点 `Δs=1e-8/1e-6/1e-4` 计算 `sigma_total`，在 `xi=-0.12,-0.10,-0.08,0.0,0.2` 上结果稳定到 `~1e-7` 量级以内；
     - 未观察到随积分细化而系统上升（发散）现象。

- 结论（就当前实现与口径）：
  1. 阈值近点确有“分母变小 -> 传播子增强”的**强敏感**，但当前窗口内未出现不可积发散；
  2. 该增强在数值上被积分后保持有限，表现为“可积尖峰/高峰”，而非“积分发散”；
  3. 是否“阈值固有”与“只在个别 `xi` 出现”两者都部分成立：
     - 阈值敏感是 mixed_P 结构固有性质；
     - 但强度是否达到异常级别取决于 `xi`（及对应平衡解）是否把系统推入近零分母的敏感窗口。

- 连续 `xi` 证据（本轮 `xi∈[-0.2,0.5]` 直扫，阈值近点 `Δs=1e-8`）：
  - 最小 `|detM^P|^2` 出现在 `xi≈-0.12`（约 `3.54e-9`），并非仅 `-0.10/-0.08` 两点特例；
  - 离开该窗口后（尤其 `xi>0`），`|detM^P|^2` 单调抬升到 `1e-6` 量级，`|D_mixed^P|^2` 同步回落到个位数，异常敏感性自然减弱。

### 4.10.7 计算链路中的“积分执行顺序”梳理（避免中间量误读）

- 当前链路中至少有四层积分/离散求积，顺序是串联而非并列：
  1. **极化与耦合层（平衡解后）**：`A` 积分 -> `G_u/G_s` -> `K_coeffs`；
  2. **截面层**：固定 `s` 下做 `t` 积分，得到 `sigma(s)`（含 blocking、`M_s/M_t/干涉`）；
  3. **平均散射率层**：在 `(p_i,p_j,cosθ_i,cosθ_j,φ)` 五维求积中使用 `f_i f_j v_rel sigma(s)`，得到 `w_ij`；
  4. **输运层**：`tauinv_i = sum_j n_j w_ij` -> `tau_i=1/tauinv_i` -> 进入输运系数动量积分（`η,σ,ζ`）。
- 因此在解释中间量时要区分：
  - `|D_mixed^P|^2` 与 `sigma(s)` 是上游局部量；
  - `w_ij`/`tauinv` 是下游聚合量，包含分布权重、相对速度、所有通道与核分项的综合结果。

### 4.10.8 回答“`-0.12->-0.10` 为何 `D` area 下降但散射率没同向下降”

- 关键点：此前“`D` area 下降”主要来自 `uubar_to_ddbar` 的 mixed_P 视角；但 `tauinv_u` 用到的是**所有主通道的加权和**，不能由单通道单分项直接代表。

- 直接证据（`u` 物种，来自 `t190_window_channel_diag.csv`）：
  - `tauinv_u: 0.5041 -> 0.5009`（仅小幅下降）；
  - `uubar_to_uubar` 通道贡献 **上升**：`0.1472 -> 0.1534`（rate `0.5758 -> 0.6154`）；
  - `uubar_to_ddbar` 通道贡献 **下降**：`0.1152 -> 0.1139`（rate `0.4507 -> 0.4568`；贡献受密度变化耦合后净降）。

- 与全链路分解对应（`t190_xi_window_full_chain_summary.csv`）：
  - 对 `uubar_to_ddbar`：`area_sigma` 比值约 `0.833`（明显下降）；
  - 对 `uubar_to_uubar`：`area_sigma` 比值约 `0.986`（接近持平）；同时该过程 `M_s` 比值约 `1.031`，部分抵消了 mixed_P 下滑。

- 因而 `-0.12->-0.10` 段不存在真正矛盾：
  1. 某个 `D` 分项（尤其 `mixed_P`）下降，不等于总 `w_ij` 必然同幅下降；
  2. 同一段里其他通道与分项（如 `uubar_to_uubar` 的 `M_s` 变化）可补偿；
   3. 最终 `tauinv_u` 只小幅下降，传到 `sigma_over_T` 就表现为该段回落而非剧烈跃迁。

### 4.11 T190 图中 `tau_s`/`tau_sbar` 在 `xi∈[-0.5,-0.3]` 的异常波动：是否与 `tau_u` 主异常同因

- 你提出的对照问题是关键：`multi_y_tau_u_tau_ubar_tau_s_tau_sbar_vs_xi.png` 里，`tau_s` 与 `tau_sbar` 在 `[-0.5,-0.3]` 的波动是否与前文 `tau_u` 在 `[-0.10,-0.08]` 的异常同一机制造成。

- 先给结论（证据见下）：
  1. **不属于同一“主导机制”**；
  2. 两者都发生在“`sigma(s) -> w_ij -> tauinv -> tau` 串联链路”上，但主导通道与传播子结构不同；
  3. `tau_s/tau_sbar` 该窗口的主异常由 `usbar_to_usbar` 通道贡献突增主导，而不是 `u` 通道案例中的 `mixed_P(detM)` 分母增强主导。

- 复跑口径（与前文一致，新增局部窗口重跑）：
  - 主文件：`D:\Desktop\Temp\relaxtime_t190_window\t190_window_xi_m05_m03_main.csv`
  - 通道诊断：`D:\Desktop\Temp\relaxtime_t190_window\t190_window_xi_m05_m03_channel_diag.csv`
  - 配置保持与前文一致（`T=190, muB=0, finite_15, tau_n_sigma_points=6, sigma_grid_n=128`）。

- 观测层事实（`tau_s=tau_sbar`，`muB=0` 下粒子-反粒子对称）：
  - `xi=-0.50 -> -0.48`: `tau_s 1.294 -> 1.316`（平缓上升）；
  - `xi=-0.48 -> -0.46`: `1.316 -> 1.143`（回落）；
  - `xi=-0.46 -> -0.44`: `1.143 -> 0.540`（显著下跌）；
  - `xi=-0.44 -> -0.42`: `0.540 -> 0.464`（继续下跌，形成低谷）；
  - `xi=-0.42 -> -0.40`: `0.464 -> 0.929`（快速反弹）。
  - 对应 `tauinv_s`：`0.875 -> 1.853 -> 2.154 -> 1.077`（在 `-0.44/-0.42` 区间出现尖峰）。

- 通道层主导项（`species=s`）显示：
  - `usbar_to_usbar` 贡献是绝对主导项，且在异常段出现突增：
    - `xi=-0.46`: `0.571`
    - `xi=-0.44`: `1.555`
    - `xi=-0.42`: `1.862`
    - `xi=-0.40`: `0.791`
  - 同期 `us_to_us` 贡献仅缓慢变化（约 `0.176 -> 0.174 -> 0.170 -> 0.167`），`ssbar_to_ssbar` 也平缓变化（约 `0.057 -> 0.055 -> 0.054 -> 0.052`）。
  - 因此 `tauinv_s` 的尖峰不是多通道同步跃迁，而是 **单主导通道 `usbar_to_usbar` 的局部增强**。

- 密度-率分解进一步确认“不是密度驱动”：
  - `usbar_to_usbar` 的密度因子 `n_u` 在该段单调下降：
    - `0.4341 (-0.46) -> 0.4224 (-0.44) -> 0.4107 (-0.42) -> 0.3986 (-0.40)`；
  - 但其 rate 同期显著上升后回落：
    - `0.658 -> 1.840 -> 2.267 -> 0.993`；
  - 故异常主要由 **rate 核上游变化** 触发，而非密度平移。

- 传播子/振幅结构对照（解释“为何不同于 `tau_u` 案例”）：
  - `usbar_to_usbar` 在映射表里为 `qqbar` 过程，`s` 道为 `K/σ_K` simple，`t` 道有 mixed；
  - 但在异常段对应的分解中，`M_s` 面积占比接近 1（`~99.9%`），`M_t` 仅 `~0.1%` 量级；
  - 且其 `s` 道 `|D_s^P|^2` 由 simple 部分主导、mixed 近零（`area_abs_DsP_mixed ≈ 0`）。
  - 这与 `tau_u` 主异常段（`uubar_to_ddbar/uubar_to_uubar` 里 `mixed_P` 与 `|detM^P|^2` 分母效应主导）**机制上不同**。

- 统一口径总结（避免“同图同因”的误判）：
  1. 两类异常都在同一串联系统中传播（上游振幅/传播子变化会传到 `tau`）；
  2. 但 `tau_u` 主异常可定位到 `mixed_P(detM)` 敏感窗口；
  3. `tau_s/tau_sbar` 在 `[-0.5,-0.3]` 的主要异常由 `usbar_to_usbar` 的 `s` 道 simple（K 类）强度重分配主导；
  4. 因此二者属于“同链路、不同主导子机制”；若把“同构”定义在“**分母接近零导致传播子放大**”这一数学结构层面，则仍可继续检验是否同构（见 4.11.1）。

### 4.11.1 深挖 `s` 道 simple：是否也存在“分母近零 -> 传播子放大”的同构机制

- 你的判断方向是对的：`usbar_to_usbar` 的异常虽然不是 `mixed_P(detM)` 主导，但其主导 simple K 传播子仍是分式结构
  `D_K^P = 2K_{4567}^+ / (1 - 4K_{4567}^+ \Pi_{us}^P)`；
  因此需要直接检查分母 `|1 - 4K\Pi|` 是否在异常段逼近零。

- 扫描结果（`T=190, muB=0`, 过程 `usbar_to_usbar`, `Δs∈[1e-8,2.0]`）：
  - `xi=-0.50`: `min|1-4K\Pi| ≈ 4.30e-2`，`max|D_K^P|^2 ≈ 8.94e1`；
  - `xi=-0.48`: `3.63e-2`，`1.26e2`；
  - `xi=-0.46`: `2.69e-2`，`2.31e2`；
  - `xi=-0.44`: `8.19e-3`，`2.52e3`（窗口内最强放大）；
  - `xi=-0.42`: `2.17e-2`，`3.63e2`；
  - `xi=-0.40`: `5.62e-2`，`5.47e1`。

- 阈值近点直接证据（`Δs=1e-8/1e-6/1e-4`）：
  - `xi=-0.44`：
    - `|1-4K\Pi| ≈ 8.19e-3`（三组 `Δs` 下稳定在同量级）；
    - `|D_K^P|^2 ≈ 2.52e3`；
    - `sigma_total(usbar_to_usbar) ≈ 1.59e2`。
  - `xi=-0.42`：`|1-4K\Pi| ≈ 2.17e-2`，`|D_K^P|^2 ≈ 3.63e2`，`sigma ≈ 2.83e1`。
  - `xi=-0.40`：`|1-4K\Pi| ≈ 5.62e-2`，`|D_K^P|^2 ≈ 5.47e1`，`sigma ≈ 5.29e0`。

- 同构性结论（精确定义）：
  1. 若“同构”指 **小分母触发的非线性放大**，则 `tau_u` 与 `tau_s` 两类异常是同构的；
  2. 但分母对象不同：
     - `tau_u` 主异常：`|D_mixed^P|^2 \propto 1/|detM^P|^2`（mixed）；
     - `tau_s` 主异常：`|D_K^P|^2 \propto 1/|1-4K_{4567}^+\Pi_{us}^P|^2`（simple）；
  3. 因此可归纳为“**同一类数学放大机制在不同传播子分支上的实现**”，而不是“完全不同机理”。

- 与 4.11 的关系修正：
  - 4.11 的“不是同一主导机制”仍成立（主导通道/传播子分支不同）；
  - 但从更高层机制抽象看，二者可统一到“阈值敏感区分母接近零导致的传播子增强”。

- 对你提出的“`abs(x)` + `x` 过零”表述，当前证据支持把机制写得更数学化：
  1. `usbar_to_usbar`（simple K）分母可写作
     `den_simple = 1 - 4K_{4567}^+\Pi_{us}^P = x_s + i y_s`，传播子强度 `|D_K^P|^2 \propto 1/|den_simple|^2`；
  2. 在 `xi∈[-0.5,-0.3]`、阈值最近点 `Δs=1e-8`，`y_s` 始终是 `O(10^-6)` 小量，而 `x_s=Re(den_simple)` 由正到负跨零：
     - `x_s(-0.44)≈+8.19e-3`，`x_s(-0.42)≈-2.17e-2`；
  3. 因 `|den_simple|=sqrt(x_s^2+y_s^2)` 且 `|y_s|<<|x_s|`，在该窗口近似表现为 `|x_s|` 的“V 型”最小值结构（先逼近 0 再远离 0）；
  4. 这正对应 `|D_K^P|^2` 在 `xi≈-0.44` 的峰值放大与两侧快速回落。
  5. 需要严格限定：当前证据支持“`x(\xi)` 在过零邻域可一阶近似线性”，**不支持**把整个窗口写成严格 `y=kx (k>0)` 的全局线性模型。

- `tau_u` 主异常是否也同构（按同一数学视角）——结论是“是，同构但分母对象不同”：
  1. `u` 通道主异常对应 `den_mixed = detM^P = x_u + i y_u`，且 `|D_mixed^P|^2 \propto 1/|den_mixed|^2`；
  2. 在 `xi∈[-0.2,0]`、`Δs=1e-8`，`y_u` 同样远小于 `x_u`，`x_u=Re(detM^P)` 也发生过零：
     - `x_u(-0.12)≈+5.95e-5`，`x_u(-0.10)≈-6.22e-5`，`x_u(-0.08)≈-1.76e-4`；
  3. 因此 `|detM^P|` 也呈“先近零后远离零”的阈值敏感窗口结构，触发 `|D_mixed^P|^2` 的非线性放大/塌缩；
  4. 这与 simple K 的 `|1-4K\Pi|` 机制在数学上同构，但物理承载分支分别是 `mixed_P` 与 `simple K`。

### 4.11.2 局部线性拟合补充：把“`abs(x)` 过零”从定性提升到局部定量

- 为避免把“近零点一阶展开”误写成“全窗口严格线性”，对分母实部做了局部线性拟合（均取 `Δs=1e-8` 阈值最近点）：
  1. simple-K 分支（`process=usbar_to_usbar`, `xi∈[-0.48,-0.40]`）：
     - 拟合 `x_s(\xi)=a_s\xi+b_s` 得 `a_s≈-1.3963`, `b_s≈-0.6104`, `R^2≈0.9921`；
     - 由拟合推得过零点 `\xi_{0,s}=-b_s/a_s≈-0.4372`，与观测到的符号翻转区间 `[-0.44,-0.42]` 一致。
  2. mixed-P 分支（`process=uubar_to_ddbar`, `xi∈[-0.14,-0.06]`）：
     - 拟合 `x_u(\xi)=a_u\xi+b_u` 得 `a_u≈-5.91e-3`, `b_u≈-6.46e-4`, `R^2≈0.9985`；
     - 由拟合推得过零点 `\xi_{0,u}≈-0.1092`，与观测到的符号翻转区间 `[-0.12,-0.10]` 一致。
- 这组拟合支持如下更严谨表述：
  - 在异常窗口内，分母实部 `x(\xi)` 经过零点且虚部次级（`|y|<<|x|`，除极近零邻域）；
  - 因而 `|den|=\sqrt{x^2+y^2}` 呈近 `|x|` 型最小值结构；
  - 在零点邻域以 `x(\xi)\approx a(\xi-\xi_0)` 一阶展开后，可得到与 `|\xi-\xi_0|` 同型的局部行为，并经 `1/|den|^2` 非线性放大到传播子/截面/`tau` 层。

- 复现脚注（本节拟合）
  1. 环境：`julia --project=.`, `T=190`, `muB=0`, 阈值最近点 `Δs=1e-8`。
  2. simple 拟合窗口：`xi=-0.48:0.02:-0.40`，分母 `den_simple=1-4K_{4567}^+\Pi_{us}^P`。
  3. mixed 拟合窗口：`xi=-0.14:0.02:-0.06`，分母 `den_mixed=detM^P`（`uubar_to_ddbar`）。
  4. 复现命令：使用 `scripts/analysis/relaxtime/t190_sigma_chain_decomposition_lib.jl` 中 `build_state_point/process_threshold_info/decompose_mixed_p_propagator_chain` 计算分母，再对 `Re(den)` 做最小二乘直线拟合（本轮拟合命令已在会话记录中保留，可直接复跑）。

### 4.11.3 `y` 的来源判定：异常窗口内“严格为 0 的数值噪声”还是“物理上可非零的小量”

- 结论先行：当前实现与公式共同支持 **后者**——`y=Im(den)` 在该窗口是“物理上允许且本次计算中确为非零的小量”，不是“理论应为 0 但被数值误差污染”。

- 理论/实现依据（为什么“可非零”是模型内生结果）：
  1. `Π` 的底层 `B0` 在代码中显式返回 `(real_part, imag_part)`，且虚部由奇点/切割条件触发（`src/relaxtime/OneLoopIntegrals.jl` 的 `tilde_B0_k_zero` 与 `tilde_B0_k_positive` 都有解析虚部项）；
  2. `polarization_aniso` 直接把 `B0_imag` 乘入 `Π_imag`（`src/relaxtime/PolarizationAniso.jl`），因此 `Π` 一般是复数；
  3. simple 分母 `den_simple=1-4K\Pi` 与 mixed 分母 `detM(M00,M08,M88)` 都由复数 `Π` 传播而来，所以 `Im(den)` 在一般参数点并不要求严格为 0。

- 数值证据（用于排除“纯数值噪声”解释）：
  1. `usbar_to_usbar`（simple）在 `xi=-0.48,-0.46,-0.44,-0.42,-0.40` 的 `Δs=1e-8` 上，`den_im` 分别约
     `-3.19e-6,-3.73e-6,-4.38e-6,-5.17e-6,-6.12e-6`，随 `xi` 平滑单调变化，不呈随机抖动；
  2. 在 `xi=-0.44` 固定点，`Δs: 1e-8 -> 1e-6 -> 1e-4` 时 `den_im` 约
     `-4.38e-6 -> -4.38e-5 -> -4.38e-4`，呈结构化尺度响应；
  3. `uubar_to_ddbar`（mixed）在 `xi=-0.10`，`detM_im` 约
     `-8.04e-8 (Δs=1e-8)`, `-8.04e-7 (Δs=1e-6)`, `-8.04e-6 (Δs=1e-4)`，同样是稳定的系统性变化；
  4. `detM` 改用“仅保留 `Π` 实部”重算后，`|detM|` 仍只发生极小修正（例如 `xi=-0.10`：`6.22231985e-5 -> 6.22231466e-5`），说明虚部是**真实但次级**，而非噪声主导。

- 对“应不应该为 0”的精确口径：
  - 在本模型这条计算链里，`y` 不是“理论上必须严格为 0”的量；
  - 只有在特定运动学/参数点（例如某些切割关闭条件）才可能出现 `Im(Π)=0`；
  - 本异常窗口对应点位下，`y` 处于“非零但小”的区间，因而 `|den|` 的主导仍由 `x` 控制，这与 4.11.1/4.11.2 的 `|x|` 近似并不矛盾。

### 4.11.4 “虚部路径”再核验（按公式条件逐条落位，补强证据）

- 用户关切点是“文档已说明虚部可以严格为 0，那么两类异常区是不是其实都应为 0”。本轮按公式分支把路径显式分成 `k=0` 与 `k>0`，并对异常区逐点核验。

- 先对齐公式判据（来自 `OneLoopIntegral_B0.md` / `OneLoopIntegral_B0_Aniso_Collinear.md`）：
  1. `k=0` 路径：
     - 虚部条件是 `E0 ∈ [m, ΛE]`（等价于阶跃函数 `Θ(ΛE-E0)Θ(E0-m)` 触发）；
     - 若 `E0` 不在区间内，虚部严格为 0。
  2. `k>0` 路径：
     - 虚部条件是角变量极点进入物理区间（可写成 `x_pole∈[-1,1]` 或等价 `Θ(2pk-|...|)`）；
     - 若极点不在区间内，虚部严格为 0。

- 本异常区实际走的是哪条路径：
  - 对本报告两类主异常通道（`usbar_to_usbar` 与 `uubar_to_ddbar`）在阈值最近点 `Δs=1e-8` 的核验显示：`k_s` 全部为 `0.0`，即都走 `k=0` 路径，而不是 `k>0` 路径。
  - 这是由当前分析口径使用 `s` 道质心系动量给出（`calculate_cms_momentum(...,:s,...)` 下 `k=0`）决定的。

- 新增可复现证据（脚本+产物）：
  - 脚本：`scripts/analysis/relaxtime/t190_imag_path_evidence.jl`
  - 明细：`D:\Desktop\Temp\relaxtime_t190_window\t190_imag_path_evidence_detail.csv`
  - 汇总：`D:\Desktop\Temp\relaxtime_t190_window\t190_imag_path_evidence_summary.csv`

- 核验结果 1（`s` 类异常区，simple 分支：`usbar_to_usbar`, `xi=-0.46,-0.44,-0.42,-0.40`）：
  1. `k=0` 且每个 `xi` 的四项 `tilde_B0` 里有 2 项满足 `E0∈[m,ΛE]`（`pole_fraction=0.5`），另外 2 项严格不触发；
  2. 对应 `Π_{us}^P` 虚部非零（约 `4.55e-6 -> 7.36e-6`），并传到 `den_simple=1-4KΠ` 的 `den_im`（约 `-3.73e-6 -> -6.12e-6`）；
  3. 说明该异常区虚部不是噪声，也不是“应为 0 却算成非 0”，而是 `k=0` 分支中 `E0` 命中积分区间后的合法虚部贡献。

- 核验结果 2（`u` 类异常区，mixed 分支：`uubar_to_ddbar`, `xi=-0.12,-0.10,-0.08`）：
  1. `k=0` 且 `Π_{uu}^P` 对应四项中 2 项触发 `E0∈[m,ΛE]`（`pole_fraction=0.5`），因此 `Π_{uu}^P` 虚部非零（约 `1.49e-5 ~ 1.61e-5`）；
  2. `Π_{ss}^P` 对应四项全部不触发（`pole_fraction=0`），因此 `Π_{ss}^P` 虚部严格为 `0`；
  3. 尽管 `Π_{ss}^P` 可严格为 0，mixed 分母 `detM^P=M00*M88-M08^2` 的虚部仍非零（约 `-7.79e-8 ~ -8.27e-8`），来源是 `Π_{uu}^P` 非零虚部经矩阵组合传播到 `detM`。

- 对“是否存在完全为 0 的可能性”的最终口径（与公式文档一致）：
  1. 存在，且本轮已在同一异常分析链内直接观测到：`Π_{ss}^P` 在上述 `u` 类窗口就是严格 0（`k=0` 且未满足 `E0` 区间条件）；
  2. 但“某一子通道虚部为 0”不等于“最终分母虚部必须为 0”；
  3. 两类异常区最终分母虚部（`den_simple_im` 与 `detM_im`）在本口径下都非零，因此前文“`y` 为内生小非零量、不是数值噪声”的结论得到更强、可复现且按公式分支对齐的支持。

- 复现脚注（本节新增）：
  1. `julia --project=. scripts/analysis/relaxtime/t190_imag_path_evidence.jl`
  2. 读取 `t190_imag_path_evidence_summary.csv` 看 `k_min/k_max` 与 `pole_fraction`；
  3. 读取 `t190_imag_path_evidence_detail.csv` 看每个 `term` 的 `has_pole,E0,Emin,Emax,term_im`，可逐项映射到 `k=0` 公式的阶跃条件。

## 5) 结论与置信链条

- 结论 A（问题1）：T150 首点缺失主要由“求解入口未优先触发 `Models.solve` 鲁棒链路 + 无 sidecar 可见性”共同导致；该问题已被修复。
  - 证据链：
    1. 代码层：`solve_models_equilibrium(...)` 已改为 `Models.solve` 优先、`solve_constraint` 回退。
    2. 测试层：`tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl` 覆盖入口优先路径与 sidecar 字段。
    3. 运行层：`T=150,muB=0,xi=-0.5` 最小复跑成功，主 CSV 存在结果行，failed sidecar 无失败记录。
  - 置信度：High（代码-测试-运行三层闭环）。

- 结论 B（问题2，限定为 T190 `xi∈[-0.2,0]` 的 `sigma_over_T` 异常段）：异常折点可追溯到
  **`xi -> 平衡解(序参量/质量) -> sigma(s)阈值邻域收缩 -> 通道贡献突变 -> tauinv 跃迁 -> tau 跃迁 -> sigma/T 折点`**，
  且主导变化发生在 `uubar_to_ddbar` / `uubar_to_uubar` 两条通道。
  - 证据链（基于窗口重跑 `D:\Desktop\Temp\relaxtime_t190_window\`）：
    1. 观测层：`sigma/T` 在 `xi=-0.10 -> -0.08` 发生主要跃升（约 `+2.10e-3`），是图中“异常感”核心。
    2. tau 层：同区间 `tau_u` 同步大幅增加（约 `+0.373`），`tauinv_u` 同步下降（约 `-0.0788`）。
    3. 通道层（u 物种主导项，按降幅排序）：`uubar_to_ddbar` 贡献约 `0.1139 -> 0.0714`（占总降幅约 `53.9%`），`uubar_to_uubar` 约 `0.1534 -> 0.1247`（约 `36.4%`）；两者合计约 `90.3%`，因此被选为主分析通道。
    4. 率与密度分解：该突变同时包含
       - 轻度密度变化（`n_u≈n_ubar` 从约 `0.2493` 到 `0.2434`），
       - 以及更显著的通道率变化（例如 `uubar_to_ddbar` 率约 `0.4568 -> 0.2935`，`uubar_to_uubar` 率约 `0.6154 -> 0.5124`）。
    5. 物理背景层：`m_u,m_s` 在该窗口单调平滑上升，无同幅突跳；因此异常主要是“通道率重分配+比值敏感放大”而非相变不连续。
    6. `sigma` 分项层：反事实结果显示 `t` 区间宽度几乎不变、blocking 仅次级影响，相比之下 `K_coeffs` 变化对两主通道 `sigma` 面积有显著影响（约 18%~31% 量级），支持“振幅/传播子链路主导”的判断。
     7. 振幅层：`|M|^2` 分解显示两主通道均由 `s` 道主导（约 86%~96%）；A→B 时 `s` 道总量显著下降。
     8. s道根因层（公式对齐）：`M_s` 下降主要由 `|D_s^P|^2` 下降驱动，而 `s12^-*s34^-` 运动学因子并未同步下降（反而略升），因此“异常毛刺”的主导来源定位到 **P通道传播子（`K^+`/`Π^P`）链路**。
     9. P通道绝对强度层（新增）：`|D_s^P|^2` 的绝对面积在 A→B 显著下移（`uubar_to_ddbar`: `71.106->20.633`; `uubar_to_uubar`: `92.086->38.428`），且下移主要来自 `mixed_P` 分量（`79.870->27.905`），`simple` 分量仅小幅变化（`1.726->1.626`）。
       10. mixed_P 深层公式分解（新增）：在 `|D_mixed^P|^2 = 4(detK^+)^2|J^TMJ'|^2/|detM^P|^2` 中，主导变化来自 `|detM^P|^2` 上移（`B/A≈1.126`）；`detK^+` 近乎不变（`≈0.998`），`|J^TMJ'|^2` 仅轻微下降（`≈0.990`）。进一步把 `detM^P=M00*M88-M08^2` 拆为三项后可见：`|M00*M88|^2` 仅小幅上升（`≈1.026`），而 `|M08^2|^2`（`≈1.121`）与交叉项 `-2Re(M00*M88*(M08^2)^*)`（`≈1.184`）抬升更强，构成分母增强主导来源。再下钻到 `Π/K` 项后可见 `Π_uu^P` 相关项在 `M08/M88` 中系统性上升（约 `≈1.014`），而 `K08^+` 本身下降（`≈0.966`），因此主导方向来自“极化响应组合”而非“耦合常数单独放大”。
      11. 连续 `xi` 直扫校验（新增）：`|detM^P|^2` 在 `-0.10 -> -0.08` 的阈值最近点出现约 `8.03x` 上升，同时 `|D_mixed^P|^2` 同点约 `0.123x` 下跌，验证“近阈值小分母敏感放大”确实发生在异常窗口，而非仅 A/B 两点偶然。
      12. 全窗口深拆校验（本轮新增）：`xi∈[-0.2,0]` 的逐点分解显示，小回跳段（`-0.14 -> -0.12`）与主跃升段（`-0.10 -> -0.08`）都可由同一 mixed_P/detM 链路解释，只是相邻段响应方向与幅度不同；并且与 `xi>0` 对照可明确区分“敏感窗口非单调”与“正区间缓变单调”。
      13. 新增“同构机制”统一表述：把分母写成 `x+iy` 后，在异常窗口内均表现为 `|y|<<|x|` 且存在“`|x|` 局部最小”的小分母结构（可表现为过零，也可表现为同号近零逼近），从而通过 `1/|den|^2` 放大传播子；`tau_u` 用 `den=detM^P`，`tau_s` 用 `den=1-4K_{4567}^+\Pi_{us}^P`。
   - 置信度：High（已完成“主导通道识别 -> 阈值邻域 -> 5D核分项 -> sigma分项反事实 -> 振幅分项(s/t/干涉) -> s道传播子/运动学再分解”的逐层闭环；通道选择依据为降幅排序而非先验指定）。

- 问题B收尾判定（当前版本）：
  - 可以收尾到“机制定位”层：`mixed_P` 异常主因已定位为 `detM^P` 分母增强，且证据链可复现；
  - 不建议继续把 `Π` 做更深的一阶项以下拆解作为当前阻塞项（除非目标改为“解析灵敏度建模/参数重整”），因为现有证据已足以解释异常来源与传导方向。

- 结论 C（T200）：存在与 T190 同类的比值敏感放大现象，但本轮深挖证据链重点在 T190 异常窗口；对 T200 目前置信度低于 T190。
  - 置信度：Medium（范围统计充分，通道级链路尚未做同深度重跑分解）。

## 6) 风险与后续

- 风险 1：当前 sidecar 记录的是失败可见性，不自动给出“可修复建议”；后续排障仍需结合日志与点位重跑。
- 风险 2：比值类指标对离散步长、插值与阈值邻域更敏感，跨版本对比时应固定网格与积分配置。
- 建议后续：
  1. 在固定 `T=190/200` 增加更细 `xi` 局部重采样（仅敏感区）验证毛刺是否收敛。
  2. 对 T200 复制本次 T190 的通道级诊断流程（含 `--channel-diagnostics-output`），补齐对称证据链。
  3. 在 `sigma(s)` 阈值邻域输出关键通道（`uubar_to_ddbar`,`uubar_to_uubar`）的采样点与插值残差，验证“率突降”是否来源于阈值邻域积分/插值敏感区。
