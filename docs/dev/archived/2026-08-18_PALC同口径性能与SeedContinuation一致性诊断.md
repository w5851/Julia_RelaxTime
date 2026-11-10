---
title: "PALC 同口径性能与 SeedContinuation 一致性诊断"
archived: true
original: "docs/dev/active/2026-08-18_PALC同口径性能与SeedContinuation一致性诊断.md"
archived_date: "2026-11-10"
task_status: "completed"
---

以下为原始内容（保留，以便审阅与历史参考）：

---

# PALC 同口径性能与 SeedContinuation 一致性诊断

状态：`research / completed`。[PR #340](https://github.com/w5851/Julia_RelaxTime/pull/340) 的全部 13 项检查通过，2026-11-10 已 squash 合并到 `main`（`26671f2a`），对应本地和远程分支已回收。作者接受前两阶段并授权通用 PALC 接入与约一小时 Actions 对照；两批累计 61.10 runner-minutes。较长 PNJL Tρ 默认启用 PALC，Tμ 默认保持 legacy。具体计时、分支差异和未来可选的验证范围见[分析](../../analysis/pnjl/algorithmic_feasibility/palc_matched_20261108_analysis/README.md)。

## 1. 背景与目标

本任务始于历史 PALC 与 SeedContinuation 计时口径不一致。现已补齐同点、冷暖调用和真实求解工作量比较；新结果支持较长同密度路径上的延拓收益，短窗口需要考虑 anchor 成本。

本任务建立可比较的性能和语义 parity 证据。通用路径现已接入 `run_tmu_scan`、`run_trho_scan`、相图生产入口与 CLI；生产对照支持对标准 PNJL tensor 积分、默认 selector、至少 21 点的 Tρ 路径启用 PALC。Tμ 旧策略完整对照未在预算内完成，默认保持 legacy；未再追加预算计算。

## 2. 范围与非范围

### 范围

- 对同一 `FixedAsymmetricRho` 点集比较 PALC、当前 pressure-governed `solve_multi`/SeedContinuation 和必要的 dense reference。
- 分离 cold/warm 启动成本，固定模型、节点、容差、seed policy、branch governance 和后处理口径。
- 记录 anchor solves、continuation branch points、residual/Jacobian 调用、Newton iterations、fallback/rescue、失败尝试和 runner-minutes。
- 验证 `SeedContinuation` 在小 rho grid 上与当前 `trho_asymmetric` production-like scan 的结果语义一致。
- 第二阶段增加显式可选的单次候选遍历，验证各合格根均可到达选择器、种子顺序与路径方向的敏感性；固定密度下同时报告化学势与自由能。

### 非范围

- BifurcationKit 0.8.0 作为通用 PALC 的根环境依赖；扫描默认策略依据新实验结果决定。
- `Models.solve`、`Models.solve_multi` 的单点合同和当前 C0/C1/C2 产物保持既有状态。
- 不在本机运行 PNJL solver；正式数值诊断只允许通过明确授权的 GitHub Actions pilot。
- 不降低物理/数值容差，不改变 phase branch governance，不启动 C3/O1 或 reference promotion。

## 3. 本阶段结果与后续决策

- 两窗口 `SeedContinuation` 与组合生产脚本 equilibrium helper 的结果一致。
- T=120 的原步长 SeedContinuation 在 rho=0.375 选中低压力根；PALC、独立逐点 `solve_multi` 与半步长参考一致选中高压力根。各方法重复结果稳定，残差均通过；该差异属于候选可见性/路径依赖，不应靠放宽容差处理。
- T=130 首轮独立逐点 `solve_multi` 在 rho=0.75 选中低压力根，完整重跑五方法均为同一高压力根；跨运行敏感性反例已保留，具体环境/舍入影响未定位。
- 首轮 T=130 PALC 依赖预编译超时已修复；完整重跑两个窗口成功，PALC 与半步长参考均 5/5 一致。本地未运行 PNJL solver。
- PALC warm 中位耗时 0.818 / 0.380 s，SeedContinuation 1.874 / 1.115 s；130 MeV 同结果成本比 2.93。120 MeV 分支不同，不给生产替换意义的加速比。
- 首阶段计算与验证已获作者认可；下一阶段聚焦内层多种子根被最小残差选择提前舍弃的问题，用可选算法与现有路径比较，再决定后续接入。

## 首阶段工作量（2026-11-08）

- 复用 `T=120 MeV, rho/rho0=0.30:0.025:0.40` 和 `T=130 MeV, rho/rho0=0.75:0.025:0.85`，每窗口五个目标点；`xi=0, u/d=0.876, s_target=0, p_num=8, t_num=4`。
- 所有方法采用相同积分节点、`residual_norm_max=1e-6` 和最大压力选择；记录当前内部求解容差，最终用相同显式残差、物理条件和密度误差检验。该 pilot 比较离散方程求解，不据此声称积分收敛。
- 对比 PALC（包含 anchor 发现、双向延拓及目标 rho 校正）、SeedContinuation、逐点 `solve_multi`、生产脚本 equilibrium helper，以及半步长 SeedContinuation 参考；参考的额外计算成本单独显示。
- 每种方法在独立 Julia 进程加载相同依赖后，执行一次首次调用和两次相同工作量的 warm 调用；加载耗时单列，首次调用仍包含尚未发生的 JIT。所有重复重新发现 anchor/求解，不复用前一轮结果。
- FixedAsymmetricRho 当前在外层候选治理内再次遍历种子，必须同时统计外层请求和实际内层 NLsolve；PALC 计入有限差分 Jacobian 内部残差及失败 Newton 尝试。
- 原 spike 的最近 rho 节点对比只作 seed 来源；性能 pilot 的结果比较必须先校正到目标 rho，且将校正成本纳入 PALC。
- 生产式 parity 复用 `combined_meson_density_equilibrium.jl`，覆盖 downstream 介子步骤成功时的 equilibrium 序列；介子步骤失败时生产脚本保留旧 seed 的行为不在本 pilot 模拟。
- 同点 parity 门槛为 natural-unit 状态距离 `1e-4`、压力差 `1e-6 fm^-4`，与 `1e-6` 显式残差和密度约束同时检查；另输出连续差值供复核。所有方法固定 tensor quadrature，adaptive-quadrature 计数不适用。

## 4. 诊断任务

- [x] 设计固定 workload：相同 `T/rho/xi`、`p_num/t_num`、quadrature policy、iterations、容差、seed 和 pressure selection。
- [x] 每方法进程分别测首次调用/warm，加载和依赖预编译单列；首次调用含剩余 JIT。
- [x] 通过 Actions 采集 residual/Jacobian/Newton/attempt/fallback/nonconverged telemetry；固定 tensor quadrature 不涉及自适应积分。
- [x] 对照 PALC、SeedContinuation/`solve_multi` 和半步长参考的结果、分支数、压力选择和残差。
- [x] 采用停止条件：已有分支差异时保持 isolated diagnostic，不提升生产默认；后续改进先解决候选根保留。

## 5. 验收标准与 DoD

- [x] 同工作量配置、源码 SHA 和两窗口结果汇总保存在 `data/outputs/results/analysis/palc_matched_20261108/pilot.json`。
- [x] 每个 backend 有首次调用/warm wall time、runner-minutes、anchor/branch/continuation 计数和 solver telemetry。
- [x] `SeedContinuation` 与组合扫描 equilibrium helper 通过两窗口 parity；生产算法和其它模式语义未改。
- [x] 结果区分分支发现、离散方程精度和性能，首轮预编译失败与候选根反例均已记录。
- [x] 作者接受首阶段结果并授权下一阶段候选搜索优化和有限验证；生产默认切换须以新结果为依据。

## 第二阶段验收

- [x] 单次候选遍历保留全部合格根，统一交给既有 selector；21 点 NLsolve 次数从 2526 降至 476。
- [x] 两处反例通过种子排列、扫描方向与独立进程核对；低压力根的质量层级失败明确保留。T120 跨区化学势回弯未被平滑删除。
- [x] 完成 local/long/higher 各两批及 T120 跨区单批；同点、半步长和正反向均一致。
- [x] 24×8 下 single-pass Seed 约快 7 倍；PALC 再有约 1.9–4.4 倍长路径收益。保留可选 single-pass 和隔离 PALC 能力，暂不接入通用生产 backend。
- [x] 第二阶段数值归并至 `single_pass_pilot.json`，当前状态在台账与 PR 同步。

## 6. 风险与回退

- BifurcationKit 编译和初始化成本可能主导短测；回退为同进程 warm 对照并单独报告 startup。
- 多分支 anchor discovery 可能让 PALC 工作量天然高于单分支 scan；回退为 diagnostic oracle，不替换当前路径。
- parity 或 branch selection 漂移时，保留当前 pressure-governed production scan，不修改现有 C0/C1/C2 产物。

## 第三阶段：生产入口对照

- [x] 三路线统一路径与候选校正；固定密度保持八变量联立约束，全部合格根进入原有最大压力选择。
- [x] 五温度、两 ξ、23 个密度或 31 个 μq 的中等网格；24×8 tensor 节点，iterations=300，同一残差/物理判据。
- [x] 两类 Tρ 的 PALC/single_pass 各首次+warm，三策略共 230 点同点对照；FixedRho 全部吻合，旧 asym 的五个质量层级失败点保留为差异证据。
- [x] Tμ 的 PALC 两批各首次+warm，single_pass 完整首次 310 点一致；single_pass warm 在 232/310 点触及预算，legacy 未开始。因此 Tμ 默认不晋升，缺项列为未来按需扩展。
- [x] T120、ξ=0 的回弯保留，FixedRho Maxwell 三策略一致；production 相图单切片三策略均完成。
- [x] 两批累计 61.10 runner-minutes 后结束实验；一份 `production_comparison.json` 保留代表路径、耗时、实际误差和不同根，重复 CSV 不入库。
- [x] API、CLI/TOML、CSV/callback、完整文件 resume 验证完成；新增默认入口和部分续扫 CI 检查。
- [x] PR #340 全部 13 项检查通过并已合并；integration core 323/323，含 T120 默认 PALC、零密度、回弯/Maxwell 和部分续扫；对应本地和远程分支已回收。

比较的是离散方程求解与扫描算法；这组中等网格不承担积分收敛或高精度相图参考晋升。
