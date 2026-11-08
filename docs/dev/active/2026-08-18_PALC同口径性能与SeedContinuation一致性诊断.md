# PALC 同口径性能与 SeedContinuation 一致性诊断

状态：`research / review`。2026-11-08 已完成两窗口 Actions pilot、计数开关数值 parity 和结果复核，见 [PR #338](https://github.com/w5851/Julia_RelaxTime/pull/338)。建议保留隔离诊断能力；候选根选择敏感性在进入生产 backend 前另行处理。完整结果见[分析](../../analysis/pnjl/algorithmic_feasibility/palc_matched_20261108_analysis/README.md)。

## 1. 背景与目标

PALC multi-branch spike 已在两个 `FixedAsymmetricRho` 多根场景中展示分支发现、延拓和显式 pressure selection 能力，但现有 `125.4 s cold / 1.37 s warm` 计时来自不同 workload，且缺少 matched dense/seed baseline 与 PNJL solver telemetry。因此当前只能保留 `ready_for_opt_in_diagnostic_backend`，不能宣称 PALC 加速或替代 production path。

本任务只负责建立可比较的性能和语义 parity 证据，决定 PALC 是否值得继续作为可选诊断 backend。

## 2. 范围与非范围

### 范围

- 对同一 `FixedAsymmetricRho` 点集比较 PALC、当前 pressure-governed `solve_multi`/SeedContinuation 和必要的 dense reference。
- 分离 cold/warm 启动成本，固定模型、节点、容差、seed policy、branch governance 和后处理口径。
- 记录 anchor solves、continuation branch points、residual/Jacobian 调用、Newton iterations、fallback/rescue、失败尝试和 runner-minutes。
- 验证 `SeedContinuation` 在小 rho grid 上与当前 `trho_asymmetric` production-like scan 的结果语义一致。

### 非范围

- 不把 BifurcationKit 加入根 `Project.toml`。
- 不修改 `Models.solve`、`Models.solve_multi`、`run_trho_scan` 默认语义或当前 C0/C1/C2 配置。
- 不在本机运行 PNJL solver；正式数值诊断只允许通过明确授权的 GitHub Actions pilot。
- 不降低物理/数值容差，不改变 phase branch governance，不启动 C3/O1 或 reference promotion。

## 3. 本阶段结果与后续决策

- 两窗口 `SeedContinuation` 与组合生产脚本 equilibrium helper 的结果一致。
- T=120 的原步长 SeedContinuation 在 rho=0.375 选中低压力根；PALC、独立逐点 `solve_multi` 与半步长参考一致选中高压力根。各方法重复结果稳定，残差均通过；该差异属于候选可见性/路径依赖，不应靠放宽容差处理。
- T=130 首轮独立逐点 `solve_multi` 在 rho=0.75 选中低压力根，完整重跑五方法均为同一高压力根；跨运行敏感性反例已保留，具体环境/舍入影响未定位。
- 首轮 T=130 PALC 依赖预编译超时已修复；完整重跑两个窗口成功，PALC 与半步长参考均 5/5 一致。本地未运行 PNJL solver。
- PALC warm 中位耗时 0.818 / 0.380 s，SeedContinuation 1.874 / 1.115 s；130 MeV 同结果成本比 2.93。120 MeV 分支不同，不给生产替换意义的加速比。
- 本阶段计算与验证完成，待审阅诊断结论。下一阶段可聚焦内层多种子根被最小残差选择提前舍弃的问题，不把该改进混入本批生产默认。

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
- [ ] 作者审核后才能决定是否创建 opt-in PALC backend PR；本任务本身不升级 production。

## 6. 风险与回退

- BifurcationKit 编译和初始化成本可能主导短测；回退为同进程 warm 对照并单独报告 startup。
- 多分支 anchor discovery 可能让 PALC 工作量天然高于单分支 scan；回退为 diagnostic oracle，不替换当前路径。
- parity 或 branch selection 漂移时，保留当前 pressure-governed production scan，不修改现有 C0/C1/C2 产物。
