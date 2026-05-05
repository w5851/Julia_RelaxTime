---
title: sysimage 产品化 Phase C：package 化 C3 第三批 solver_ad 残余收敛任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3第三批solver-ad残余收敛任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C3 第三批 solver_ad 残余收敛任务单

更新日期：2026-05-04

当前状态：进行中；solver-focused trace/probe 脚本已落地，并完成首轮 solver-only probe 测量。

> 目的：把当前 transport-point focused residual 中最顽固的 `solver_ad` 桶单独抽离出来，围绕 `Models.solve` 内部的 `ForwardDiff.Tag`、`DifferentiationInterface.prepare_jacobian`、`NLSolversBase.OnceDifferentiable` 等高阶 AD specialization 做更窄的 trace / probe / precompile 收敛，而不是继续依赖宽 capability 叠加。

---

## 1. 背景

- [x] C3 第二批已经确认：`transport_point_api` 本身没有继续压低 `solver_ad` residual
- [x] `ConfigLoader` 双实例问题已修复，`config_loader` 桶已从 `30` 压到 `18`
- [x] 当前 with-sys focused residual 仍稳定包含：
  - `solver_ad = 43`
  - `transport_api = 20`
  - `config_loader = 18`
  - `quadrature = 4`
- [x] 已验证把 `gap_solver_ad` 对齐到 `Models.solve(...)`，以及在 probe sysimage 中叠加 `gap_solver_ad + ad_shape_stabilization + transport_point_api`，并没有进一步压低 `solver_ad` residual

结论：

- [x] 下一步需要从“profile/capability 级别预热”转向“solver AD 具体签名级别预热”

---

## 2. 目标

- [ ] C3-7 为 `solver_ad` 残余建立独立 trace focus 口径
- [ ] C3-8 设计只瞄准 `Models.solve` / `_solve_constraint_fixedmu` AD 链的窄 probe
- [x] C3-9 评估哪些 residual 属于可被 precompile 收敛的稳定签名，哪些属于 runtime 形状派生
- [x] C3-10 在不改动物理语义的前提下，争取把 `solver_ad` 桶从当前 `43` 继续压低

---

## 3. 本批范围

### 3.1 包含

- [x] 新增或扩展 `scripts/perf/relaxtime/` 下的 solver-focused trace 分析脚本
- [x] 新增或扩展 `scripts/dev/` 下的 solver-focused probe sysimage 构建脚本
- [x] 在 `src/models/precompile/registry.jl` 中补更窄的 solver AD probe（若 evidence 支持）
- [x] 产出 `build/trace/solver_ad_*` 或等价命名的 residual summary / diff 证据

### 3.2 不包含

- [ ] 不改动 `Models.solve` 的物理语义、收敛判据或 unified solver 合同
- [ ] 不调整数值容差来制造“性能改善”
- [ ] 不在本批处理 `transport_api` wrapper 或 `quadrature` 桶，除非它们阻塞 `solver_ad` 归因
- [ ] 不回滚现有 `transport_point_api` / `ConfigLoader` 收敛成果

---

## 4. 拟实施步骤

### 4.1 先把 `solver_ad` 残余从混合 focus 里单独切出来

- [x] 细化当前 trace diff 脚本，增加只匹配以下对象的 focus 口径：
  - `ForwardDiff.Tag`
  - `ForwardDiff.Dual`
  - `ForwardDiff.Partials`
  - `DifferentiationInterface.prepare_jacobian`
  - `NLSolversBase.OnceDifferentiable`
  - `Models.solve`
  - `Models._solve_constraint_fixedmu`
  - `Models.Conditions`
- [x] 输出按签名族聚类的桶分布，而不是只给总量

验收：

- [x] 能回答“43 条里哪几族最多、是否主要集中在某一两个 kwcall / jacobian prep / dual conversion 路径”

首轮结论：

- [x] 当前 with-sys solver focus 共 `42` 条
- [x] 共同 residual 桶分布：
  - `forwarddiff_tag = 38`
  - `prepare_jacobian = 2`
  - `solve_kwcall = 2`
- [x] `no-sys` 相比 `with-sys` 多出的 `37` 条主要是：
  - `forwarddiff_misc = 30`
  - `forwarddiff_partials = 5`
  - `once_differentiable = 1`
  - `solve_kwcall = 1`

判断：

- [x] `solver_ad` 桶中“容易被 probe 吃掉”的部分主要已经是 `ForwardDiff misc / partials / OnceDifferentiable` 这层
- [x] 目前最顽固的残余已经收缩到 `ForwardDiff.Tag` 主导，而不是泛泛的 AD 噪声

### 4.2 建立更窄的 solver probe

- [x] 从 transport-point workload 中剥离一个只覆盖 `Models.solve(...)` 的最小 workload
- [x] 尽量固定：
  - `FixedMu()`
  - `p_num/t_num`
  - `iterations`
  - `seed_guess`
  - `xi`
  - `residual_norm_max`
- [x] 避免混入 transport 后处理、config loader、quadrature 之外的噪声路径

验收：

- [x] 能单独比较 `solver-only no-sys / with-sys` trace residual

已产出：

- [x] `scripts/dev/build_solver_ad_probe_sysimage.jl`
- [x] `scripts/perf/relaxtime/measure_solver_ad_sysimage_residual.jl`
- [x] `scripts/perf/relaxtime/analyze_solver_ad_trace_focus.jl`
- [x] `build/solver_ad_probe/JuliaRelaxTimeSolverAdProbe.dll`
- [x] `build/trace/solver_ad_solver_probe_summary.json`

首轮测量：

- [x] `no-sys`: `67176.2 ms`, trace `498`, solver focus `79`
- [x] `with-sys`: `55893.5 ms`, trace `143`, solver focus `42`
- [x] 改善：`11282.7 ms`，trace 减少 `355`，solver focus 减少 `37`

### 4.3 针对稳定签名做预热试验

- [ ] 判断是否需要把下列对象显式纳入 probe：
  - `Models.solve(...)` 的 `Core.kwcall`
  - `Models._solve_constraint_fixedmu(...)` 的 `Core.kwcall`
  - `DifferentiationInterface.prepare_jacobian(...)`
  - `NLSolversBase.OnceDifferentiable(..., ForwardDiff.Chunk{5})`
  - 关键 `ForwardDiff.Dual / Partials` 构造签名
- [ ] 只在 evidence 支持时再决定是否抽成新的 capability 或合并进既有 capability

验收：

- [ ] 至少完成一轮“加某组窄 probe 前后”的 residual diff

当前判断：

- [x] 现有 solver-only probe 已经证明：专门的 solver workload 可以稳定压掉 `37` 条 solver extra residual
- [x] 下一步应继续细分 `ForwardDiff.Tag` 主导的 `42` 条共同 residual，而不是先把 solver probe 合并回 registry capability

新增收敛证据：

- [x] residual-only probe 已落地：
  - `scripts/dev/build_solver_residual_probe_sysimage.jl`
  - `scripts/perf/relaxtime/measure_solver_residual_sysimage_residual.jl`
- [x] `build/solver_residual_probe/JuliaRelaxTimeSolverResidualProbe.dll` 已构建
- [x] `build/trace/solver_residual_residual_probe_summary.json` 已产出
- [x] residual-only 首轮测量：
  - `no-sys`: `43863.3 ms`, trace `428`, focus `65`
  - `with-sys`: `41297.6 ms`, trace `114`, focus `33`
  - 改善：`2565.7 ms`，trace 减少 `314`，focus 减少 `32`
- [x] residual-only diff 结论：
  - `common = 33`, `only_a = 32`, `only_b = 0`
  - `common` 主要为：
    - `forwarddiff_tag = 29`
    - `prepare_jacobian = 2`
    - `solve_kwcall = 1`
    - `forwarddiff_misc = 1`
  - `only_a` 主要为：
    - `forwarddiff_misc = 30`
    - `forwarddiff_partials = 2`

判断更新：

- [x] 更窄的 `Conditions.build_residual! + ForwardDiff.jacobian` probe 可以稳定吃掉 `ForwardDiff` 杂项和 partials
- [x] 剩余共同 residual 已进一步集中到 `Models.Conditions.var"#gap_conditions..."` 与外层 jacobian tag 嵌套
- [x] 基于该证据，已把 `solver_residual_ad` 独立 capability 接入 `src/models/precompile/registry.jl` 的 `scan/core/full/all` profile，用于验证对主 sysimage focused residual 的实际收益
- [x] 主 sysimage 复测已完成，但 focused residual 未继续下降：
  - `solver_ad_main_after_solver_residual_with_sys.jl` 与 `solver_ad_solver_probe_with_sys.jl` 逐行一致
  - solver-focused `with-sys` 仍为 `42`
  - tag family 仍为：
    - `residual_140 = 26`
    - `conditions_gap = 6`
    - `thermo_probe_256_257 = 4`
    - `thermo_probe_258_259 = 4`
- [x] transport-focused `with-sys` residual 也未变化：
  - `transport_point_main_after_solver_residual_with_sys.jl` 与 `transport_point_probe_sharedcfg_with_sys.jl` 逐行一致
  - focused common 仍为 `87`
  - 桶分布仍为：
    - `solver_ad = 43`
    - `transport_api = 20`
    - `config_loader = 20`
    - `quadrature = 4`
- [x] 当前结论：`solver_residual_ad` 能稳定收敛“额外的 ForwardDiff 杂项/partials”，但不能继续压低主 sysimage 的剩余 nested-tag 主干

### 4.4 回写结论

- [x] 更新本批活动文档
- [x] 更新总任务单的 C3 进度
- [x] 明确下一步是继续打 `solver_ad`，还是切回 `transport_api`

---

## 5. 风险与判断标准

### 5.1 风险

- [ ] 风险 R1：`solver_ad` residual 中部分是 ForwardDiff / ADTypes / NLSolversBase 的运行时派生签名，难以完全通过静态 precompile 吃掉
- [ ] 风险 R2：为了压 residual 而把 probe 写得过宽，重新混入 transport-point 其它噪声路径
- [ ] 风险 R3：probe 如果依赖过深内部细节，后续 package 化边界变化时容易失效

### 5.2 判断标准

- [ ] 优先接受“trace residual 明显下降”的改动
- [ ] wall-time 只在排除本地 package cache 干扰后才作为主要证据
- [ ] 若 residual 不下降，但能把“不可收敛 residual”与“可收敛 residual”边界说清，也算有效产出

---

## 6. 验证计划

- [ ] `julia --project=. scripts/dev/check_precompile_profile_coverage.jl`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [ ] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`
- [x] solver-focused trace / probe 脚本实跑至少一轮

已完成：

- [x] `julia --project=. scripts/perf/relaxtime/analyze_solver_ad_trace_focus.jl build\trace\transport_point_probe_sharedcfg_no_sys.jl build\trace\transport_point_probe_sharedcfg_with_sys.jl`
- [x] `julia --project=. scripts/dev/build_solver_ad_probe_sysimage.jl`
- [x] `julia --project=. scripts/perf/relaxtime/measure_solver_ad_sysimage_residual.jl .\build\solver_ad_probe\JuliaRelaxTimeSolverAdProbe.dll solver_probe`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_solver_ad_trace_focus.jl build\trace\solver_ad_solver_probe_no_sys.jl build\trace\solver_ad_solver_probe_with_sys.jl`
- [x] `julia --project=. scripts/dev/build_solver_residual_probe_sysimage.jl`
- [x] `julia --project=. scripts/perf/relaxtime/measure_solver_residual_sysimage_residual.jl .\build\solver_residual_probe\JuliaRelaxTimeSolverResidualProbe.dll residual_probe`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_solver_ad_trace_focus.jl build\trace\solver_residual_residual_probe_no_sys.jl build\trace\solver_residual_residual_probe_with_sys.jl`
- [x] `julia --project=. scripts/dev/check_precompile_profile_coverage.jl`
- [x] `julia --project=. scripts/dev/build_sysimage.jl`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\build\JuliaRelaxTime.dll main_after_solver_residual`
- [x] `julia --project=. scripts/perf/relaxtime/measure_solver_ad_sysimage_residual.jl .\build\JuliaRelaxTime.dll main_after_solver_residual`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_solver_ad_trace_focus.jl build\trace\solver_ad_solver_probe_with_sys.jl build\trace\solver_ad_main_after_solver_residual_with_sys.jl`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl build\trace\transport_point_probe_sharedcfg_with_sys.jl build\trace\transport_point_main_after_solver_residual_with_sys.jl`

---

## 7. DoD

- [x] 已有独立的 `solver_ad` focused trace/probe 方案
- [x] 已明确 `ForwardDiff.Tag / prepare_jacobian / OnceDifferentiable` 中哪些是主 residual 来源
- [x] 已完成至少一轮窄 probe 前后对比
- [x] 已给出下一步是继续压 `solver_ad` 还是转回 `transport_api` 的判断
