---
title: sysimage 产品化 Phase C：package 化 C3 第四批 closure/tag 形状治理任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3第四批closure-tag形状治理任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C3 第四批 closure/tag 形状治理任务单

更新日期：2026-05-04

当前状态：已完成首轮 closure/tag 形状治理验证；已确认收益可传导到主 sysimage 真实稳定路径。

> 目的：把 `solver_ad` 残余从“可被更窄 workload 吃掉的 ForwardDiff 杂项”继续推进到“`residual!#140 / _gap_conditions_with_model` 的 closure/tag 形状治理”，优先尝试通过更稳定的 helper/调用形状减少 nested-tag 分裂，而不是继续扩大 precompile capability 集。

---

## 1. 背景

- [x] C3 第三批已完成 solver-only probe 与 residual-only probe
- [x] `solver_residual_ad` capability 已接入 `scan/core/full/all`
- [x] 主 sysimage 复测后，focused residual 未继续下降
- [x] 当前稳定残余主干为：
  - `residual_140 = 26`
  - `conditions_gap = 6`
  - `thermo_probe_256_257 = 4`
  - `thermo_probe_258_259 = 4`

结论：

- [x] 继续“加 workload / capability”已不是高收益方向
- [x] 下一步应直接检查 solver closure 捕获形状、`Base.Pairs`、`solver_backend`、以及 `Conditions` 内部闭包/tag 嵌套来源

---

## 2. 目标

- [x] C3-11 定位 `residual!#140` 与 `#129#130` 的真实源码锚点
- [x] C3-12 盘点 solver residual 闭包当前捕获了哪些 shape-sensitive 对象
- [x] C3-13 评估是否能用更稳定 helper 降低 `Base.Pairs` / kwcall / nested-tag 分裂
- [x] C3-14 在不改动物理语义与 unified solver 合同的前提下，完成一轮 closure/tag 形状治理试验

---

## 3. 本批范围

### 3.1 包含

- [x] `src/models/` 中 solver 入口与 `Conditions` 内部 helper 的形状治理
- [ ] 如有必要，新增仅服务于 solver residual 稳定化的 script-private perf/trace 脚本
- [x] 对 focused trace family 的前后对比

### 3.2 不包含

- [ ] 不改动物理模型、收敛判据、容差、种子治理语义
- [ ] 不回滚 C2/C3 前批已经验证有效的 capability / workflow 收敛
- [ ] 不处理与当前活跃物理开发冲突的 phase / physics 主线改动

---

## 4. 拟实施步骤

### 4.1 找到 residual/tag 的源码锚点

- [x] 精确定位：
  - `Models.var"#residual!#140"`
  - `Models.var"#129#130"`
  - `Models.Conditions.var"#_gap_conditions_with_model..."`
- [x] 判断这些匿名函数分别对应：
  - `solve(...)` 入口 residual
  - `solve_constraint(FixedMu)` 入口 residual
  - `Conditions` 内部 pressure/gradient 闭包

验收：

- [x] 能把 trace family 逐一映射到具体源码位置与调用层级

### 4.2 盘点 shape-sensitive 捕获项

- [x] 检查 solver residual 闭包当前是否捕获：
  - `Base.Pairs`
  - `solver_backend`
  - `model`
  - `T_fm / mu_fm`
  - `GapParams`
  - `mu_vec`
- [x] 判断哪些捕获项会导致 tag / kwcall specialization 分裂

验收：

- [x] 输出“哪些是必要捕获，哪些可下沉到更稳定 helper 参数”的判断

### 4.3 试做形状治理

- [x] 优先考虑：
  - 把 residual 闭包改为更稳定的具名 helper 调用
  - 减少 `Base.Pairs` 直接进入 `ForwardDiff.Tag` 相关签名
  - 减少 `Conditions` 内部 pressure closure 的额外捕获层次
- [x] 保持 `Models.solve / solve_constraint / build_residual!` 对外语义不变

验收：

- [x] 至少完成一轮治理前后 focused trace diff

### 4.4 回写结论

- [x] 更新本批活动文档
- [x] 更新总任务单 C3 进度
- [x] 明确后续是继续细化 closure/tag，还是转入 Phase D / 其他任务

---

## 5. 风险与判断标准

### 5.1 风险

- [x] 风险 R1：匿名函数编号不稳定，trace family 与源码锚点映射可能需间接确认
- [x] 风险 R2：若 helper 抽取得过头，可能意外触碰 unified solver 合同
- [x] 风险 R3：残余中有一部分就是 ForwardDiff 运行时 shape 派生，无法完全消除

### 5.2 判断标准

- [x] 优先接受 focused trace family 的明确收敛
- [x] 若 residual 未下降，但能把“不可再靠 precompile 收敛”的原因说明白，也算有效产出

---

## 6. 验证计划

- [ ] `julia --project=. scripts/dev/check_precompile_profile_coverage.jl`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] 针对 solver focused workload 的 trace diff 至少一轮

---

## 7. DoD

- [x] 已把主要 residual tag family 映射到具体源码锚点
- [x] 已说明当前 closure/tag 分裂的关键来源
- [x] 已完成至少一轮形状治理试验与对比
- [x] 已给出下一步是否继续深挖 solver tag 主干的判断

---

## 8. 首轮实施结果

### 8.1 代码侧形状治理

- [x] `src/models/solver/runtime/GapSolver.jl`
  - 以具名 callable helper 替代匿名 residual / omega closure：
    - `GapResidualInplace`
    - `GapOmegaFn`
    - `_materialize_gap_kwargs`
- [x] `src/models/solver/spec/Conditions.jl`
  - 以具名 callable helper 替代匿名 pressure closure：
    - `GapPressureFn`
    - `_gap_gradient`
- [x] 保持 unified solver 合同不变，仅治理 AD tag/closure 形状

### 8.2 solver-focused probe 结果

- [x] `build/solver_ad_probe/JuliaRelaxTimeSolverAdProbe.dll`
- [x] `build/trace/solver_ad_solver_probe_shape_governance_summary.json`
- [x] 结果：
  - `no-sys`: `496 -> 77`
  - `with-sys`: `142 -> 41`
  - 相比治理前的 probe with-sys focus `42`，本轮下降到 `41`

### 8.3 主 sysimage 真实稳定路径复测

- [x] `julia --project=. scripts/perf/relaxtime/measure_solver_ad_sysimage_residual.jl .\\build\\JuliaRelaxTime.dll main_shape_governance`
  - `build/trace/solver_ad_main_shape_governance_summary.json`
  - `with-sys focus_trace_lines: 41`
  - 对比上一轮 `main_after_solver_residual` 的 `42`，主 sysimage 同步下降 `1`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\\build\\JuliaRelaxTime.dll main_shape_governance`
  - `build/trace/transport_point_main_shape_governance_summary.json`
  - summary 中 `focus_trace_lines` 仍为 `20`，该口径是 transport-point 专用 focus lens，比较的是 no-sys / with-sys 同一脚本内的命中情况，不等于 with-sys residual common 总量
  - 真实 residual 变化需看 with-sys trace diff；对比上一轮 `transport_point_main_after_solver_residual_with_sys.jl`，focused residual 从 `87` 降到 `86`

### 8.4 新旧 with-sys trace diff 结论

- [x] solver 主路径 diff：
  - `focus_a=42`
  - `focus_b=41`
  - 旧匿名 family（`residual!#140` / `_gap_conditions_with_model...` / `#129#130`）已整体切换为具名 helper family
  - 新 family 以 `GapResidualInplace / GapOmegaFn / GapPressureFn / OnceDifferentiable{GapResidualInplace,...}` 为主
- [x] transport_point 主路径 diff：
  - `focus_a=87`
  - `focus_b=86`
  - `solver_ad` 桶从 `40` 降到 `39`
  - 说明 closure/tag 形状治理不只在 solver-focused probe 中生效，也传导到了真实 transport 稳定路径

### 8.5 当前判断

- [x] 本批证明：从“继续加 workload/capability”转向“closure/tag 形状治理”是有效方向
- [x] 但单批收益目前仍是边际压降（主 sysimage 两条真实路径都只下降 `1`）
- [x] 下一步若继续留在 C3，更合适的方向是：
  - 继续细化 `config_loader` 与 transport workflow 残余桶
  - 或继续做更窄的 solver/thermo helper 形状治理
- [x] 若当前 package 化目标已满足，也可以把 C3 记为“已证明可继续，但收益转入精修”，转入 Phase D 的 launcher/app 化准备
