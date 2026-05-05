---
title: sysimage 产品化 Phase C：package 化 C3 第五批 config-defaults 与冷启动指标任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseC_package化_C3第五批config-defaults与冷启动指标任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C3 第五批 config-defaults 与冷启动指标任务单

更新日期：2026-05-05

当前状态：已完成首轮实现与复测；本批主要产出是“冷启动 wall time 主指标化”与一轮 `TransportWorkflow` config-defaults 缓存收敛试验。

> 目的：在不触碰物理语义与 unified solver 合同的前提下，继续收窄 transport workflow 侧的冷启动开销，并把 C3 的收益表达从“trace 数量”提升为“冷启动 wall time 为主、trace 为辅”。

---

## 1. 本批目标

- [x] C3-15 把冷启动耗时纳入 sysimage perf 对比主指标
- [x] C3-16 评估 `TransportWorkflow` 默认 TOML 配置的重复加载是否还能收敛
- [x] C3-17 对当前主 sysimage 真实路径给出 wall time + trace 的联合判断

---

## 2. 改动范围

### 2.1 代码

- [x] [src/models/workflow_apps/TransportWorkflow.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/models/workflow_apps/TransportWorkflow.jl)
  - 新增 `transport_workflow_defaults_cache`
  - 将 `prefer_energy_aniso` 与 `a_builder_config` 的默认 TOML 读取收敛为一次 `load_config(...)`
  - 保持原有 `_default_prefer_energy_aniso_from_toml` / `_default_a_builder_config_from_toml` 对外调用形状不变
- [x] [src/models/pnjl/workflows/TransportWorkflow.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/models/pnjl/workflows/TransportWorkflow.jl)
  - 同步收敛历史兼容实现，避免双轨偏移
- [x] [tests/unit/models/test_transport_workflow.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/tests/unit/models/test_transport_workflow.jl)
  - 补充新 cache 字段的 reset / populated 校验

### 2.2 工具

- [x] [scripts/perf/relaxtime/compare_sysimage_summary_walltime.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/perf/relaxtime/compare_sysimage_summary_walltime.jl)
  - 对两份 sysimage summary JSON 输出：
    - `no-sys wall_ms`
    - `with-sys wall_ms`
    - `improvement wall_ms`
    - 对应 trace/focus 变化

---

## 3. 验证

- [x] `julia --project=. -e 'include("tests/unit/models/test_transport_workflow.jl")'`
- [x] `julia --project=. -e 'include("tests/integration/relaxtime/test_transport_workflow_toml_prefer_energy_aniso_smoke.jl")'`
- [x] `julia --project=. scripts/dev/build_sysimage.jl`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\build\JuliaRelaxTime.dll main_config_defaults_cache`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl build\trace\transport_point_main_shape_governance_with_sys.jl build\trace\transport_point_main_config_defaults_cache_with_sys.jl`
- [x] `julia --project=. scripts/perf/relaxtime/compare_sysimage_summary_walltime.jl build\trace\transport_point_main_shape_governance_summary.json build\trace\transport_point_main_config_defaults_cache_summary.json`

---

## 4. 结果

### 4.1 冷启动 wall time

- [x] 当前复测得到：
  - `no-sys`: `51954.7 ms`
  - `with-sys`: `35318.9 ms`
  - `improvement`: `16635.7 ms`

### 4.2 相比上一轮 `main_shape_governance`

- [x] wall time 对比脚本输出：
  - `no-sys: 82747.0 -> 51954.7 ms`
  - `with-sys: 45637.2 -> 35318.9 ms`
  - `improvement: 37109.7 -> 16635.7 ms`

### 4.3 trace / residual

- [x] `transport_point_main_shape_governance_with_sys` vs `transport_point_main_config_defaults_cache_with_sys`
  - `focus_a=86`
  - `focus_b=86`
  - `common=86`
  - `only_a=0`
  - `only_b=0`
- [x] 结论：当前 `transport_point` 主 workload 的 with-sys focused residual 完全未变

---

## 5. 判断

- [x] `TransportWorkflow` config-defaults 收敛本身是合理的工程整理：
  - 默认配置从“两次 `load_config`”收成“一次 `load_config` + 两处字段复用”
  - 对 profile cache/reset 语义无回归
- [x] 但对当前 `transport_point` 主 workload 而言，本批**没有形成可归因的 residual 收敛**
- [x] 当前观察到的 wall time 下降不能直接归因到本批代码改动，因为：
  - 新旧 with-sys focused trace 完全一致
  - 冷启动绝对时间会受宿主 package cache / 依赖预编译状态影响
- [x] 因此，从证据强度排序：
  1. `wall_ms` 更直观，适合作为主展示指标
  2. `trace diff` 仍然是“这次改动是否真的击中 residual family”的归因证据

---

## 6. 下一步建议

- [x] C3 若继续推进，优先选择“当前主 workload 能实际命中的残余家族”：
  - transport kwargs 分流/清洗相关匿名迭代器与 generator family
  - `transport_api` 桶而非未命中的 config-defaults 分支
- [x] 之后所有 C3 perf 汇报，默认采用：
  - 主指标：`no-sys / with-sys / improvement` 三组 `wall_ms`
  - 辅指标：`trace_lines / focus_trace_lines`
  - 归因：必要时补 `with-sys trace diff`
