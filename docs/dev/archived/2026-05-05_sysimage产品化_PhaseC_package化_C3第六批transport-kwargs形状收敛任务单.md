---
title: sysimage 产品化 Phase C：package 化 C3 第六批 transport-kwargs 形状收敛任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseC_package化_C3第六批transport-kwargs形状收敛任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C3 第六批 transport-kwargs 形状收敛任务单

更新日期：2026-05-05

当前状态：已完成首轮实现、验证与主 sysimage 复测；本批确认命中了 `transport_point` 主 workload 的真实 residual family。

> 目的：针对 `TransportWorkflow.solve_transport_from_equilibrium(...)` 内 `transport_kwargs` 的匿名 `Filter/Generator/merge` 路径做 shape 收敛，在不改动物理语义的前提下进一步压缩 `transport_api` 桶残余，并继续以冷启动 wall time 为主展示指标、以 trace diff 为归因证据。

---

## 1. 目标

- [x] C3-18 定位 `transport_api` 桶中仍未收敛的 kwargs/generator family
- [x] C3-19 用更稳定的 materialize helper 替代匿名 `Filter/Generator/merge`
- [x] C3-20 复测主 sysimage 的 `transport_point` 冷启动与 with-sys trace diff

---

## 2. 代码改动

- [x] [src/models/workflow_apps/TransportWorkflow.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/models/workflow_apps/TransportWorkflow.jl)
  - 新增：
    - `_select_namedtuple_keys`
    - `_drop_namedtuple_keys`
    - `_finalize_transport_kwargs`
  - 将：
    - `_extract_transport_integration_kwargs`
    - `_drop_transport_integration_keys`
    - `_drop_transport_provider_keys`
    从匿名 `pairs(kwargs)` 过滤链收敛为显式 helper
  - `solve_transport_from_equilibrium(...)` 改为统一走 `_finalize_transport_kwargs(...)`
- [x] [src/models/pnjl/workflows/TransportWorkflow.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/models/pnjl/workflows/TransportWorkflow.jl)
  - 同步历史兼容实现，避免主线与兼容层偏移
- [x] [tests/unit/models/test_transport_workflow.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/tests/unit/models/test_transport_workflow.jl)
  - 新增 helper 级断言：
    - integration kwargs 提取
    - provider 保留/剔除
    - extra kwargs 保留

---

## 3. 验证

- [x] `julia --project=. -e 'include("tests/unit/models/test_transport_workflow.jl")'`
- [x] `julia --project=. -e 'include("tests/integration/relaxtime/test_transport_workflow_toml_prefer_energy_aniso_smoke.jl")'`
- [x] `julia --project=. scripts/dev/build_sysimage.jl`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\build\JuliaRelaxTime.dll main_transport_kwargs_materialize`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl build\trace\transport_point_main_config_defaults_cache_with_sys.jl build\trace\transport_point_main_transport_kwargs_materialize_with_sys.jl`
- [x] `julia --project=. scripts/perf/relaxtime/compare_sysimage_summary_walltime.jl build\trace\transport_point_main_config_defaults_cache_summary.json build\trace\transport_point_main_transport_kwargs_materialize_summary.json`

---

## 4. 结果

### 4.1 冷启动 wall time

- [x] 本轮 summary：
  - `no-sys`: `69507.9 ms`
  - `with-sys`: `43000.9 ms`
  - `improvement`: `26507.0 ms`

### 4.2 相比上一轮 `main_config_defaults_cache`

- [x] wall time 对比：
  - `no-sys: 51954.7 -> 69507.9 ms`
  - `with-sys: 35318.9 -> 43000.9 ms`
  - `improvement: 16635.7 -> 26507.0 ms`

说明：

- [x] 单次 wall time 绝对值出现波动，不能单独作为“本批代码收益”的充分证据
- [x] 本批更可信的归因证据来自 with-sys trace diff

### 4.3 with-sys trace diff

- [x] `transport_point_main_config_defaults_cache_with_sys` vs `transport_point_main_transport_kwargs_materialize_with_sys`
  - `focus_a=86`
  - `focus_b=84`
  - `common=74`
  - `only_a=12`
  - `only_b=10`

- [x] 旧 family 被消掉的主要是：
  - `Base.Generator{Base.Iterators.Filter{...}}`
  - `Base.merge(..., Generator(Filter(...)))`
  - `_extract/_drop transport *_kwargs` 对应匿名 family

- [x] 新 family 仍有少量 helper 相关残余，但总量已下降 `2`

### 4.4 transport-point focus lens

- [x] `focus_trace_lines`：
  - `no-sys: 20 -> 18`
  - `with-sys: 20 -> 18`

这说明本批不只是替换了 trace family 名称，而是确实减少了 `transport_point` 主 workload 当前 focus lens 命中的编译项。

---

## 5. 判断

- [x] 本批是 C3 中一次“真正命中 transport 主路径 residual family”的收敛
- [x] 这次命中的不是 solver_ad，而是 `transport_api` 桶中的 kwargs/materialize 路径
- [x] `transport_api` 桶仍未清零，但已从匿名筛选链进一步收窄
- [x] 后续若继续留在 C3，应优先沿着新出现的 helper family 继续细化，而不是回到未命中的 config-defaults 分支

---

## 6. 指标口径

- [x] 默认主指标：`wall_ms`
  - `no-sys`
  - `with-sys`
  - `improvement`
- [x] 默认辅指标：
  - `trace_lines`
  - `focus_trace_lines`
- [x] 默认归因证据：
  - `with-sys trace diff`

补充：

- [x] 对于单次冷启动 wall time，后续若需要做批间强结论，建议增加重复测量/中位数口径，避免宿主 package cache 波动误导判断。
