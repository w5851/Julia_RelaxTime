---
title: sysimage 产品化 Phase C：package 化 C3 第二批 transport-point 收益测量任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3第二批transport-point收益测量任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C3 第二批 transport-point 收益测量任务单

更新日期：2026-05-04

当前状态：第二批已完成；已补充 focused probe sysimage 与 transport-point residual/timing 测量脚本，并完成 residual trace diff 与 `ConfigLoader` 单实例复用验证。

> 目的：不只停留在“新增 capability 已接入 precompile registry”的结构结论，而是进一步用 trace / wall-time 测量验证 `transport_point_api` 对冷启动残余编译的实际影响。

---

## 1. 目标

- [x] C3-4 为 transport-point 热路径补 focused sysimage probe
- [x] C3-5 固定测量 workload，比较 `no-sys / with-sys`
- [x] C3-6 判断 `transport_point_api` 的收益是“显著新增”还是“显式化既有覆盖”

---

## 2. 本批范围

### 2.1 包含

- [x] 新增 `scripts/dev/build_transport_point_probe_sysimage.jl`
- [x] 新增 `scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl`
- [x] 产出 `build/trace/transport_point_*_summary.json`

### 2.2 不包含

- [x] 不修改现有主 sysimage release workflow
- [x] 不在本批强行重构 trace budget 门槛

---

## 3. 实现结果

### 3.1 新测量能力

- [x] focused probe sysimage 构建脚本
- [x] transport-point workload 的 residual trace + wall-time 测量脚本
- [x] 支持指定 sysimage 路径，便于比较主 sysimage 与 probe sysimage

### 3.2 实测结论

#### 3.2.1 focused probe sysimage

测量文件：

- [x] `build/trace/transport_point_probe_sysimage_summary.json`

结果：

- [x] `no-sys`: `155649.6 ms`, residual trace `578`
- [x] `with probe sysimage`: `131477.6 ms`, residual trace `195`
- [x] 改善：`24172.1 ms`，residual trace 减少 `383`

#### 3.2.2 当前主 sysimage 复测

测量文件：

- [x] `build/trace/transport_point_main_sysimage_summary.json`

结果：

- [x] `no-sys`: `137180.0 ms`, residual trace `578`
- [x] `with main sysimage`: `93271.7 ms`, residual trace `195`
- [x] 改善：`43908.3 ms`，residual trace 减少 `383`

### 3.3 阶段判断

- [x] `transport_point_api` 对 transport-point workload 的 `no-sys -> with-sys` 路径有明确收益
- [x] 但从 residual trace 角度看，focused probe sysimage 与当前主 sysimage 都落在 `195` 行，未显示出新的 residual-trace 压降

解释：

- 当前 capability 更像是把已有的 transport-point 热路径覆盖显式化、稳定化
- 仅从本轮 evidence 看，还不能证明它相对现有主 sysimage 带来了额外的残余编译压缩

#### 3.2.3 `ConfigLoader` 单实例复用后复测

测量文件：

- [x] `build/trace/transport_point_probe_sharedcfg_summary.json`

结果：

- [x] `no-sys`: `139893.5 ms`, residual trace `568`
- [x] `with probe sysimage`: `113725.7 ms`, residual trace `185`
- [x] 改善：`26167.8 ms`，residual trace 减少 `383`
- [x] focused summary 口径从 `32` 降到 `20`

### 3.4 residual focus lines 细分定位

新增脚本：

- [x] `scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl`

细分 diff 结论：

- [x] `main_with_sys` 与 `probe_with_sys` 的 focused residual 完全相同
  - `focus_a=97`, `focus_b=97`, `common=97`, `only_a=0`, `only_b=0`
- [x] 说明 `transport_point_api` 没有继续压低这批 focused residual

按桶统计，当前 with-sys 共同 residual 主要分为：

- [x] `solver_ad = 43`
- [x] `config_loader = 30`
- [x] `transport_api = 20`
- [x] `quadrature = 4`

与 `no_sys` 相比，sysimage 已压掉的 focused extra residual 主要是：

- [x] `solver_ad = 37`
- [x] `config_loader = 3`

定位判断：

- [x] `solver_ad` 残留主要不是简单 capability 缺失，而是 `Models.solve` / `_solve_constraint_fixedmu` / ForwardDiff dual/partials / `kwcall` 这一层高阶 AD specialization 仍在运行时落下 residual
- [x] `config_loader` 残留集中在 `Constants_PNJL.ConfigLoader` 与 `Models.TransportWorkflow.ConfigLoader` 的 canonical-json / generator / lock 路径
- [x] `transport_api` 残留集中在 `Models.solve_transport_from_equilibrium` 与 `TransportIntegrationConfig` / provider 装配相关 wrapper

因此，“那 32 行为什么没有继续下降”的更准确表述是：

- [x] 原先粗粒度 focus 口径低估了 residual 的细分规模；细化后实际共同 residual 是 97 行
- [x] 这些 residual 并不是 probe capability 新旧不同导致的差异项，而是主 sysimage 和 probe sysimage 共同未完全吃掉的一层 runtime specialization / config-loading / wrapper residual

### 3.5 `ConfigLoader` 双实例根因与修复结果

- [x] 已确认 `Constants_PNJL` 与 `Models.TransportWorkflow` 分别本地 `include(src/config/ConfigLoader.jl)`，导致 trace 中同时出现 `Constants_PNJL.ConfigLoader` 与 `Models.TransportWorkflow.ConfigLoader`
- [x] 已将 `TransportWorkflow` 改为复用 `Main.Constants_PNJL.load_config`
- [x] `main_with_sys` vs `probe_sharedcfg_with_sys`：
  - `common=85`
  - `only_main=12`，且全部落在旧 `Models.TransportWorkflow.ConfigLoader` residual
  - `only_probe=2`，转为 `Constants_PNJL.ConfigLoader` 的布尔 canonical-json 路径

修复后 with-sys 共同 residual 桶分布：

- [x] `solver_ad = 43`
- [x] `config_loader = 18`
- [x] `transport_api = 20`
- [x] `quadrature = 4`

结论：

- [x] 这一批的实质收益来自共享 helper / 单实例模块复用，不是继续追加 capability 名称
- [x] `config_loader` 桶已从此前的 `30` 压到 `18`
- [x] 下一批优先级应转向 `solver_ad`

---

## 4. 风险与备注

- [x] 尝试全量重建主 sysimage 时出现超时；中断后本地 `build/JuliaRelaxTime.dll` 一度变成 0 字节
- [x] 已从已有 release zip 恢复 `build/JuliaRelaxTime.dll` 与 metadata，避免本地稳定入口继续指向损坏产物

说明：

- focused probe sysimage 结论可信，可用于判断 capability 本身是否有价值
- 若要评估“主 sysimage 产品化产物”的精确增量收益，后续仍需要一次成功的全量主 sysimage 重建

---

## 5. 验证

- [x] `julia --project=. scripts/dev/build_transport_point_probe_sysimage.jl`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\build\transport_point_probe\JuliaRelaxTimeTransportPointProbe.dll`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\build\JuliaRelaxTime.dll`
- [x] `julia --project=. scripts/perf/relaxtime/measure_transport_point_sysimage_residual.jl .\build\transport_point_probe\JuliaRelaxTimeTransportPointProbe.dll probe_sharedcfg`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl build\trace\transport_point_probe_sharedcfg_no_sys.jl build\trace\transport_point_probe_sharedcfg_with_sys.jl`
- [x] `julia --project=. scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl build\trace\transport_point_main_with_sys.jl build\trace\transport_point_probe_sharedcfg_with_sys.jl`
- [x] `julia --project=. scripts/dev/check_precompile_profile_coverage.jl`

---

## 6. 下一步建议

- [x] 若继续推进 C3，应优先针对 residual focus lines 做更细的 trace diff，而不是继续机械增加 capability 名称
- [ ] 下一批优先考虑：
  - 先针对 `solver_ad` 桶做更窄、更稳定的 ForwardDiff / kwcall precompile probe
  - 再判断 `transport_api` wrapper residual 是否值得继续收敛
- [ ] 若要把 evidence 闭环到产品化主 sysimage，需要安排一次可控的全量重建窗口

---

## 7. DoD

- [x] 已有 focused probe sysimage 与对应测量脚本
- [x] 已完成至少一轮 `no-sys / with-sys` 实测
- [x] 已明确当前 capability 的收益边界与下一步重点
