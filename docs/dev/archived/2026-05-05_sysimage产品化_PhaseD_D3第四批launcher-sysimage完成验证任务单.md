---
title: sysimage 产品化 Phase D：D3 第四批 launcher sysimage 完成验证任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D3第四批launcher-sysimage完成验证任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase D：D3 第四批 launcher sysimage 完成验证任务单

更新日期：2026-05-05

当前状态：已完成收口；已将构建参数压到明显偏验证导向的轻量档，并确认在 10 分钟窗口内仍停留在 `create sysimage start` 之后，尚未拿到 `sys.dll` / `jrt-launcher(.exe)` 落盘证据。该批现作为 D3 主收口的最终判据文档，而非继续追短窗口重跑的进行中文档。

> 目的：把 D3 当前剩余问题严格收敛为“完成条件验证”，不再回到入口边界、包壳结构或命名契约讨论。

---

## 1. 前置结论

- [x] D3 首批已完成 launcher package 壳与 `phase` 解释器链路验证
- [x] D3 第二批已确认结构性失败来自 launcher env 的 `Pkg` 依赖/manifest，同步后已修复
- [x] D3 第三批已把构建路线收敛到“手工 bundle + repo sysimage 作为 base sysimage 的极薄 launcher sysimage”
- [x] 当前剩余问题只在 launcher sysimage 的完成时长与最终产物落盘证据

---

## 2. 本批目标

- [x] D3c-1 明确 launcher sysimage 构建的最小完成判据
- [x] D3c-2 选定一组更轻的本地 PoC 编译参数
- [x] D3c-3 以可观测方式运行一次较长窗口验证
- [ ] D3c-4 拿到以下至少一项硬证据：
  - `build/launcher-app/lib/julia/sys.dll` 落盘
  - `build/launcher-app/bin/jrt-launcher(.exe)` 落盘
- [x] D3c-5 把结果同步回 D3 主任务单与 backlog

---

## 3. 本批范围

- [ ] `scripts/dev/build_launcher_app.jl`
- [ ] `scripts/dev/diagnose_launcher_app_build.jl`
- [ ] `build/launcher-app/` 产物验证
- [ ] D3 文档状态回写

非范围：

- [ ] 不新增 launcher 子命令
- [ ] 不再重做 package 壳结构
- [ ] 不推进 installer / release 分发

---

## 4. 关键策略

- [x] 优先通过“更轻编译参数”降低 launcher sysimage 完成时间
- [ ] 保留 repo sysimage 作为 base sysimage 的主路线不变
- [x] 用日志/阶段标记确保能区分“仍在编译”与“已走到产物生成”

---

## 5. 验收标准

- [x] 已给出 launcher sysimage 构建的明确完成判据
- [x] 已完成一轮较长窗口验证
- [x] 已明确本地 PoC 路线是否能在可接受窗口内产出 `sys.dll` / `jrt-launcher(.exe)`

---

## 6. DoD

- [ ] 已拿到首份 launcher app 最终产物证据，或
- [x] 已能明确证明当前参数下仍不可接受，并给出下一步唯一合理的收口方向

---

## 7. 完成判据

- [x] `build/launcher-app/lib/julia/sys.dll` 落盘，视为 launcher sysimage 已完成
- [x] `build/launcher-app/bin/jrt-launcher(.exe)` 落盘，视为可执行包装层已生成
- [x] 若日志长期停留在 `create sysimage start` 且以上两项均未出现，则判定当前参数下仍不可接受

## 8. 本批轻量参数验证

- [x] 已将 [scripts/dev/build_launcher_app.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/dev/build_launcher_app.jl) 与 [scripts/dev/diagnose_launcher_app_build.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/dev/diagnose_launcher_app_build.jl) 调整为：
  - `base_sysimage = build/JuliaRelaxTime.dll`
  - `include_transitive_dependencies = false`
  - `sysimage_build_args = \`-O0 --compile=min --strip-metadata --strip-ir\``
- [x] 轻量参数版日志已证明：
  - 前序 `bundle_*` 步骤可在秒级完成
  - 日志可稳定推进到 `create sysimage start`
  - 10 分钟窗口内仍未出现 `create sysimage done`
  - 同时 `build/launcher-app-diagnose/lib/julia/` 下仍未看到 `sys.dll`
- [x] 因此当前结论是：
  - D3 结构、命名契约、launcher env 依赖已基本排障完毕
  - 当前剩余阻塞是 launcher sysimage 编译 wall time 过长
  - 继续在本地短窗口内反复尝试已接近低收益

## 9. 下一步建议

- [x] 下一步不宜继续在当前线程做更多短窗口本地重跑
- [x] 更合理的收口方向二选一：
  - 方向 A：接受 D3 当前定位为“仓库内 launcher PoC + 产物链未在本机短窗口收口”，转入后续阶段
  - 方向 B：单独开一个长时构建/离线验证批次，只做一次完整等待，目标是拿最终 `sys.dll` / `jrt-launcher(.exe)` 证据
  - 2026-05-05 已采纳方向 A：D3 主线按“仓库内 launcher PoC 已成立，最终可执行产物长时验证另列”收口；若后续仍需要最终产物证据，应另开独立长时批次
