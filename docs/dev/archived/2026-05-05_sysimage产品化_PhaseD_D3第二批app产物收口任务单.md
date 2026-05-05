---
title: sysimage 产品化 Phase D：D3 第二批 app 产物收口任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D3第二批app产物收口任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase D：D3 第二批 app 产物收口任务单

更新日期：2026-05-05

当前状态：已完成收口；结构性失败已修复，且该批结论已被后续第三/第四批验证固化。当前保留结论是：D3 的真实剩余阻塞位于 `create_sysimage(...)` 的长耗时编译阶段，而不是 launcher env / manifest / executables 契约。

> 目的：把 D3 首批中“解释器内 launcher PoC 已通，但 app bundle 最终 launcher 产物未收口”的剩余问题单独拉平，避免继续混在主 PoC 任务单里。

---

## 1. 前置事实

- [x] D3 首批已完成 `phase_cli_support` 复用层
- [x] D3 首批已完成最小 `app/launcher/` package 壳
- [x] D3 首批已完成 launcher `help` / `phase --help` / `phase smoke` 解释器链路验证
- [x] `scripts/dev/build_launcher_app.jl` 已能开始产出 `build/launcher-app/` bundle 目录
- [x] 当前缺口只剩：自定义 launcher 可执行文件未在验证窗口内落盘

---

## 2. 本批目标

- [x] D3b-1 定位 `create_app(...)` 实际停留阶段
- [x] D3b-2 判断问题属于“纯耗时过长”还是“构建策略/产物生成异常”
- [x] D3b-3 给出最小收口方案：
  - 要么让自定义 launcher 可执行文件真正落盘
  - 要么把 D3 路线进一步明确拆成“仓库内 app / 半独立 app”
- [x] D3b-4 把结论同步回 D3 主任务单与 backlog

---

## 3. 排查范围

- [ ] `scripts/dev/build_launcher_app.jl`
- [ ] `app/launcher/`
- [ ] 本地 `PackageCompiler.create_app(...)` 行为
- [ ] `build/launcher-app/` 实际产物布局

非范围：

- [ ] 本批不扩到 `transport-scan` / `server`
- [ ] 本批不回头改动物理求解主线
- [ ] 本批不直接推进 installer / auto-update

---

## 4. 关键问题

### 4.1 停留点定位

- [x] `create_app(...)` 是否卡在 `create_sysimage(...)`
- [x] 是否已经走到 `create_executable_from_sysimg(...)`
- [x] `executables=["jrt-launcher" => "julia_main"]` 是否被正确采纳

### 4.2 产物策略判断

- [x] 若问题只是耗时过长，判断是否通过更薄 app env / 更小 precompile 面即可收口
- [ ] 若问题与当前 launcher 运行时 `Pkg.activate(repo_root)` 路线冲突，明确是否应把 D3 当前 PoC 改称“仓库内 app”而非“app 分发 PoC”

---

## 5. 验收标准

- [x] 已有明确证据说明 `create_app(...)` 卡在哪一步
- [x] 已有下一步最小收口方案，而不是继续盲目长时间重跑
- [x] D3 主任务单与 backlog 已同步新的阶段判断

---

## 6. DoD

- [x] 本批已把 D3 剩余问题收束成可执行结论
- [x] 后续若还需下一批，应已明确是“产物生成问题”还是“产品定位问题”

---

## 7. 本批诊断结论

- [x] 本地 `PackageCompiler.create_app(...)` 关键路径已读：`bundle_* -> create_sysimage(...) -> create_executable_from_sysimg(...)`
- [x] 通过拆步脚本 `scripts/dev/diagnose_launcher_app_build.jl` 已确认：
  - `bundle_julia_libraries / bundle_julia_libexec / bundle_julia_executable / bundle_artifacts / bundle_project / bundle_preferences / bundle_cert` 均能正常完成
  - 当前停留点在 `create sysimage start` 之后
  - 说明自定义 launcher 可执行文件未落盘，并不是 `executables` 参数失效，而是还没有走到 `create_executable_from_sysimg(...)`
- [x] 首个结构性失败已定位并修复：
  - 失败原因：`JRTLauncherApp` 使用了 `Pkg`，但 [app/launcher/Project.toml](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/app/launcher/Project.toml) 未声明 `Pkg` 依赖，且 launcher manifest 未同步
  - 已修复：补 `Pkg` 依赖，并在构建/诊断脚本中先 `Pkg.resolve()` 再 `Pkg.instantiate()`
- [x] 当前策略判断：
  - D3 本地 PoC 构建脚本已改为 `incremental=true`
  - 这更符合当前“先证明 app 产物链可收口”的 PoC 目标，而不是追求首轮就做完全独立、从零构建的 app
- [x] 第三批已进一步把构建脚本切到“手工 bundle + 以仓库现有 `build/JuliaRelaxTime.dll` 作为 base sysimage 编译极薄 launcher sysimage”的路线
  - 结论：这一路线已能稳定推进到 `build/launcher-app/lib/julia/` 目标路径创建阶段，但当前验证窗口内仍停留在长耗时 sysimage 编译，尚未拿到 `sys.dll` 与 `jrt-launcher(.exe)` 落盘证据

## 8. 下一步建议

- [x] 下一批若继续 D3，应只做一件事：围绕 `create_sysimage(...)` 的 wall time 与完成条件做更长窗口或更轻参数验证，目标是拿到“sysimage 完成”与“launcher exe 落盘”的第一份硬证据
  - 2026-05-05 D3 第四批已执行该方向；结论是即使切到 `base_sysimage + include_transitive_dependencies=false + -O0 --compile=min --strip-metadata --strip-ir` 的轻量参数组合，10 分钟窗口内仍停留在 `create sysimage start` 之后，尚未拿到 `sys.dll` / `jrt-launcher(.exe)` 落盘证据
  - 2026-05-05 阶段决策更新：接受 D3 以“仓库内 launcher PoC 已成立，最终可执行产物长时验证另列”收口，因此本批诊断也同步转为已完成收口状态
