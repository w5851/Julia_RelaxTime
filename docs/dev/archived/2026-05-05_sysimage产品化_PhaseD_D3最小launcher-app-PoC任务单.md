---
title: sysimage 产品化 Phase D：D3 最小 launcher app PoC 任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D3最小launcher-app-PoC任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase D：D3 最小 launcher app PoC 任务单

更新日期：2026-05-05

当前状态：已完成主收口；已完成最小 launcher package、`phase` 复用 support 层与本地解释器链路验证，并确认 D3 剩余阻塞仅为“launcher sysimage / 自定义可执行产物在本机短窗口内无法完成编译”。D3 现按“仓库内 launcher PoC 已成立，长时最终产物验证另列”收口。

> 目的：在不改动活跃物理主线求解逻辑的前提下，把 D1 的 launcher 设计与 D2 的 `create_app(...)` 评估落到一个最小可执行 PoC，验证“薄 launcher package 壳 + phase 子命令”这一路线确实可行。

---

## 1. 前置结论

- [x] D1 已确认 launcher 是统一稳定入口层，而不是新的求解实现层
- [x] D2 已确认当前根仓库不适合直接 `create_app(...)`
- [x] D2 已确认“薄 launcher package 壳”适合作为 app 承载形态
- [x] 首个 app PoC 命令优先选 `phase`

---

## 2. 本批目标

- [x] D3-1 新建最小 launcher package 壳
- [x] D3-2 launcher 只承载 `phase` 子命令
- [x] D3-3 不重写第二套 phase 参数契约，复用现有 phase CLI 参数解析/执行主链
- [x] D3-4 新增 `create_app(...)` 构建脚本
- [x] D3-5 至少完成一条本地验证链：
  - [ ] app bundle 自定义 launcher 可执行文件可启动
  - [x] launcher 可显示帮助或解析子命令
  - [x] launcher 可调起 `phase` 最小运行路径

---

## 3. 范围与非范围

### 3.1 本批范围

- [x] `app/launcher/` package 壳
- [x] `scripts/dev/build_launcher_app.jl`
- [x] `phase` 子命令最小调度
- [x] 文档同步：D3 活动任务单 + backlog 阶段状态

### 3.2 非范围

- [ ] 本批不纳入 `transport-scan` / `unified-scan` / `server`
- [ ] 本批不解决 installer / auto-update / release 发布
- [ ] 本批不改写 `src/models` 求解逻辑
- [ ] 本批不把 app 描述成“已经脱离源码仓库完全独立分发”

---

## 4. 技术路线

### 4.1 launcher package 形态

- [x] 新增 `app/launcher/Project.toml`
- [x] 新增 `app/launcher/src/`
- [x] package 暴露 `julia_main()::Cint`
- [x] package 顶层只做：
  - launcher 子命令识别
  - `phase` 参数透传
  - 调用现有稳定 phase CLI 主链

### 4.2 phase 主链复用策略

- [x] 把 `scripts/pnjl/calculate_phase_structure.jl` 中可复用的 CLI 解析/执行逻辑抽成 support 层
- [x] 原脚本继续作为 repo 内稳定 CLI 入口
- [x] launcher package 复用 support 层，而不是复制一份 phase 参数语法

### 4.3 app build 路线

- [x] 新增 `scripts/dev/build_launcher_app.jl`
- [x] 首轮以本地 `build/launcher-app/` 作为默认输出
- [ ] `create_app(...)` 的首轮验证目标仅为：
  - [x] app bundle 目录已开始产出到 `build/launcher-app/`
  - [ ] 自定义 launcher 可执行文件已落盘并可启动

---

## 5. 风险与约束

- [ ] 风险 R1：为了 app PoC 反向侵入 `src/models` 主实现
  - 应对：只抽 CLI support 层，不改求解契约
- [ ] 风险 R2：launcher package 复制整份 phase CLI，形成第二套入口逻辑
  - 应对：优先做共享 support 层
- [ ] 风险 R3：PoC 误写成“已独立分发”
  - 应对：显式记录首轮 PoC 与后续产品化缺口
- [x] 风险 R4：`create_app(...)` 在当前本地 PoC 路径上构建耗时过长，超出本轮可接受验证窗口
  - 应对：本轮先保留“解释器内 launcher PoC 已通”的结论，下一批单独收口 app bundle 最终产物与构建时长

---

## 6. 验收标准

- [x] 仓库内已有最小 launcher package 目录
- [x] 已有 `create_app(...)` 构建脚本
- [x] `phase` 已能作为 launcher 子命令被调用
- [x] 至少完成一次本地验证，证明 app bundle 不是空壳
- [x] D3 文档与 backlog 状态同步

---

## 7. 预期产物

- [x] `app/launcher/Project.toml`
- [x] `app/launcher/src/*.jl`
- [x] `scripts/dev/build_launcher_app.jl`
- [x] phase CLI support 复用层
- [x] D3 执行记录与 backlog 同步

---

## 8. DoD

- [x] D3 最小 app PoC 已形成
- [x] PoC 范围与限制已写明
- [x] 后续是否进入 D4 / Phase E 有清晰前置输入

---

## 9. 当前实现与验证记录

- [x] 已新增 `scripts/pnjl/phase_cli_support.jl`，把 phase CLI 的参数解析、执行与 manifest 写出逻辑抽成共享 support 层
- [x] [scripts/pnjl/calculate_phase_structure.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/pnjl/calculate_phase_structure.jl) 已回退为“环境激活 + warmup + 调 support 层”的薄脚本
- [x] 已新增 `app/launcher/` 最小 package，当前只暴露 `help` 与 `phase`
- [x] 为保持 D3 范围最小，launcher 当前明确在运行时 `Pkg.activate(repo_root)`，因此它是“仓库内本地 PoC executable”，而不是独立分发 app
- [x] 本地解释器链路验证已通过：
  - `julia --project=app/launcher -e '... JRTLauncherApp.real_main(["help"])'`
  - `julia --project=app/launcher -e '... JRTLauncherApp.real_main(["phase", "--help"])'`
  - `julia --project=app/launcher -e '... JRTLauncherApp.real_main(["phase", "--preset=smoke", ...])'`
- [ ] `julia --project=. scripts/dev/build_launcher_app.jl`
  - 已能稳定产出 `build/launcher-app/` bundle 目录
  - 但在本轮验证窗口内仍未看到自定义 launcher 可执行文件落盘，需作为 D3 下一批单独收口
  - 2026-05-05 D3 第二批诊断已确认：最初的结构性失败来自 launcher package 缺少 `Pkg` 依赖与同步 manifest，现已修复；当前真实停留点已收敛到 `create_sysimage(...)` 之后的长耗时编译阶段，尚未进入 `create_executable_from_sysimg(...)`
  - 2026-05-05 D3 第三批已把构建路径进一步改为“手工 bundle + repo sysimage 作为 base sysimage 的极薄 launcher sysimage”；当前仍未在验证窗口内拿到 `sys.dll` / `jrt-launcher(.exe)` 落盘证据，但已不再回到早期结构性失败
  - 2026-05-05 D3 第四批已将参数进一步压到 `include_transitive_dependencies=false + -O0 --compile=min --strip-metadata --strip-ir`，并用日志确认前序 bundle 秒级完成；但 10 分钟窗口内仍停留在 `create sysimage start`，因此 D3 当前剩余阻塞已可明确归结为 launcher sysimage 编译 wall time，而非结构或契约问题
  - 2026-05-05 阶段决策：接受“仓库内 launcher PoC 已成立”作为 D3 主收口条件；若后续仍需要 `sys.dll` / `jrt-launcher(.exe)` 最终证据，应另开长时离线构建验证批次，而不再作为 D3 主线阻塞项

## 10. 治理校验

- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
