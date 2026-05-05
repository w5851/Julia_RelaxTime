---
title: sysimage 产品化 Phase D：D2 `create_app(...)` 适配性评估任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D2_create_app适配性评估任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase D：D2 `create_app(...)` 适配性评估任务单

更新日期：2026-05-05

当前状态：已完成首轮评估；已给出“适合 / 不适合 / 约束 / PoC 方向”结论，可直接进入 D3 最小 launcher app PoC。

> 目的：评估 `PackageCompiler.create_app(...)` 是否适合作为 launcher 层的承载形态，并明确当前仓库是“直接 app 化”还是“先做薄 package 壳再 app 化”。

---

## 1. 评估输入

- [x] D1 已完成 launcher 职责边界与首批白名单设计
- [x] 仓库当前已有 `PackageCompiler` 依赖与 `create_sysimage(...)` 使用经验
- [x] 当前仓库仍是 include-driven 结构，而不是标准 `src/PackageName.jl` package 形态

本轮依据：

- [x] [Project.toml](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/Project.toml)
  - 已声明 `PackageCompiler = "2"`
- [x] [scripts/dev/build_sysimage.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/dev/build_sysimage.jl)
  - 当前已有稳定 `create_sysimage(...)` 产线
- [x] 本地 PackageCompiler 源码
  - `create_app(package_dir, app_dir)` 要求 `package_dir` 是 package
  - 默认入口要求提供 `julia_main()::Cint`

---

## 2. 关键事实

### 2.1 `create_app(...)` 的本地 API 约束

- [x] 本地 `PackageCompiler.create_app` 只有 package-dir 入口
- [x] 它期望：
  - `package_dir` 含 package 元数据（`name` / `uuid`）
  - package 内提供 `julia_main()::Cint`
  - executable 默认映射到 package 名与 `julia_main`

证据来源：

- [x] 本地 PackageCompiler 源码 `src/PackageCompiler.jl`
  - `create_app(package_dir::String, compiled_app::String; ...)`
  - 文档明确要求“package includes a function `julia_main()::Cint`”

### 2.2 当前仓库形态

- [x] 当前仓库根目录虽有 `Project.toml`，但主实现不是标准 package 布局
- [x] 仓库当前以：
  - include-driven 源码组织
  - 稳定 CLI 白名单
  - `run_with_sysimage.*` + `bootstrap_sysimage.*`
  为主
- [x] 当前稳定入口主要是脚本，而不是 package 导出的 `julia_main`

### 2.3 当前已经具备的有利条件

- [x] 已有 sysimage 构建产线
- [x] 已有稳定 CLI 白名单
- [x] 已有 launcher 设计
- [x] 已有 Phase C 对稳定 API / 脚本边界的收敛成果

---

## 3. 评估结论

### 3.1 直接对仓库根目录运行 `create_app(...)`

- [x] 结论：`不适合`

原因：

- [x] 当前仓库根目录不是为 `create_app(...)` 准备的标准 app package
- [x] 缺少单一的 package 级 `julia_main()::Cint`
- [x] 当前稳定入口是多脚本白名单，不是“一个 package = 一个 app main”
- [x] 若强行直接 app 化，会把 D1 刚确定的 launcher / wrapper / bootstrap / stable CLI 分层重新打乱

### 3.2 以“薄 launcher package 壳”承载 `create_app(...)`

- [x] 结论：`适合`

建议形态：

- [x] 新增一个极薄的 launcher package（例如放在 `app/launcher/` 或等价目录）
- [x] 该 package 负责：
  - 提供 `julia_main()::Cint`
  - 解析 launcher 子命令
  - 选择并调度当前稳定 CLI 白名单
  - 复用当前 wrapper / bootstrap 语义或把其关键逻辑内嵌到 launcher package 中
- [x] 该 package 不负责：
  - 重写物理业务逻辑
  - 替代 `src/models` / 稳定脚本主体

### 3.3 对当前主线的总体判断

- [x] `create_app(...)` 不是当前仓库“零改造直接落地”的工具
- [x] 但它很适合作为 D3 的承载目标，只要先引入薄 launcher package 壳
- [x] 因此 D2 的最终结论是：
  - `仓库根目录直接 app 化：不推荐`
  - `launcher package 壳 app 化：推荐 PoC`

---

## 4. 主要约束

### 4.1 结构约束

- [x] 需要一个标准 package 目录
- [x] 需要显式 `julia_main()::Cint`
- [x] 需要把 launcher 命令面绑定到 package 内部入口

### 4.2 入口约束

- [x] app 不能继续依赖“用户自己记住底层脚本路径”
- [x] app 应只暴露 launcher 白名单能力
- [x] 未进入白名单的研究/分析型脚本不应混入 app 首批命令面

### 4.3 sysimage / bootstrap 约束

- [x] 现有 wrapper 的 mismatch policy 是 launcher 语义的一部分
- [x] app 化后需要明确：
  - app 是否内嵌自己的 sysimage
  - app 是否仍依赖外部 `build/JuliaRelaxTime.*`
  - app 是否替代 bootstrap，还是 bootstrap 仍独立存在

判断：

- [x] D3 首轮 PoC 更适合做“自带 app bundle，不再依赖仓库根目录 build/ 下 sysimage”的路线
- [x] bootstrap 在 app 路线下应逐步退化为“仓库开发环境方案”，而不是 app 用户默认动作

---

## 5. D3 推荐路线

### 5.1 最小 PoC 目标

- [x] 只做一个 launcher package 壳
- [x] 只暴露 1-2 个最稳定能力
  - `phase`
  - `transport-scan` 或 `unified-scan`
- [x] 验证：
  - `create_app(...)` 能成功产出 app bundle
  - app 能正确解析子命令并调起底层稳定能力
  - app 不需要用户直接接触底层脚本路径

### 5.2 首轮不做的事

- [x] 不把全部稳定白名单一次性都塞进 app
- [x] 不把 service/server 能力作为首个 app PoC 核心
- [x] 不在 D3 首轮同时解决所有 release / installer / auto-update 问题

### 5.3 推荐首个 app 命令

- [x] 优先 `phase`

原因：

- [x] README / Quickstart 已把它作为最小复现实验主链
- [x] CLI 已足够薄
- [x] 产物与成功判据明确
- [x] 对 launcher / app 演示价值最高

---

## 6. 风险

- [x] 风险 R1：为了迁就 `create_app(...)`，反向把仓库主实现改造成 package 结构，导致当前稳定主线被打乱
- [x] 风险 R2：launcher package 壳若过厚，会重新复制 wrapper / stable CLI 逻辑
- [x] 风险 R3：D3 若一次纳入过多子命令，PoC 会退化成“又一层复杂脚本集合”

应对：

- [x] 仅新增薄壳 package，不重构主仓库主体
- [x] 业务参数继续透传，避免重复实现
- [x] D3 首轮只做最小白名单命令集

---

## 7. 最终结论

- [x] `PackageCompiler.create_app(...)` 对 Julia_RelaxTime `当前根仓库形态`：`不直接适合`
- [x] `PackageCompiler.create_app(...)` 对 `薄 launcher package 壳`：`适合`
- [x] 下一步应进入 D3：
  - 先做 launcher package 壳
  - 再以 `phase` 为首个命令验证最小 app PoC

---

## 8. DoD

- [x] 已说明为什么仓库根目录不能直接作为 `create_app(...)` 输入
- [x] 已说明为什么 launcher package 壳是更合适承载形态
- [x] 已给出 D3 最小 PoC 的目标与范围
- [x] 已可直接进入 D3 实施任务单

---

## 9. 治理校验

- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [x] 两项校验均已通过，可将 D2 结论同步回 backlog
