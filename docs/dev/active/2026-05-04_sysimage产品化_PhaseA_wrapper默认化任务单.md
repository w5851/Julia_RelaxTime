# sysimage 产品化 Phase A：wrapper 默认化任务单

更新日期：2026-05-04

当前状态：首轮收口已完成；已补跨平台 wrapper 与稳定 CLI 默认入口映射，待后续按需要补跨平台实机验证。

> 目的：把稳定 CLI 的默认启动方式统一到 sysimage wrapper，并补齐 Windows 与 Linux/macOS 的最小入口约定，确保用户和零知识 agent 在同机环境下能稳定吃到冷启动优化。

---

## 1. 目标

- [x] A1 盘点所有稳定 CLI，定义哪些必须默认经 wrapper 启动
- [x] A2 为 Windows 用户入口补齐统一命令示例
- [x] A3 评估并补 `run_with_sysimage.sh`

---

## 2. 范围与非范围

### 2.1 本期范围

- [ ] `scripts/dev/run_with_sysimage.ps1`
- [ ] `scripts/dev/run_with_sysimage.sh`
- [ ] `README.md`
- [ ] `docs/guides/QUICKSTART.md`
- [ ] `docs/guides/scripts/README.md`

### 2.2 非范围

- [ ] 不进入 sysimage release / 下载分发
- [ ] 不改动稳定 CLI 主体求解实现
- [ ] 不进入 package / app / service 后续阶段

---

## 3. 任务分解

### A1 稳定 CLI 白名单映射

- [x] 梳理当前稳定 CLI 白名单中哪些入口默认应经 wrapper 启动
  - 验收：形成按平台区分的入口映射表

### A2 Windows 入口统一

- [x] 在 README / guides 中把 Windows 稳定入口主示例统一到 `run_with_sysimage.ps1`
  - 验收：稳定 CLI 主示例不再默认写成裸 `julia --project=.`

### A3 POSIX wrapper

- [x] 评估 Linux/macOS 是否需要同类 wrapper
  - 验收：给出“补 / 不补”的明确结论
- [x] 若补，新增 `scripts/dev/run_with_sysimage.sh`
  - 验收：语义与 PowerShell wrapper 对齐，至少支持 `strict` 与 `build if missing`

---

## 4. 验证

- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [ ] 若补 `.sh` wrapper，至少做一次脚本语法检查
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [ ] 若补 `.sh` wrapper，至少做一次脚本语法检查
  - 说明：当前 Windows 主机无 `bash` 可执行文件，尚未完成本机 POSIX shell 语法跑测；脚本已按 POSIX `sh` 语法编写，待 Linux/macOS 或带 bash 环境复核。

---

## 5. DoD

- [x] 稳定 CLI 默认入口已经明确走 wrapper
- [x] Windows 文档示例统一
- [x] Linux/macOS 路径有明确结论与落地
- [x] 该阶段可独立归档，不阻塞后续 Phase B
