# sysimage 产品化 Phase B：release workflow 与资产打包任务单

更新日期：2026-05-04

当前状态：B5/B6 首轮已收口；本地打包脚本与 GitHub Actions release workflow 已落地，metadata / 命名契约已接到可发布资产链路。

> 目的：在不改变当前稳定 CLI 契约的前提下，补齐 sysimage release 资产打包脚本与 GitHub Actions workflow，让 metadata 命名规则真正产出可发布资产。

---

## 1. 目标

- [x] B5 实现本地 release 资产打包脚本
- [x] B6 实现 GitHub Actions release workflow

---

## 2. 范围与非范围

### 2.1 本期范围

- [x] `scripts/dev/package_sysimage_release.jl`
- [x] `.github/workflows/sysimage-release.yml`
- [x] release 资产打包输出目录约定
- [x] release / workflow_dispatch 触发策略

### 2.2 非范围

- [x] 不修改当前物理主线求解逻辑
- [x] 不引入 app / service 层能力
- [x] 不在本期处理 release note 模板细节

---

## 3. 设计要点

### 3.1 资产内容契约

- [x] 压缩包内至少包含：
  - `JuliaRelaxTime.<ext>`
  - `JuliaRelaxTime.sysimage.json`
- [x] 压缩包外补充同名 `.sha256`

### 3.2 打包输出目录

- [x] 本地默认输出到 `build/releases/`
- [x] 资产名直接复用 metadata 中的 `release_asset_name`

### 3.3 workflow 触发

- [x] `workflow_dispatch`
  - packaging-only 或 publish-to-release 两种模式
- [x] `release.published`
  - 对既有 release 自动补传平台资产

### 3.4 首版平台矩阵

- [x] `windows-latest` / `x86_64`
- [x] `ubuntu-latest` / `x86_64`
- [x] `macos-latest` / `aarch64`

---

## 4. 任务分解

### B5 资产打包脚本

- [x] 读取 `build/JuliaRelaxTime.sysimage.json`
- [x] 校验本地 sysimage 与 metadata 存在
- [x] 生成 release archive
- [x] 生成 `.sha256`

### B6 release workflow

- [x] matrix build sysimage
- [x] matrix package asset
- [x] upload workflow artifacts
- [x] 条件发布到 GitHub Release

---

## 5. 验证

- [x] `julia --project=. scripts/dev/package_sysimage_release.jl --overwrite`
- [x] 检查 `build/releases/`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`

---

## 6. DoD

- [x] metadata 命名契约已接到可执行打包脚本
- [x] GitHub Actions 可按平台矩阵产出 sysimage 资产
- [x] release workflow 可上传对应 release assets

---

## 7. 实现说明

- `scripts/dev/package_sysimage_release.jl`
  - 读取 `build/JuliaRelaxTime.sysimage.json`
  - 按 metadata 中的 `release_asset_name` / `release_archive_format` 打包到 `build/releases/`
  - 额外生成同名 `.sha256`
- `.github/workflows/sysimage-release.yml`
  - `workflow_dispatch` 支持仅打包，或指定 `release_tag` 后发布到 GitHub Release
  - `release.published` 会对既有 release 自动补传 sysimage 资产
  - 首版矩阵固定为 `windows-latest`、`ubuntu-latest`、`macos-latest`
