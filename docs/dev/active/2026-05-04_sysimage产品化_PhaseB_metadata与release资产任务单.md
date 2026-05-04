# sysimage 产品化 Phase B：metadata、bootstrap 与 release 资产任务单

更新日期：2026-05-04

当前状态：B1/B2/B3/B4 首轮已收口；metadata 首版、release 资产命名 / 平台矩阵、bootstrap 获取脚本、wrapper 版本策略已形成可执行契约，可作为后续 release workflow 的直接输入。

> 目的：定义 sysimage metadata 契约、预构建 sysimage release 资产的命名 / 平台矩阵、bootstrap 获取入口，以及 wrapper 的版本不匹配策略，为后续 release workflow 与默认稳定 CLI 提供统一契约。

---

## 1. 目标

- [x] B1 定义 sysimage metadata 契约
- [x] B2 设计 GitHub Release 资产命名与平台矩阵
- [x] B3 设计 bootstrap 下载脚本
- [x] B4 设计“版本不匹配时”的回退策略

---

## 2. 范围与非范围

### 2.1 本期范围

- [x] `scripts/dev/build_sysimage.jl` 产出 metadata 字段扩展
- [x] metadata 字段说明
- [x] release 资产命名规范
- [x] 首版平台矩阵
- [x] bootstrap 脚本契约与最小实现
- [x] wrapper mismatch policy 契约与最小实现

### 2.2 非范围

- [x] 不实现 GitHub Actions release 工作流
- [x] 不直接修改 wrapper 自动联网下载逻辑
- [x] 不在本期引入 app / service 入口

---

## 3. metadata 契约（Phase B 首版）

### 3.1 必填字段

- `generated_at`
- `julia_version`
- `git_commit`
- `sysimage_path`
- `precompile_script`
- `platform_family`
- `platform_os`
- `platform_arch`
- `artifact_basename`
- `release_asset_name`
- `release_archive_format`

### 3.2 字段语义

- `platform_family`
  - 取值：`windows | linux | macos`
- `platform_os`
  - 原始平台名，建议直接记录 Julia / host 识别到的 OS 语义
- `platform_arch`
  - 如 `x86_64 | aarch64`
- `artifact_basename`
  - 本地 sysimage 文件名，如 `JuliaRelaxTime.dll`
- `release_asset_name`
  - 发布到 GitHub Release 的压缩资产名
- `release_archive_format`
  - `zip | tar.gz`

### 3.3 本地 metadata 与 release 资产的关系

- metadata 文件既服务于本地 wrapper 判断，也服务于未来 release 资产下载。
- 本地 `sysimage_path` 保留绝对路径，方便当前机器直接消费。
- release 相关字段不要求本地真实存在同名压缩包，但命名必须已稳定。

---

## 4. release 资产命名规范（首版）

统一命名：

`JuliaRelaxTime-sysimage-<platform_family>-<platform_arch>-julia<julia_version>.<archive_ext>`

示例：

- `JuliaRelaxTime-sysimage-windows-x86_64-julia1.12.5.zip`
- `JuliaRelaxTime-sysimage-linux-x86_64-julia1.12.5.tar.gz`
- `JuliaRelaxTime-sysimage-macos-aarch64-julia1.12.5.tar.gz`

命名原则：

- [x] 带上平台族
- [x] 带上架构
- [x] 带上 Julia 版本
- [x] 不把 git commit 放进资产文件名
  - commit 应留在 metadata 中，而不是放进 release 文件名导致资产名频繁变化

---

## 5. 平台矩阵（首版建议）

### 5.1 第一批支持

| platform_family | platform_arch | archive | sysimage ext | 优先级 |
|---|---|---|---|---|
| `windows` | `x86_64` | `zip` | `.dll` | P0 |
| `linux` | `x86_64` | `tar.gz` | `.so` | P0 |
| `macos` | `aarch64` | `tar.gz` | `.dylib` | P1 |
| `macos` | `x86_64` | `tar.gz` | `.dylib` | P2 |

### 5.2 暂不纳入首版

- `windows/aarch64`
- `linux/aarch64`

原因：

- 当前优先保证主开发与常见用户平台
- 先降低 release 资产与 CI 复杂度

---

## 6. 任务分解

### B1 metadata 契约

- [x] 扩展 `build_sysimage.jl` metadata 字段
  - 验收：本地 `build/JuliaRelaxTime.sysimage.json` 含首版契约字段
- [x] 固定平台族与 release 资产命名规则
  - 验收：metadata 中可直接导出 release 资产名

### B2 release 资产命名 / 平台矩阵

- [x] 将资产命名规范和平台矩阵写入任务单
  - 验收：未来可直接照此落 CI / release

### B3 bootstrap 获取脚本

- [x] 新增 `scripts/dev/bootstrap_sysimage.ps1`
  - 验收：Windows / PowerShell 下可基于 `platform_family + platform_arch + julia_version` 解析 release 资产 URL
- [x] 新增 `scripts/dev/bootstrap_sysimage.sh`
  - 验收：Linux / macOS 下可基于同一命名规则解析 release 资产 URL
- [x] 固定 release 资产最小内容契约
  - 验收：压缩包内至少包含 `JuliaRelaxTime.<ext>` 与 `JuliaRelaxTime.sysimage.json`

推荐命令：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1
```

```bash
sh scripts/dev/bootstrap_sysimage.sh
```

bootstrap 规则：

- 当前平台、架构、Julia 版本三者共同决定目标资产名
- 默认走 GitHub `releases/latest/download/...`
- 如需固定 release，可显式传 `-ReleaseTag <tag>` 或 `--release-tag=<tag>`
- 当前脚本只负责“获取并解包到 `build/`”，不负责在 wrapper 中隐式联网

### B4 版本不匹配策略

- [x] Windows wrapper 支持 `-MismatchPolicy fallback|strict|rebuild`
- [x] POSIX wrapper 支持 `--mismatch-policy=fallback|strict|rebuild`
- [x] 保留旧参数 `-BuildIfMissing` / `--build-if-missing` 作为 `rebuild` 别名

策略定义：

- `fallback`
  - 默认策略；若 sysimage 缺失或不兼容，则退回普通 `julia --project=.` 运行
- `strict`
  - 要求必须存在兼容 sysimage；否则直接失败
- `rebuild`
  - 若 sysimage 缺失或不兼容，则本地调用 `scripts/dev/build_sysimage.jl` 重建

兼容性判定首版：

- `julia_version` 必须匹配当前 `julia --version`
- `platform_family` 必须匹配当前平台族
- `platform_arch` 必须匹配当前架构

示例：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 -MismatchPolicy strict scripts/pnjl/calculate_phase_structure.jl --preset=smoke
```

```bash
sh scripts/dev/run_with_sysimage.sh --mismatch-policy=rebuild scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true
```

---

## 7. 验证

- [x] 运行 `julia --project=. scripts/dev/build_sysimage.jl`
- [x] 检查 `build/JuliaRelaxTime.sysimage.json`
- [x] `powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1 -DryRun -Overwrite`
- [x] `bash scripts/dev/bootstrap_sysimage.sh --dry-run --overwrite`（本机通过 Git Bash 验证）
- [x] `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 -MismatchPolicy fallback scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`

---

## 8. DoD

- [x] metadata 契约首版落地
- [x] release 资产命名规范首版落地
- [x] 平台矩阵首版落地
- [x] bootstrap 获取入口首版落地
- [x] wrapper mismatch policy 首版落地
- [x] 可作为后续 bootstrap / release workflow 的直接输入
