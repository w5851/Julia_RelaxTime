生成项目依赖图（简要说明）

目的：自动化生成 `src/` 内的模块/文件依赖图，并输出为 Mermaid Markdown，放在 `docs/architecture/dependencies.md`。

快速使用：

```powershell
# 在项目根执行（Windows）
julia --project=. scripts/dev/gen_deps.jl
```

生成依赖审计报告：

```powershell
julia --project=. scripts/dev/analyze_deps.jl
```

文档一致性检查：

```powershell
julia --project=. scripts/dev/check_docs_consistency.jl
```

归档元信息补齐（历史文档一次性修复）：

```powershell
julia --project=. scripts/dev/backfill_archived_frontmatter.jl
```

Benchmark 阈值检查：

```powershell
julia --project=. scripts/dev/check_benchmark_thresholds.jl
```

散射截面热点基准（legacy vs fast_path）：

```powershell
julia --project=. scripts/dev/benchmark_total_cross_section_hotpath.jl
```

总截面固定点 baseline 导出：

```powershell
julia --project=. scripts/dev/export_total_cross_section_fixedpoint_baseline.jl
```

PNJL 扫描固定点 baseline 导出：

```powershell
julia --project=. scripts/dev/export_pnjl_scan_fixedpoint_baseline.jl
```

unit skip/deprecated 门禁检查：

```powershell
julia --project=. scripts/dev/check_unit_skip_policy.jl
```

active 文档治理检查（命名 + 归档触发）：

```powershell
julia --project=. scripts/dev/check_active_docs_governance.jl
```

PNJL 迁移门禁检查（限制 src/pnjl 新增核心实现）：

```powershell
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
```

PNJL 裁剪波次门禁（删除前快照 + 删除后 smoke + 失败回滚提示）：

```powershell
julia --project=. scripts/dev/run_prune_wave_gate.jl
```

PNJL 下线前检查表生成（`src/models/scans/*` + `src/pnjl/PNJL.jl`）：

```powershell
julia --project=. scripts/dev/generate_pnjl_decommission_checklist.jl --base HEAD --head HEAD
```

产物路径：
- `outputs/results/pnjl_decommission_checklist_<timestamp>.md`

执行台账追加（append-only，不回读历史正文）：

```powershell
julia --project=. scripts/dev/append_exec_log.jl \
	--task-file docs/dev/active/2026-02-26_多重派发重构与PNJL迁移下线开发任务单.md \
	--batch "Batch N2" \
	--goal "移除旧入口并完成冒烟校验" \
	--code-change "删除 src/pnjl 兼容路径调用点" \
	--cmd "julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD" \
	--artifact "outputs/results/pnjl_prune_wave_snapshot_20260226_101500.txt" \
	--result "通过" \
	--mainline "N2"
```

说明：
- 若未提供 `--log-file`，脚本会优先根据 `--task-file` 自动推导同目录“执行台账”；
- 若台账不存在会自动创建标准骨架后再追加；
- 仅做末尾追加（append-only），默认不读取台账历史正文。

检查表新增内容：
- `PNJL` 顶层 `scans/*` 默认 include 审计（排除白名单）
- `PNJL` scan loader 可见性清单（`load_*scan*!`）
- `PNJL.solve/PNJL.solve_multi` 外部调用方审计（排除 `src/pnjl/**`）

说明：`run_prune_wave_gate.jl` 会在执行前输出 `models_invokelatest_allowlist.toml` 的增删摘要（若当前 diff 有变更），并将摘要写入快照文件，便于统一入口审阅。

同时会额外生成独立审计产物：
- `outputs/results/models_invokelatest_allowlist_delta_<timestamp>.txt`
- `outputs/results/pnjl_scan_default_include_audit_<timestamp>.txt`

PR 范围裁剪门禁示例：

```powershell
julia --project=. scripts/dev/run_prune_wave_gate.jl --base origin/main --head HEAD
```

失败自动回滚（仅 working-tree 模式）：

```powershell
julia --project=. scripts/dev/run_prune_wave_gate.jl --auto-rollback
```

PR 范围检查示例：

```powershell
julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base origin/main --head HEAD
```

检查项包含：
- `docs/guides/**/*.md` 中历史路径（如 `test_unit/`、`julia server.jl`、`.\\start.bat`、`doc/domain-knowledge/`）
- 当 `README.md` 标注“修复中”时，guides 中是否出现“系统完全就绪/已知问题: 无”等绝对化状态词
- 当 `README.md` 标注“已验证可用”时，guides 中是否残留“修复中”表述
- `src/models/**/*.jl` 新增 `Base.invokelatest(...)` 时必须命中迁移门禁白名单（防止新增分散 world-age 边界）
- 输出 `models-invokelatest-audit` 漂移提示（`observed` vs `allowlist_baseline`）以辅助审计白名单漂移
- 白名单单一来源：`config/ci/models_invokelatest_allowlist.toml`（门禁脚本与文档均以该文件为准）
- 当白名单文件在当前 diff 中变更时，门禁输出 `models-invokelatest-allowlist-diff`（added/removed 摘要）辅助 PR 审阅
- 禁止在 `src/**` 新增 `PNJL.run_tmu_scan/run_trho_scan` 运行时调用（要求统一走 `Models.run_tmu_scan/run_trho_scan`）
- 禁止在 `src/tests/scripts` 新增 `PNJL.TmuScan/PNJL.TrhoScan` 子模块直连依赖（要求统一走 `PNJL.run_*` 顶层入口）
- 禁止在 `src/tests/scripts` 新增对 `src/pnjl/PNJL.jl` 的直接 include 路径依赖（要求统一走 `src/models/Models.jl` 入口）
- 输出 `pnjl-scan-runtime-dependency-audit`，审计 `src/**`（排除 `src/pnjl/**`）中残留 `PNJL.run_*` 调用点
- 禁止在 `src/pnjl/PNJL.jl` 新增 `scans/*` 默认 include（避免默认耦合回流；当前仅允许门禁白名单项）
- 输出 `pnjl-scan-default-include-audit`，审计 `src/pnjl/PNJL.jl` 中 `scans/*` 默认 include 现状（排除白名单）；非 0 将触发门禁失败

可复现安装（推荐）：

```powershell
npm install
```

手动渲染 SVG（可选）：

```powershell
npm run deps:render
```

输出文件：`docs/architecture/dependencies.md`（脚本覆盖）
辅助文件：`docs/architecture/dependencies.mmd`、`docs/architecture/dependencies.svg`

注意：当前脚本为 MVP：
- 只解析 `include("...")` 和 `using .Module` / `import .Module`（带点号的内部模块引用）
- 忽略第三方包依赖（通过 `Pkg.dependencies()` 可单独导出）
- 在发现循环依赖时会在输出中列出强连通分量

后续增强建议：
- 将 `web/js` 的 ES module import 解析加入图
- 在 CI 中加入自动检查并阻止新增循环依赖
- 使用 `docs/architecture/dependency_rules.md` 维护目录级依赖矩阵

---

## 归档开发文档

自动归档 `docs/dev/active` 中的文档到 `docs/dev/archived`，并添加元信息头部。

### 交互式归档

```powershell
julia --project=. scripts/dev/archive_docs.jl -i
```

### 批量归档指定文件

```powershell
julia --project=. scripts/dev/archive_docs.jl file1.md file2.md
```

### 使用自定义日期归档

```powershell
julia --project=. scripts/dev/archive_docs.jl -d 2026-01-15 file1.md
```

### 验证已归档文件格式

```powershell
julia --project=. scripts/dev/archive_docs.jl -c
```

### 预览归档操作（不实际执行）

```powershell
julia --project=. scripts/dev/archive_docs.jl --dry-run file1.md
```

**功能说明**：
- 自动添加 YAML 元信息头部（title、archived、original、archived_date）
- 自动提取文档标题（从第一个 # 标题）
- 自动添加日期前缀到文件名（如果尚未存在）
- 从 active 目录移动到 archived 目录
- 支持批量操作和交互式选择
