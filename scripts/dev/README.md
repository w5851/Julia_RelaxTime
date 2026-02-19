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

检查项包含：
- `docs/guides/**/*.md` 中历史路径（如 `test_unit/`、`julia server.jl`、`.\\start.bat`、`doc/domain-knowledge/`）
- 当 `README.md` 标注“修复中”时，guides 中是否出现“系统完全就绪/已知问题: 无”等绝对化状态词
- 当 `README.md` 标注“已验证可用”时，guides 中是否残留“修复中”表述

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
