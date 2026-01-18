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
