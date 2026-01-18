# 安装与环境复现（多语言）

本项目包含 Julia / Python / JavaScript 三类环境。建议根据实际使用场景选择安装范围。

## 1) Julia 环境（核心）

要求：已安装 Julia（版本以 `Project.toml` / `Manifest.toml` 为准）。

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

验证（任选其一）：

```powershell
julia --project=. scripts/server/server_full.jl
```

或：

```powershell
julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --help
```

## 2) Python 环境（脚本/分析用，可选）

Python 依赖来自 `requirements.txt`，适用于 `scripts/` 与数据处理流程。

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

> 如果你已有 Python 环境管理器（conda/uv/poetry），可自行替换上述步骤。

## 3) JavaScript 工具链（文档/可视化工具，可选）

用于 Mermaid 依赖图 SVG 渲染（`@mermaid-js/mermaid-cli`）。已添加到 `package.json` 作为开发依赖。

```powershell
npm install
```

生成依赖图与 SVG：

```powershell
julia --project=. scripts/dev/gen_deps.jl
```

若 `mmdc` 在安装时因浏览器下载失败，可用本机浏览器替代：

```powershell
$env:PUPPETEER_SKIP_DOWNLOAD="true"
$env:PUPPETEER_EXECUTABLE_PATH="C:\Program Files\Google\Chrome\Application\chrome.exe"

npm install
```

## 4) 依赖图渲染说明（补充）

- Mermaid 源文件：`docs/architecture/dependencies.mmd`
- SVG 输出：`docs/architecture/dependencies.svg`
- 文档入口：`docs/architecture/dependencies.md`

> 生成脚本会优先使用全局 `mmdc`，否则使用本地 `node_modules/.bin/mmdc`。
