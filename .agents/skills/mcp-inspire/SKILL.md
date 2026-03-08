---
name: mcp-inspire
description: Workflows for using the local mcp-inspire MCP server in VS Code (INSPIRE search, sorting, arXiv PDF download, and handoff to PDF utilities).
license: MIT
metadata:
  author: local
  version: "1.0.0"
---

# MCP Inspire

## When to apply
- 需要在 VS Code / Copilot 里通过 MCP 工具检索 INSPIRE 文献。
- 需要把检索结果里的 arXiv PDF 批量下载到本地，供后续阅读/摘录。
- 需要在“可自动下载的开放资源”与“出版社落地页（可能付费墙）”之间做合规分流。

## Core workflow（单篇）
1. `inspire_search`：用关键词/条件检索。
2. 选定 `recid` 后用 `inspire_get_literature` 或 `inspire_list_open_files` 查看详情与可下载候选。
3. 若存在 `arxiv.pdf` 或 INSPIRE files：用 `download_open_file` 下载。
4. 若只有 DOI/出版社页：用 `resolve_publisher_link` 得到 `doi.org`（可选 `PROXY_PREFIX` 代理），必要时用 `open_url` 手动打开。

## Batch workflow（批量）
- 优先用 `inspire_batch_download_arxiv_pdfs`：
  - 输入：`q`（必填），可选 `sort`（例如按被引/时间等，具体取决于 INSPIRE API 支持），以及 `limit/pages` 控制规模。
  - 输出：每条记录的 `recid/arxivId/url/savedTo` 或跳过原因。

建议的安全习惯：
- 先用小 `limit`（比如 3-5）跑通流程，再逐步加大。
- 需要大批量时，把下载目录指向一个临时文件夹，避免污染项目目录。

## Sorting notes
- `inspire_search` 支持可选 `sort` 参数（透传到 INSPIRE API）。
- “相关性”一般是默认排序；若要“被引量/最新”，需要 INSPIRE 侧支持对应的 `sort` 值。

## Handoff to PDF utilities
- 下载完成后，可用 PDF MCP 工具读取元信息/抽取文本（通常要求绝对路径）。

## Boundary
- 不做出版社付费墙自动化登录/爬取。
- `download_open_file` 仅允许下载白名单开放域名（INSPIRE files、arXiv）。
