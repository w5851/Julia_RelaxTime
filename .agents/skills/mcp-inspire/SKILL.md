---
name: mcp-inspire
description: Use when Codex needs to retrieve scholarly literature through INSPIRE-first workflows, including INSPIRE search/sorting, open-access arXiv PDF download, compliant triage between open files and publisher landing pages, and supervised Playwright institution-login handoff for legally authorized full-text access.
---

# MCP Inspire

先走开放获取路径；只有在用户明确说明自己拥有合法机构权限时，才进入受监督的出版社访问流程。

## Core policy
- 优先 INSPIRE / arXiv / 出版社公开文件，不要先走付费墙。
- 只有在用户已声明拥有合法机构权限时，才使用 Playwright 自动化“打开页面、导航到机构登录、登录后下载用户本来就有权限下载的 PDF”。
- 不要绕过付费墙，不要尝试规避机构认证，不要引导使用盗版镜像站。
- 默认不要替用户保存账号密码、验证码、cookie 导出文件或长期凭据。若用户明确要求在其本人机器上自动化机构登录，可接受“仅本机可解密”的本地密文方案，例如 Windows DPAPI；不要把明文密码写入 repo、skill 文档或可共享文件。短信、邮箱验证码、MFA、人机验证仍应优先保留人工接管。

## Open-access workflow（单篇）
1. `inspire_search`：用关键词/条件检索。
2. 选定 `recid` 后用 `inspire_get_literature` 或 `inspire_list_open_files` 查看详情与可下载候选。
3. 若存在 `arxiv.pdf` 或 INSPIRE files：用 `download_open_file` 下载。
4. 若只有 DOI/出版社页：用 `resolve_publisher_link` 得到 `doi.org`（可选 `PROXY_PREFIX` 代理），先记录 DOI、期刊、年份、落地页。
5. 只有在无开放文件且用户具备合法机构权限时，才切换到下述 Playwright 流程。

## Batch workflow（批量开放下载）
- 优先用 `inspire_batch_download_arxiv_pdfs`：
  - 输入：`q`（必填），可选 `sort`（例如按被引/时间等，具体取决于 INSPIRE API 支持），以及 `limit/pages` 控制规模。
  - 输出：每条记录的 `recid/arxivId/url/savedTo` 或跳过原因。

建议的安全习惯：
- 先用小 `limit`（比如 3-5）跑通流程，再逐步加大。
- 需要大批量时，把下载目录指向一个临时文件夹，避免污染项目目录。

## Institution-access workflow（Playwright，受监督）
前提：用户已经说明自己拥有合法机构权限，例如学校/单位 VPN、图书馆代理、Shibboleth/OpenAthens、校园统一身份认证。

1. 先完成开放获取分流，确认该条记录没有可直接下载的公开文件。
2. 用 `resolve_publisher_link` 拿到 DOI 或出版社落地页。
3. 若本机存在 `D:\Desktop\paper\dev\publisher_gateway_session.py`，优先使用它而不是临时拼 Playwright CLI 命令；它已经把“西安交大图书馆网关直达常用物理出版社 -> 站内搜 DOI -> 点击 PDF -> 复用当前会话 cookie 下载真实 PDF”固化成单命令和批量命令。
4. 若本机已由用户配置 `publisher_gateway_session.py init-secret` 生成的 DPAPI 密文文件，可允许脚本自动填充学校登录页的账号密码。
5. 若仍出现验证码、邮箱验证、MFA 或人机验证，暂停自动化，由用户手动完成。
6. 用户完成认证后，脚本继续自动执行 DOI 搜索、PDF 按钮点击和保存。
7. 只下载当前账号本来就有权访问的文件，并把 PDF、DOI、来源站点、保存路径记到本地目录。

建议输出目录：
- `output/playwright/literature-access/`：截图、trace、网络日志等临时调试产物。
- `output/literature/<publisher-or-recid>/`：PDF 与题录元信息。

执行细节和命令模板见 [references/institution-access.md](references/institution-access.md)。

### XJTU gateway shortcut
当用户明确说明自己拥有西安交通大学图书馆权限，且目标主要是物理期刊（APS / IEEE / IOP / Springer / Wiley / Nature / Taylor & Francis）时，优先使用下面这条受监督路线。`ScienceDirect` 当前不属于默认稳定自动下载集合，只有在开放版本缺失且用户明确指定时才单独尝试：

1. 先检查 INSPIRE files / arXiv / 出版社开放版本；只有确认没有开放 PDF 时，才进入下面的机构脚本步骤。
2. 选择对应 `provider`。
3. 默认运行带本机 DPAPI 密文凭据参数的命令，例如：
   `python D:\Desktop\paper\dev\publisher_gateway_session.py fetch-doi --provider aps --doi 10.1103/PhysRevD.101.054019 --session-name xjtu-gateway --auto-login --secret-file D:\Desktop\paper\dev\local_secrets\xjtu_library_login.json`
4. 若仍进入验证码/MFA 页，再人工接管；完成页面认证后脚本自动继续，无需默认等待终端输入。
5. 脚本自动进入数据库官网、站内搜索 DOI、点击 `PDF`/`Download PDF`、并保存到 `D:\Desktop\paper\dev\outputs\gateway_fetch\`。

批量场景：

1. 准备一个 DOI 文本文件，每行一个 DOI，允许 `#` 注释行。
2. 默认运行带本机 DPAPI 密文凭据参数的批量命令，例如：
   `python D:\Desktop\paper\dev\publisher_gateway_session.py fetch-doi-batch --provider aps --doi-file D:\Desktop\paper\dev\outputs\gateway_fetch\batch_dois.txt --session-name xjtu-gateway --auto-login --secret-file D:\Desktop\paper\dev\local_secrets\xjtu_library_login.json`
3. 保持同一浏览器进程和持久会话，避免每篇都重新登录。

当前经验规则：
- `aps` 已验证可作为默认稳定分支。
- `ieee` 已验证可用。2026-05-06 的实测结果是：脚本可通过西安交大代理域名直接进入 IEEE 站内搜索结果页、匹配 DOI 对应文章、点击 `PDF`，并从 `stampPDF/getPDF.jsp` 保存真实 PDF 到本地。
- `sciencedirect` 不要默认走机构自动下载分支；先查 INSPIRE / arXiv。若用户明确指定 `sciencedirect`，先提醒其该站在当前环境下存在浏览器/IP 风控，且 2026-05-06 的实测结果是：手工 `Chrome` 可过，`Edge` 手工与 Playwright 会话会落到 `IP blocked` 拦截页。

## Sorting notes
- `inspire_search` 支持可选 `sort` 参数（透传到 INSPIRE API）。
- “相关性”一般是默认排序；若要“被引量/最新”，需要 INSPIRE 侧支持对应的 `sort` 值。

## Handoff to PDF utilities
- 下载完成后，可用 PDF MCP 工具读取元信息/抽取文本（通常要求绝对路径）。

## Boundary
- `download_open_file` 仅允许下载白名单开放域名（INSPIRE files、arXiv）。
- 机构访问仅限“用户已合法授权 + 用户手动完成敏感认证 + 自动化只做页面导航和下载”。
- 若页面条款、验证码、DRM 或下载限制表明不应自动下载，则停止并改为人工访问建议。
