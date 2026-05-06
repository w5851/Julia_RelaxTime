# Institution Access via Playwright

只在用户已明确拥有合法机构权限时使用本流程。

## 目标
- 保持 INSPIRE-first：先证明没有开放文件，再访问出版社站点。
- 用 Playwright 做“导航、快照、点击、下载”。
- 默认仅把验证码、MFA、人机验证留给用户手动输入；若用户已在本机配置 DPAPI 密文凭据，则账号密码可由脚本自动填充。

## XJTU gateway script first
若本机存在 `D:\Desktop\paper\dev\publisher_gateway_session.py`，优先使用该脚本而不是手工拼 Playwright CLI 步骤。该脚本已经内置西安交大图书馆数据库网关直达链接，并支持：

- `init-secret`：将账号密码写入仅当前 Windows 用户可解密的 DPAPI 密文文件
- `fetch-doi`：单篇 DOI 自动获取
- `fetch-doi-batch`：批量 DOI 自动获取
- `serve` + `cmd ...`：调试或半自动控制

当前已内置的 `provider`：

- `aps`
- `ieee`
- `iop`
- `springer`
- `sciencedirect`
- `wiley`
- `nature`
- `taylor-francis`

## Current status snapshot
- `aps`：已打通。当前推荐路线是 `publisher_gateway_session.py fetch-doi --provider aps ...`，脚本可自动登录学校认证页、进入 APS、点击 PDF，并把最终 PDF 保存到 `D:\Desktop\paper\dev\outputs\gateway_fetch\`。
- `ieee`：已打通。当前推荐路线是 `publisher_gateway_session.py fetch-doi --provider ieee ...`，脚本会通过西安交大代理域名直接进入 IEEE 站内 DOI 搜索结果页，匹配文章后点击 `PDF`，并从 `stampPDF/getPDF.jsp` 保存真实 PDF 到 `D:\Desktop\paper\dev\outputs\gateway_fetch\`。
- `sciencedirect`：当前只完成到“登录后站内标题搜索并定位到 `View PDF`/`pdfft` 链接”。在 2026-05-06 的实测中，手工 `Chrome` 可通过，而 `Edge` 手工与 Playwright 会话都会落到 Elsevier 的 `IP blocked` / 拦截页，因此当前不要把 `sciencedirect` 视为稳定可用的自动下载分支。
- 对 `sciencedirect`，默认策略应改为：先查 INSPIRE / arXiv / 出版社开放版本；只有在开放版本缺失且用户明确要求时，才把它当成人工补充访问任务，而不是默认批量自动化目标。

## Session and artifact setup
- 优先使用单独 session，避免不同出版社 cookie 串扰。
- 在仓库内把 Playwright 调试产物放到 `output/playwright/literature-access/`。
- 把最终 PDF 和元信息放到 `output/literature/<label>/`，其中 `<label>` 优先使用 `recid`、DOI slug 或出版社缩写。

脚本默认产物目录：

- 会话 profile：`D:\Desktop\paper\dev\outputs\playwright_sessions\<session-name>`
- 状态快照：`D:\Desktop\paper\dev\outputs\gateway_session\<session-name>`
- PDF 输出：`D:\Desktop\paper\dev\outputs\gateway_fetch\`

## Recommended commands

初始化本机密文凭据：

```powershell
python D:\Desktop\paper\dev\publisher_gateway_session.py init-secret `
  --username 你的学号 `
  --secret-file D:\Desktop\paper\dev\local_secrets\xjtu_library_login.json
```

可选：也可额外传 `--password ...` 非交互写入，但不推荐把该命令保存在 shell 历史里。

单篇 DOI：

```powershell
python D:\Desktop\paper\dev\publisher_gateway_session.py fetch-doi `
  --provider aps `
  --doi 10.1103/PhysRevD.101.054019 `
  --session-name xjtu-gateway `
  --auto-login `
  --secret-file D:\Desktop\paper\dev\local_secrets\xjtu_library_login.json
```

这就是默认推荐命令；只有在用户明确要求不用本机密文凭据时，才退回不带 `--auto-login --secret-file ...` 的变体。

批量 DOI：

```powershell
python D:\Desktop\paper\dev\publisher_gateway_session.py fetch-doi-batch `
  --provider aps `
  --doi-file D:\Desktop\paper\dev\outputs\gateway_fetch\batch_dois.txt `
  --session-name xjtu-gateway `
  --auto-login `
  --secret-file D:\Desktop\paper\dev\local_secrets\xjtu_library_login.json
```

`doi-file` 约定：

- 每行一个 DOI
- 允许空行
- 允许 `#` 开头注释

示例：

```text
# APS sample
10.1103/PhysRevD.101.054019
10.1103/PhysRevLett.123.456789
```

## Minimal supervised loop
1. `resolve_publisher_link` 获取 DOI 或落地页。
2. 先检查 INSPIRE files / arXiv / 出版社开放版本；只有确认没有开放 PDF 时，才继续后面的机构脚本步骤。
3. 优先把 DOI 映射到合适的 `provider`，再运行默认带 `--auto-login --secret-file D:\Desktop\paper\dev\local_secrets\xjtu_library_login.json` 的 `fetch-doi` 或 `fetch-doi-batch`。
4. 脚本自动打开西安交大图书馆数据库网关直达链接。
5. 如配置了 DPAPI 密文凭据，脚本自动填充学校登录页的账号密码并提交。
6. 如后续仍进入验证码、MFA、人机验证页面，暂停并让用户手动完成；完成页面认证后脚本自动继续，无需默认等待终端回车。
7. 脚本自动进行站内 DOI 或标题搜索。
8. 脚本自动点击 `PDF` / `Download PDF` / `View PDF` 一类入口。
9. 脚本从当前会话页面提取真实 PDF 资源地址，并复用当前 cookie 直接保存 PDF。
10. 若目标文件名已存在，脚本自动追加 `-2`、`-3` 等后缀，避免批量任务被同名文件阻断。

补充说明：
- 上述第 8 步目前对 `aps` 已实测成功。
- 上述第 7-9 步目前对 `ieee` 也已实测成功；IEEE 的 PDF 页面不是裸 `.pdf`，而是 `stamp/stamp.jsp` 中嵌套 `stampPDF/getPDF.jsp`，当前脚本已兼容该结构。
- 对 `sciencedirect` 不要默认假定第 7-8 步可用；若 `View PDF` 后进入 `pdfft?...` 但页面标题仍是 `ScienceDirect`、且只出现帮助/远程访问/支持链接，应判定为出版社侧拦截，停止自动下载并回退到开放获取分流或人工处理。

## Command template
```bash
export CODEX_HOME="${CODEX_HOME:-$HOME/.codex}"
export PWCLI="$CODEX_HOME/skills/playwright/scripts/playwright_cli.sh"

"$PWCLI" --session springer-<label> open "https://doi.org/..." --headed
"$PWCLI" --session springer-<label> snapshot
# click institution-login ref from latest snapshot
# pause for user manual SSO / MFA
"$PWCLI" --session springer-<label> snapshot
# click PDF ref from latest snapshot
"$PWCLI" --session springer-<label> screenshot
```

只有在 `publisher_gateway_session.py` 不存在、或当前站点不在其内置 provider 范围内时，才退回这套原始 Playwright CLI 手工流程。

## User handoff points
以下步骤默认应交给用户本人：
- 短信验证码、邮件验证码、令牌或 MFA
- Captcha / 人机验证
- 任何带有“仅限本人使用”“同意条款”“机构使用政策”的确认页

若用户明确要求并理解风险，可把“输入账号名/密码”改为本机 DPAPI 密文自动填充；但不要把明文密码写入 repo、skill 文件或可同步目录。

## Publisher-page heuristics
- 常见入口文案：`Sign in via your institution`、`Access through your institution`、`Shibboleth`、`OpenAthens`、`Institutional login`
- 常见成功信号：页面显示 `Access provided by ...`、`PDF` 下载按钮可见、全文标签从锁定变为可读
- 若 DOI 重定向到聚合页，优先停在出版社正式落地页，再进行机构登录

## Stop conditions
- 用户未声明合法授权
- 页面明确禁止自动下载或存在明显 DRM 限制
- 认证流程要求不可转交的敏感操作且无法在用户监督下完成
- 目标站点只允许馆藏浏览、不提供 PDF 下载

## Recommended result record
至少记录：
- 标题
- DOI / recid
- 出版社域名
- 访问方式（open / institution）
- 下载时间
- 保存路径
- 是否需要后续人工补充元数据

## Batch guidance for zero-knowledge agents
- 先按 DOI / 期刊名把记录分组到具体 `provider`，不要把不同出版社混在同一个批处理命令里。
- 每个 `provider` 复用同一个 `session-name`，这样一次人工登录后可连续抓取多篇。
- 批量开始前先用 1 篇样例 DOI 跑通，确认该出版社页面结构仍兼容。
- 若某篇失败，记录 DOI 与报错并继续后续条目；不要因为单条失败中断整批。
