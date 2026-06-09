---
name: formal-production-artifact
description: 在 Julia_RelaxTime 中生产可入库的正式高精度数值产物。适用于 formal production, 正式产物, 高精度产物, 收敛性测试, convergence, production audit, README/manifest, plot_manifest, data/outputs/results, data/outputs/figures。必须先完成收敛性证据，再用证据支持的参数重跑正式产物。
---

# formal-production-artifact

## Purpose

将“生产高精度正式产物”收束成可审计的闸门式流程：

- 先证明参数精度足够；
- 再优先通过已注册的 GitHub Actions workflow 用通过收敛性测试的参数重跑正式产物；
- 最后只提交正式数据、正式图像、审计文档和必要脚本改动。

本 Skill 负责产物生产治理，不负责发明新物理口径、重写稳定入口或替代数值回归 Skill。若项目中已有稳定
`workflow_dispatch` 生产入口，应优先使用 GitHub Actions 触发高精度生产，以便复用固定 runner、统一依赖、artifact
留痕和手动触发入口。

## Apply When

- 用户要求生产、重跑、入库、归档或发布正式数值产物。
- 结果需要写入 `data/outputs/results/<domain>/<topic>/<case_slug>/`。
- 图像需要写入 `data/outputs/figures/<domain>/<topic>/<case_slug>/`。
- 任务要求高精度、收敛性证据、`PRODUCTION_AUDIT`、manifest 或 `plot_manifest.json`。
- 需要把临时诊断输出升级为正式产物前的 gate。

## Prefer Another Skill When

- 目标是建立或维护测试基线：先用 `baseline-regression-governance`。
- 目标是判断 relaxtime/transport 数值漂移：配合 `transport-regression-keeper`。
- 目标只是 append-only 记录实验批次：配合 `experiment-logbook-append`。
- 目标是隔离复现文献：先用 `literature-reproduction-spike`，不要直接写正式产物。

## Hard Gates

- 不读入口就不运行：先读 `AGENTS.md`、`.github/copilot-instructions.md`、相关 `docs/api/`、`docs/guides/scripts/`、目标脚本和入口。
- 不猜测路线：优先使用当前主项目稳定入口、脚本和 documented workflow；存在已注册的 GitHub Actions 生产 workflow 时，优先使用 action 而不是本机长跑。
- 不绕过默认分支注册：新增或修改 workflow 后，必须先通过 PR 合入默认分支或确认目标 ref 上的 workflow 已被 GitHub 注册，再触发正式 production action。
- 不把 action 成功等同于正式入库：workflow success 只表示 artifact 已生成；只有通过收敛性、通道诊断、失败点和审计 gate 后，才能把 artifact 晋升为 production-grade。
- 不先产正式数据：正式 production 前必须创建 result-side `convergence/` 目录并完成收敛性比较。
- 不只看脚本成功：收敛性比较必须覆盖核心观测量、关键比值、status counts、NaN/Inf、负值、异常点和失败点。
- 收敛不通过就停止：不得生成或标注为正式产物；只能报告 `blocked` 或 `diagnostic-only`。
- diagnostic 精度必须明示：若只能支撑诊断用途，README/AUDIT 必须显式标注，不能写成 production-grade。
- 不混放产物：figure-side 只放 PNG/SVG/PDF 和 `plot_manifest.json`；CSV、README、audit、manifest 和 convergence 证据留在 result-side。
- 不提交临时试跑：提交前清理或排除 scratch、preview、failed run、临时 convergence 草稿。

## Required Scope Lock

正式运行前必须写清或在汇报中锁定：

- 物理口径与目标观测量；
- 扫描路径、参数网格、单位和输出字段；
- diagnostic policy、异常值处理和缺失值语义；
- 非目标范围；
- 目标 result 目录与对应 figure 目录；
- GitHub Actions workflow 名称、触发 ref、输入参数和预期 artifact 名称；若本机运行，必须说明不用 action 的具体原因；
- 预期验证层：smoke、regression、validation 或 governance checks。

## GitHub Actions First Policy

- 优先顺序：
  1. 已在默认分支注册且支持 `workflow_dispatch` 的正式/生产 workflow；
  2. 已在当前目标 ref 注册的专用 production workflow；
  3. 本机脚本运行，仅限 convergence gate、quick diagnostic、workflow 不存在/不可用或用户明确要求本机执行。
- 触发前检查：
  1. `gh workflow list` 或 `.github/workflows/*.yml` 中存在目标 workflow；
  2. workflow 支持所需输入（case name、run tier、物理口径、节点数、是否绘图等）；
  3. 当前代码、workflow 和脚本已经在触发 ref 上可用；
  4. case slug 不覆盖旧正式产物，除非用户明确要求替换；
  5. 预期 artifact retention 足够下载和审计。
- 触发后必须记录：
  - run URL / run id / head SHA / branch；
  - workflow 输入；
  - artifact 名称、下载目录、结果 CSV、图目录；
  - action audit note 或 `PRODUCTION_AUDIT.md` 的 verdict；
  - 后续是否需要下载、比较、晋升或重跑。
- action artifact 默认视为候选产物。若 workflow 自带 audit 写明 `diagnostic-only`，最终汇报必须沿用该边界，不能口头晋升为正式产物。

## Standard Workflow

1. Inspect current contract
   - 读取仓库指令、相关 API/guides、`.github/workflows/`、目标脚本、入口和既有输出目录。
   - 若触及 `src/models/`、`src/relaxtime/` 或 `src/simulation/`，先读对应入口和测试层。

2. Define production case
   - 选定 `<domain>/<topic>/<case_slug>`。
   - 规划：
     - `data/outputs/results/<domain>/<topic>/<case_slug>/`
     - `data/outputs/figures/<domain>/<topic>/<case_slug>/`
   - 明确哪些现有输出是临时诊断，哪些可作为 comparison input。

3. Choose execution surface
   - 优先选择可手动触发的 GitHub Actions workflow，并锁定 workflow 文件、触发 ref 和输入参数。
   - 若 workflow 是本轮新增或修改的，先完成 PR/merge，让 GitHub 注册 workflow，再从目标 ref 触发。
   - 若必须本机执行，记录 fallback 原因、命令、环境和与 action 口径的差异。

4. Run convergence gate
   - 在 result 目录下创建 `convergence/`。
   - 至少比较低/中/高三档参数，或比较相邻高精度档 `N` 与 `N + Delta`。
   - 记录每档完整命令、参数、环境、git commit、输出路径、耗时和失败点。
   - 生成机器可读比较表和简短人工摘要。
   - 若 convergence 由 action 运行，下载 artifact 后把比较表和摘要放入 result-side `convergence/`。

5. Judge convergence
   - 比较对象至少包括：
     - 核心观测量；
     - 关键比值或派生量；
     - status counts；
     - NaN/Inf counts；
     - 负值或物理域外值；
     - 异常点、失败点、边界点。
   - 给出三态 verdict：
     - `production-grade`
     - `diagnostic-only`
     - `blocked`

6. Produce formal outputs
   - 只有 verdict 为 `production-grade` 时才能重跑正式 production。
   - 使用通过收敛性测试的最高参数，或有证据支持的足够收敛参数。
   - 优先通过 GitHub Actions 生产 workflow 触发正式重跑；下载并审计 artifact 后，再决定是否导入仓库数据树。
   - 正式数据写入 result 目录，不覆盖 convergence 原始证据。
   - 正式图像写入 figure 目录。

7. Write audit artifacts
   - result 目录至少包含：
     - `README.md`
     - `PRODUCTION_AUDIT.md`
     - manifest（例如 `manifest.json` 或既有等价文件）
   - result-side 文档必须记录：
     - source CSV / data files；
     - figure 目录相对路径；
     - 图像生成命令；
     - 参数、环境、git commit；
     - GitHub Actions run URL / head SHA / artifact 名称（若使用 action）；
     - 收敛性 verdict 与残余风险；
     - smoke/regression/governance 检查结果。
   - figure 目录必须包含 `plot_manifest.json` 或既有等价索引，记录来源数据路径、hash 或相对路径、绘图参数、格式和 dpi。

8. Validate and clean
   - 运行与变更风险匹配的 smoke/regression/governance checks。
   - 若数值路径、扫描逻辑、solver 或 relaxtime workflow 受影响，显式决定是否需要 regression/validation。
   - 提交前确认只包含正式数据、正式图像、审计文档和必要脚本改动。

## Minimum Audit Fields

`PRODUCTION_AUDIT.md` 至少包含：

- production case；
- physics scope；
- non-goals；
- command log；
- convergence matrix；
- selected production parameters；
- data outputs；
- figure outputs；
- validation commands and results；
- known limitations and residual risks。

## Final Report Contract

最终汇报必须区分：

- `written/configured`：已经写入或配置的文件；
- `effective/usable`：已经通过收敛性、生产重跑和验证证明可用的正式产物；
- `not run / skipped`：未运行检查及原因；
- `residual risk`：仍需人工判断或后续补证的部分。

若 verdict 不是 `production-grade`，最终汇报不得使用“正式产物已完成”。
