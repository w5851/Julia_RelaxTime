---
name: formal-production-artifact
description: 在 Julia_RelaxTime 中生产可入库的正式高精度数值产物。适用于 formal production, 正式产物, 高精度产物, 收敛性测试, convergence, production audit, README/manifest, plot_manifest, data/outputs/results, data/outputs/figures。必须先完成收敛性证据，再用证据支持的参数重跑正式产物。
---

# formal-production-artifact

## Purpose

将“生产高精度正式产物”收束成可审计的闸门式流程：

- 先证明参数精度足够；
- 再用通过收敛性测试的参数重跑正式产物；
- 最后只提交正式数据、正式图像、审计文档和必要脚本改动。

本 Skill 负责产物生产治理，不负责发明新物理口径、重写稳定入口或替代数值回归 Skill。

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
- 不猜测路线：优先使用当前主项目稳定入口、脚本和 documented workflow。
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
- 预期验证层：smoke、regression、validation 或 governance checks。

## Standard Workflow

1. Inspect current contract
   - 读取仓库指令、相关 API/guides、目标脚本、入口和既有输出目录。
   - 若触及 `src/models/`、`src/relaxtime/` 或 `src/simulation/`，先读对应入口和测试层。

2. Define production case
   - 选定 `<domain>/<topic>/<case_slug>`。
   - 规划：
     - `data/outputs/results/<domain>/<topic>/<case_slug>/`
     - `data/outputs/figures/<domain>/<topic>/<case_slug>/`
   - 明确哪些现有输出是临时诊断，哪些可作为 comparison input。

3. Run convergence gate
   - 在 result 目录下创建 `convergence/`。
   - 至少比较低/中/高三档参数，或比较相邻高精度档 `N` 与 `N + Delta`。
   - 记录每档完整命令、参数、环境、git commit、输出路径、耗时和失败点。
   - 生成机器可读比较表和简短人工摘要。

4. Judge convergence
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

5. Produce formal outputs
   - 只有 verdict 为 `production-grade` 时才能重跑正式 production。
   - 使用通过收敛性测试的最高参数，或有证据支持的足够收敛参数。
   - 正式数据写入 result 目录，不覆盖 convergence 原始证据。
   - 正式图像写入 figure 目录。

6. Write audit artifacts
   - result 目录至少包含：
     - `README.md`
     - `PRODUCTION_AUDIT.md`
     - manifest（例如 `manifest.json` 或既有等价文件）
   - result-side 文档必须记录：
     - source CSV / data files；
     - figure 目录相对路径；
     - 图像生成命令；
     - 参数、环境、git commit；
     - 收敛性 verdict 与残余风险；
     - smoke/regression/governance 检查结果。
   - figure 目录必须包含 `plot_manifest.json` 或既有等价索引，记录来源数据路径、hash 或相对路径、绘图参数、格式和 dpi。

7. Validate and clean
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
