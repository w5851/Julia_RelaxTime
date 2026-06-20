---
name: analysis-artifact
description: 为 Julia_RelaxTime 的正式数值产物、诊断产物或研究结果生成可追溯分析包。适用于 docs/analysis 分析产物沉淀、输入证据审计、派生表和关键图生成、机制/归因边界判定、claim ledger、论文写作前的 evidence package；可覆盖 transport、relaxtime、phase structure、meson 或其他仓库内数值产物。
---

# analysis-artifact

## Purpose

把“看结果后写解释”的过程收束为可审计分析产物：输入证据、派生表、关键图、机制边界、claim ledger 和论文写作素材都能追溯到仓库文件。

本 skill 负责分析正式或诊断产物，不负责重跑 production、不晋升数值产物、不发明外部文献解释。

## Hard Rules

- 只使用仓库内证据，除非用户明确要求新增文献检索或外部对比。
- 不修改正式 result/figure 产物；分析输出默认写入 `docs/analysis/...`。
- 先验证输入，再写结论：CSV 行数、关键字段、manifest/hash、failed points、convergence gate 和 figure manifest 必须先检查。
- 区分四层表述：
  - observation：结果表/图直接支持的趋势和局部结构。
  - attribution：可分解诊断表支持的贡献来源。
  - mechanism：传播子分母、阈值、上游解分支、算法误差等因果解释；必须有定点机制诊断或等价证据。
  - paper claim：可写入论文的表述；必须绑定 evidence file、字段、点位或图。
- 证据不足时标为 `author_check`、`inconclusive` 或 `candidate`，不要用流畅文字掩盖证据缺口。
- 不把具体 case 的数值结论写进 skill；case 信息放在分析脚本参数、README 和派生表中。

## Standard Workflow

1. Scope lock
   - 明确 artifact case、物理/数值口径、观测量、输出目录和不做事项。
   - 默认生成中文分析主文档；英文 caption 或论文正文另起任务。

2. Evidence inventory
   - 读取 result-side `README.md`、`PRODUCTION_AUDIT.md`、`manifest.json`、`effective_config.json` 或等价审计文件。
   - 读取主结果 CSV、诊断 CSV、failed/error points。
   - 读取 convergence gate 摘要和 CSV。
   - 读取 figure-side `plot_manifest.json` 和关键 PNG；必要时抽查图像。

3. Input validation
   - 检查所有输入存在、row count 符合 audit、关键字段存在。
   - 检查 NaN/Inf、失败点、负贡献、重复 key。
   - 记录输入 hash、生成命令、repo HEAD 和脚本 hash。

4. Derived evidence
   - 生成全局趋势摘要：端点差、极值、turning count、单调性分类。
   - 生成收敛风险摘要：高相对差异点、阈值分档和 gate 边界。
   - 生成 attribution 摘要：按 case 合适的维度汇总主导贡献与局部变化。
   - 若要写 mechanism，生成机制候选窗口和定点诊断表；没有机制表时不要升级解释。
   - 生成 claim ledger：每条论文候选陈述绑定证据文件、字段、点位或图号。

5. Figures
   - 生成全局 overview 图、局部结构放大图和 attribution/mechanism 图。
   - 写 `figures/plot_manifest.json`，记录输入 hash、图 hash、生成命令和脚本信息。

6. Analysis document
   - 写中文主文档，至少包含：数据与口径、全局观察、收敛风险、局部结构、机制判定、论文写作候选、作者确认项。
   - 对机制解释使用逐窗口 verdict：`supported`、`supported_indirect`、`not_dominant`、`inconclusive`。
   - 对旧产物只写 supersession 和不足边界，不把旧图作为当前 production-grade 证据。

7. Validation
   - 运行生成脚本 smoke，确认派生表无 NaN/Inf、key 唯一、claim 引用文件和字段存在、图文件非空。
   - 运行仓库治理检查：`git diff --check`、`julia --project=. scripts/dev/check_docs_consistency.jl`；新增脚本时加跑 script governance checks。
   - 汇报哪些结论可直接用于论文，哪些仍需作者确认或补充诊断。

## Output Contract

分析包目录建议为：

```text
docs/analysis/<domain>/<case>_analysis/
├── README.md
├── manifest.json
├── figures/
│   ├── *_overview*.png
│   ├── *_local*.png
│   ├── *_attribution*.png
│   └── plot_manifest.json
└── tables/
    ├── global_trend_summary.csv
    ├── convergence_risk_summary.csv
    ├── attribution_summary.csv
    ├── mechanism_*summary.csv
    └── claim_ledger.csv
```

若已有 case 专用生成脚本，优先复用脚本并更新参数；否则在 `scripts/analysis/` 的合适子目录新增只读分析脚本。
