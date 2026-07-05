---
status: ready_for_archive
owner: both
created: 2026-07-04
updated: 2026-07-05
target_project: D:\Desktop\Julia_RelaxTime
source_project: D:\Desktop\paper\My_Paper\pnjl_aniso_transport
production_case_slug: figure4_phase_diagram_prod_v1
completion_scope: completed_production_artifact
---

# PRD: Production-Grade PNJL Phase Reference And Figure Assets

## 0. 任务单状态

- [x] 任务单编制完成，已同步到仓库 `docs/dev/active/`。
- [x] 人工决策已固化：full-grid `xi=-0.5:0.05:0.5`、spinodal 必算、GitHub Actions 优先但非强制、本项目只交付正式产物。
- [x] 生产边界、convergence gate、result-side / figure-side 输出和验收标准已写清。
- [x] P4-A / P4-B 已完成：现有 reference 仅作 comparison input，执行面优先使用更新后的 `PNJL Dense Reference` workflow。
- [x] P4-C 至 P4-H 已完成：正式数据、正式双联图、manifest、README、PRODUCTION_AUDIT 与外部消费索引已生成。

## 1. 背景与目标

`pnjl_aniso_transport` 英文期刊稿计划加入 Figure 4，用于为 finite `mu_B` transport results 提供独立 PNJL 相结构语境。

当前 writing package 中已有相图候选和本地可复现 review candidate，但用户判断新生成的单图信息较密集，不适合直接作为正式正文图。短期策略是使用现有相图结果作为占位；长期需要由 `D:\Desktop\Julia_RelaxTime` 产出 production-grade 正式相图数据与图像资产。

本项目侧目标限定为生产正式产物，不负责论文侧 caption、claim table、正文叙事或 `main.tex` 替换。

核心目标：

- 生产覆盖 `xi=-0.5:0.05:0.5` 的正式 PNJL phase reference 数据，以便支持未来不同图像选择。
- 计算 first-order boundary、CEP、crossover、spinodal；spinodal 数据必做，图中是否展示可选。
- 生成一张可作为外部消费示例的双联 PNJL 相图资产：`T-mu_B` + `T-rho`。
- 先完成 convergence gate，再用通过 gate 的参数生成正式产物。
- 输出 result-side 审计包、figure-side 图像与 `plot_manifest.json`。

## 2. 文档落位与生产 case

本 PRD 的仓库内 canonical 任务单为：

- `D:\Desktop\Julia_RelaxTime\docs\dev\active\2026-07-04_Figure4正式PNJL相图生产任务单.md`

正式产物采用隔离 case，不覆盖旧 `data/outputs/figures/pnjl/phase_diagrams` 预览图：

- result-side：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/`
- figure-side：`data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/`

result-side 至少包含：

- `README.md`
- `PRODUCTION_AUDIT.md`
- `manifest.json`
- `convergence/`
- `figure_assets/`
- 正式 phase reference CSV 或指向正式 reference root/tag 的索引

figure-side 只放图像文件和 `plot_manifest.json`：

- `figure4_phase_diagram_TmuB_Trho.png`
- `figure4_phase_diagram_TmuB_Trho.pdf`
- `plot_manifest.json`

## 3. 范围

- [x] 生产 `xi=-0.5:0.05:0.5` 的正式 phase reference 数据。
- [x] 计算 first-order boundary、CEP、crossover、spinodal，并记录 row counts、失败点和异常点。
- [x] 生成 `T-mu_B` 面板，而不是 quark chemical potential `mu_q` 面板。
- [x] 生成 `T-rho` 面板，用于展示密度表示下的相结构或共存区。
- [x] 在图中保留 crossover 信息，但线型、透明度和图例形式可重新设计，不要求沿用 dotted style。
- [x] spinodal 数据必须生产；正式主图中是否显示 spinodal 由绘图参数控制。
- [x] 输出 PNG/PDF，优先满足正式文档消费的阅读清晰度。
- [x] 输出 manifest，记录脚本、源数据、参数、git commit、轴变换、row counts 和图形样式。

## 4. 非范围

- [ ] 不在本 PRD 中重写 PNJL phase solver 的物理判据。
- [ ] 不把 transport 曲线用于反推 CEP 或一阶相变边界。
- [ ] 不在本 PRD 中决定投稿充分性。
- [ ] 不修改 paper 项目的 `main.tex`、caption、claim table 或 writing package manifest。
- [ ] 不负责外部论文叙事审核；本项目只交付正式数据、图像资产、manifest 和 audit。
- [ ] 不把 GitHub Actions 成功或脚本成功直接等同于 production-grade。

## 5. 现状盘点

已存在候选资产：

- `D:\Desktop\Julia_RelaxTime\data\outputs\figures\pnjl\phase_diagrams\phase_diagram_combined.png`
- `D:\Desktop\Julia_RelaxTime\data\outputs\figures\pnjl\phase_diagrams\phase_diagram_T_mu.png`
- `D:\Desktop\Julia_RelaxTime\data\outputs\figures\pnjl\phase_diagrams\phase_diagram_T_rho.png`
- `D:\Desktop\Julia_RelaxTime\data\reference\pnjl\boundary.csv`
- `D:\Desktop\Julia_RelaxTime\data\reference\pnjl\cep.csv`
- `D:\Desktop\Julia_RelaxTime\data\reference\pnjl\spinodals.csv`
- `D:\Desktop\Julia_RelaxTime\data\reference\pnjl\crossover_dense.csv`
- `D:\Desktop\Julia_RelaxTime\data\reference\pnjl\paper_p1_mott_phase_isentropic_20260531\`

writing package 中已生成 review candidate：

- `My_Paper/pnjl_aniso_transport/writing_package/attachments/phase_diagrams/generated/figure4_phase_diagram_T_muB_candidate.png`
- `My_Paper/pnjl_aniso_transport/writing_package/attachments/phase_diagrams/generated/figure4_phase_diagram_T_muB_candidate_manifest.json`

当前判断：

- 旧 `phase_diagrams/*.png` 只能作为视觉参考或占位，不作为正式生产 source。
- writing-package review candidate 可追溯，但不是 `Julia_RelaxTime` 生产管线产物。
- 现有 `data/reference/pnjl` 资产可作为 comparison input；是否作为正式 source 必须通过 schema、manifest、参数与 convergence 审计。
- `paper_p1_mott_phase_isentropic_20260531` 具有较完整 manifest，但其 `xi` 为 `-0.3,0.0,0.3`，不满足本次正式产物要求的完整 `xi=-0.5:0.05:0.5` 覆盖。

已关闭的缺口：

- [x] 旧相图横轴为 `mu` 或 `mu_q`，正式图资产已使用 `mu_B`。
- [x] 当前旧 PNG 与当前 `plot_phase_diagram.py` + reference CSV 的直接复现存在 schema/source gap；本轮已新建正式后处理入口。
- [x] `plot_phase_diagram.py` 仍偏向旧字段与单图预览用途，已保留为预览入口，不作为 Figure 4 正式脚本。
- [x] 图形密度已通过默认 `xi` subset 与 spinodal 可选图层控制。
- [x] 已完成 convergence gate，并将 reference 和图像标注为 production-grade。

## 6. 数据源与 schema 要求

- [x] 正式图必须基于正式 phase reference root/tag，而不是直接从旧 PNG 或 paper 侧 candidate 派生。
- [x] reference CSV 必须覆盖 `xi=-0.5:0.05:0.5` 的 first-order boundary、CEP、crossover、spinodal。
- [x] boundary 和 spinodal CSV 必须包含稳定排序字段，或在 figure asset 生成阶段显式生成：
  - `curve_parameter`
  - `plot_order_key`
- [x] CEP schema 必须直接包含 `muB_CEP_MeV`，或 manifest 中明确记录 `muB_CEP_MeV = 3 * muq_CEP_MeV` 的生成过程。
- [x] crossover schema 必须明确 `mu_MeV` 是 `mu_q`，并在 figure asset 中生成 `muB_MeV = 3 * mu_MeV`。
- [x] `T-rho` 面板使用的 `rho` 字段必须记录单位为 `rho/rho_0`，并说明 missing-value 语义。
- [x] figure asset 中间表应统一输出到 result-side `figure_assets/`，避免只保留 PNG/PDF。
- [x] figure asset 应支持从完整 `xi` 数据中选择绘图 subset，不要求主图展示全部 `xi`。
- [x] spinodal 图层必须能通过绘图参数开启或关闭；关闭时仍应在 result-side 保存 spinodal 数据与审计统计。

## 7. 执行面选择

优先审计并复用已注册 `workflow_dispatch`，以节省本机算力：

- `.github/workflows/pnjl-phase-diagram.yml`
- `.github/workflows/pnjl-dense-reference.yml`

触发正式 production 前必须确认：

- [ ] workflow 是否支持 `xi=-0.5:0.05:0.5`、`T/rho` 网格、crossover、boundary、spinodal 和 CEP。
- [ ] workflow 是否支持或可显式记录 `p_num`、`t_num`、`iterations` 等精度参数。
- [ ] workflow artifact 是否包含可下载的 CSV、manifest、validation report 和图像。
- [ ] 目标 workflow 已在触发 ref 上注册；若本轮需要修改 workflow，应先合入默认分支或确认 GitHub 已注册后再触发。

GitHub Actions 优先但不是硬性要求。若 workflow 不支持正式 case 的必要参数或审计字段，则允许本机运行，但必须在 `PRODUCTION_AUDIT.md` 中记录 fallback 原因、完整命令、环境、git commit，以及与 action 口径的差异。稳定 CLI 本机运行优先使用：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl ...
```

或在需要 per-`xi` 聚合 reference 时使用：

```powershell
julia --project=. scripts/pnjl/build_dense_phase_reference.jl ...
```

## 8. Convergence gate

正式 production 前必须创建：

- `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/convergence/`

建议三档：

- C0 diagnostic：低成本参数，仅用于发现明显 schema/流程问题。
- C1 full-grid candidate：覆盖 `xi=-0.5:0.05:0.5`，建议 `p_num=24,t_num=8,iterations=80`，`T_step=5 MeV`，`rho_step=0.05`，作为候选正式口径。
- C2 refined check：提高积分节点或加密关键窗口网格，可在代表性 `xi` 锚点与 CEP/低温边界邻域执行；建议锚点至少覆盖 `xi=-0.5,-0.25,0.0,0.25,0.5`。若 C2 只做局部检查，必须在 audit 中标注 coverage 和 residual risk。

比较对象至少包括：

- boundary row counts 和按 `(xi,T_MeV)` 匹配后的 `muB`/`rho_hadron`/`rho_quark` 差异。
- CEP 是否存在，以及 `T_CEP_MeV`、`muB_CEP_MeV` 漂移。
- crossover row counts、converged counts、`T_crossover_MeV` 漂移、非单调或异常点。
- spinodal row counts、低温异常点、分支排序稳定性和物理域外值。
- NaN/Inf counts、负值或物理域外值、失败点、异常跳变点。
- `curve_parameter`/`plot_order_key` 的稳定排序是否有效。

verdict 只能是：

- `production-grade`
- `diagnostic-only`
- `blocked`

只有 `production-grade` 才能进入正式图生成与外部消费。

## 9. 功能需求

### 9.1 数据内容

- [x] 正式数据覆盖 `xi=-0.5:0.05:0.5`。
- [x] 正式数据包含 first-order boundary、CEP、crossover、spinodal。
- [x] 每类数据均记录 row counts、`xi` coverage、单位、失败点和异常点。
- [x] spinodal 数据即使不进入默认图，也必须随正式产物发布和审计。

### 9.2 图像内容

- [x] 默认图像资产采用双联布局：左图 `T-mu_B`，右图 `T-rho`。
- [x] 左图至少包含 first-order boundary、CEP、crossover。
- [x] 右图至少包含 coexistence region 或 phase-boundary density representation。
- [x] 图像脚本支持 `xi` subset 参数；默认展示 subset 可从完整数据中选择，不要求展示全部 `xi=-0.5:0.05:0.5`。
- [x] spinodal 图层支持参数化开关，默认是否开启由图形密度审核决定。

### 9.3 视觉要求

- [x] 图中信息不能比当前 writing-package 单图 candidate 更密集。
- [x] crossover 必须保留，但可使用更轻的线型、颜色、透明度或独立图例策略。
- [x] first-order boundary、CEP、crossover 的图例必须清晰区分。
- [x] 若保留 spinodal，应使用低权重视觉编码，避免与 crossover/first-order boundary 混淆。
- [x] 坐标轴必须明确标注 `mu_B (MeV)`、`T (MeV)`、`rho/rho_0`。

### 9.4 可追溯输出

- [x] 输出正式图：`figure4_phase_diagram_TmuB_Trho.{png,pdf}`。
- [x] 输出 `plot_manifest.json`：记录 source CSV、row counts、`xi` values、axis conversion、plot styles、生成脚本、git commit、dpi、figure size。
- [x] 输出用于绘图的中间 CSV 或 JSON 到 result-side `figure_assets/`。
- [x] manifest 明确 `mu_B=3*mu_q` 的转换来源和适用字段。
- [x] `PRODUCTION_AUDIT.md` 区分 written/configured、effective/usable、not run/skipped、residual risk。

## 10. 实施任务

- [x] P4-0：将 PRD 同步为 repo active doc，并补齐正式 production gate。
- [x] P4-Decision：固化生产覆盖、spinodal、execution surface 和本项目交付边界。
- [x] P4-TaskSheet：任务单编制完成，可进入生产执行阶段。
- [x] P4-A：审计现有 phase reference 数据源，决定复用、规范化还是重跑。
- [x] P4-B：审计 `.github/workflows/` 与本地脚本入口，决定 action-first 或本机 fallback。
- [x] P4-C：运行 convergence gate，生成机器可读比较表和人工摘要。
- [x] P4-D：用通过 convergence 的参数生成正式 full-grid phase reference 或锁定正式 reference root/tag。
- [x] P4-E：实现或更新正式绘图脚本，支持 `T-mu_B + T-rho` 双联输出、`xi` subset、spinodal 开关、figure assets 和 `plot_manifest.json`。
- [x] P4-F：生成正式 PNG/PDF/manifest，并写入 result-side README / PRODUCTION_AUDIT。
- [x] P4-G：运行 validation/governance checks，并人工目视审核图像可读性、图例、坐标轴、line-style 语义。
- [x] P4-H：输出正式产物索引与外部消费说明，不修改 paper 项目文件。

## 11. 验收标准

### 11.1 written/configured

- [x] PRD 已落位到 `docs/dev/active/`。
- [x] production `xi` coverage 已锁定为 `-0.5:0.05:0.5`。
- [x] spinodal 已锁定为必算数据、可选图层。
- [x] production case 目录、source 选择、执行命令、脚本入口、输出文件名已写清。
- [x] 绘图脚本或后处理入口已能生成双联 `T-mu_B + T-rho` 图。
- [x] result-side 与 figure-side manifest schema 已写出。

### 11.2 effective/usable

- [x] convergence verdict 为 `production-grade`。
- [x] 正式数据覆盖 `xi=-0.5:0.05:0.5`，并包含 first-order boundary、CEP、crossover、spinodal。
- [x] `T-mu_B` 面板横轴为 `mu_B (MeV)`，不是 quark `mu_q`。
- [x] `T-rho` 面板与 `T-mu_B` 面板使用一致的 plotted `xi` subset，并且 subset 可追溯到完整正式数据。
- [x] CEP markers 来源可追溯到正式 `cep` artifact 或正式 phase summary。
- [x] Crossover line 来源可追溯到正式 crossover artifact。
- [x] Spinodal 数据来源可追溯；若图中关闭 spinodal，manifest 明确记录关闭原因或图形参数。
- [x] 图像能在目标阅读宽度下阅读，图例不遮挡核心曲线。
- [x] `plot_manifest.json` 能回答：图从哪些 CSV 生成、用了哪些 plotted `xi`、完整 source `xi` 覆盖是什么、是否做了 `mu_B=3mu_q` 转换、生成脚本和 commit 是什么。
- [x] `PRODUCTION_AUDIT.md` 记录 convergence matrix、selected production parameters、validation commands、known limitations。

## 12. 风险与回退

- 风险：完整 `xi=-0.5:0.05:0.5` 数据生产成本较高。
  - 回退：优先使用 GitHub Actions；若 action 不足以支撑参数或 artifact 审计，本机 fallback 并记录原因。
- 风险：crossover、spinodal、first-order boundary 同图显示过密。
  - 回退：正文式主图只保留 first-order boundary、CEP、crossover；spinodal 保留在数据和可选图层中。
- 风险：完整 `xi` 趋势不适合单张主图展示。
  - 回退：正式数据保留全覆盖；图像只展示 subset，另可生成 supplement / appendix 风格版本。
- 风险：旧 reference 与新 production pipeline 结果不一致。
  - 回退：以 production manifest 和最新 phase reference 为准，旧 reference 只作为 comparison input。
- 风险：convergence 不通过。
  - 回退：不得标注为正式图；只能保留为 `diagnostic-only` 或 `blocked`。
- 风险：GitHub Actions 不支持必要 production 参数。
  - 回退：先补 workflow 并等注册，或使用本机 fallback 并在 audit 中说明。

## 13. 外部消费边界

本项目完成后只返回正式产物索引与消费说明：

- result-side path
- figure-side path
- `manifest.json`
- `PRODUCTION_AUDIT.md`
- `plot_manifest.json`
- 推荐消费的正式 CSV / figure asset 文件列表

本项目不负责把产物写入 paper 项目，也不负责论文 caption、claim table、正文段落或 LaTeX 修改。

## 14. 当前用户决策记录

- 2026-07-04：新生成的 `T-mu_B` 单图候选没有明显错误，但信息可能过密，不适合直接放入正文。
- 2026-07-04：短期使用现有结果占位；正式结果后续通过 PRD 交给 `Julia_RelaxTime` 产出。
- 2026-07-04：最终目标图采用双联 `T-mu_B + T-rho`。
- 2026-07-04：crossover 是重要相结构信息，需要保留，但不强制使用 dotted 线型。
- 2026-07-04：正式产物必须先过 convergence gate，再进入正式图和外部消费。
- 2026-07-04：正式产物数据覆盖锁定为 `xi=-0.5:0.05:0.5`。
- 2026-07-04：spinodal 数据必须计算，以防未来需要；图中是否展示作为可选参数。
- 2026-07-04：优先走 GitHub Actions 以节省本机算力，但不是必须。
- 2026-07-04：本项目侧只负责产出正式产物，不处理论文侧细节。

## 15. 待人工判断点

当前没有阻塞正式产物完成的人工决策。仍建议用户后续人工判断：

- 是否接受默认主图 subset：`xi=-0.5,-0.25,0.0,0.25,0.5`。
- spinodal 是否需要另出 supplemental / appendix 版本；默认 PNG/PDF 已关闭 spinodal。
- paper 项目侧是否采用本默认双联图，或基于同一正式数据再派生更适合正文版面的裁剪/子集版本。

## 16. P4-A/P4-B 审计记录

日期：2026-07-04。

### P4-A 数据源审计结论

- 现有 `data/reference/pnjl/boundary.csv` 与 `cep.csv` 覆盖 `xi=-0.5:0.1:0.5`，不满足本任务锁定的 `xi=-0.5:0.05:0.5` full-grid 要求。
- 现有 `data/reference/pnjl/spinodals.csv` 只覆盖 `xi=0.0,0.2,0.4`，不满足 spinodal 必算且 full-grid 覆盖的要求。
- 现有 `data/reference/pnjl/crossover_dense.csv` 覆盖 `xi=-0.5:0.05:0.5`，但 manifest 标记为 `crossover_only=true`，只能作为 crossover comparison input，不能单独作为正式 full reference source。
- `data/reference/pnjl/paper_p1_mott_phase_isentropic_20260531/` 有较完整 schema、`curve_parameter` / `plot_order_key` 和 `p_num=24,t_num=8` 记录，但只覆盖 `xi=-0.3,0.0,0.3`，只能作为 schema/精度口径参考与局部 comparison input。
- P4-A 决策：不复用任何现有 reference 作为正式 source；正式 source 需要用 `scripts/pnjl/build_dense_phase_reference.jl` 重跑 full-grid reference，并把现有资产仅作为 comparison input。

### P4-B 执行面审计结论

- `PNJL Phase Diagram` workflow 已注册且可手动触发，但当前设计偏预览图：不暴露 `p_num/t_num/iterations`，输出/上传路径偏向旧 `data/reference/pnjl/*.csv` 与 `phase_diagrams/*.png`，不适合作为本任务正式 production 入口。
- `PNJL Dense Reference` workflow 已注册且可手动触发，基础方向与本任务更接近，但原 workflow 固定 `--crossover-only`，无法产出 boundary/CEP/spinodal。
- 已更新 `.github/workflows/pnjl-dense-reference.yml`：支持 full-reference 模式、`p_num/t_num/iterations` 输入、`xi_min/xi_max/xi_step` full-grid 默认、以及完整 CSV artifact 上传。
- 已更新 `scripts/pnjl/validate_dense_reference_artifact.py`：支持 full-reference artifact 校验，并兼容带 `curve_parameter` / `plot_order_key` 的 crossover schema。
- P4-B 决策：优先走更新后的 `PNJL Dense Reference` workflow；但由于 workflow 修改需要先在 GitHub 注册，正式 action 触发前必须确认目标 ref 已包含该 workflow 版本。若未注册或 GitHub Actions 不可用，则使用本机 fallback 运行 `scripts/pnjl/build_dense_phase_reference.jl`，并在 `PRODUCTION_AUDIT.md` 记录原因。
- `gh workflow view "PNJL Dense Reference" --yaml` 显示当前 GitHub 已注册版本仍是旧版：固定 `--crossover-only`，且没有 `p_num/t_num/iterations` 与 `crossover_only` 输入。因此在本轮 workflow 修改合入/注册前，不应触发 action 作为 P4-C production/convergence 执行面。

### P4-A/P4-B 自动验证记录

- `python -m py_compile scripts/pnjl/validate_dense_reference_artifact.py`：通过。
- `python -c "import yaml, pathlib; yaml.safe_load(pathlib.Path('.github/workflows/pnjl-dense-reference.yml').read_text(encoding='utf-8')); print('yaml ok')"`：通过。
- `python scripts/pnjl/validate_dense_reference_artifact.py --reference-root data/reference/pnjl --tag dense --min-crossover-rows 1 --expect-crossover-only`：通过；确认旧 `dense` artifact 为 crossover-only comparison input。
- `python scripts/pnjl/validate_dense_reference_artifact.py --reference-root data/reference/pnjl/paper_p1_mott_phase_isentropic_20260531 --tag paper_p1_mott_phase_isentropic_20260531 --min-crossover-rows 1 --expect-full-reference`：通过；确认该旧 full-reference schema 可被新 validator 接受，但 `xi` 覆盖不足。
- `julia --project=. scripts/dev/check_docs_consistency.jl`：通过。
- `julia --project=. scripts/dev/check_script_entrypoints.jl`：通过。
- `gh workflow list`：确认 `PNJL Dense Reference` 与 `PNJL Phase Diagram` 均已注册；但 `gh workflow view "PNJL Dense Reference" --yaml` 确认远端注册版本尚未包含本轮修改。
- `git diff --check`：通过。
- `julia --project=. scripts/dev/check_active_docs_governance.jl`：未通过，原因是两个既有 stale active doc（`2026-04-30_介子数密度与BU工作流任务单.md`、`2026-05-01_PNJL可选功能盘点与优先级任务单.md`）超过 60 天；与本任务新增/修改内容无直接关系。
- 本机未能运行 `bash -n` 类 shell 语法检查，因为当前 PowerShell 环境未提供 `bash`；workflow 的 `run` 段已改为显式 `shell: bash` 并移除 validator array 中不必要的反斜杠续行。

### 推荐 P4-C 起点

- P4-C 应先跑 convergence gate，不直接生成正式产物。
- 首选 action 输入：`tag=figure4_phase_diagram_prod_v1`、`crossover_only=false`、`xi_values=""`、`xi_min=-0.5`、`xi_max=0.5`、`xi_step=0.05`、`p_num=24`、`t_num=8`、`iterations=80`、`T_min=60` 或按 convergence 档显式覆盖。
- 若需要覆盖低温 first-order / spinodal 低温段，应在 convergence 档中显式比较 `T_min=30` 与 `T_min=60` 的影响。

## 17. P4-C convergence gate 执行记录

日期：2026-07-05。

### C1 full-grid candidate

- GitHub Actions run：`https://github.com/w5851/Julia_RelaxTime/actions/runs/28736257287`
- head SHA：`aca739f6bac4565dddef888abf9588b8c0fa583f`
- workflow：`PNJL Dense Reference`
- 输入：`tag=figure4_phase_diagram_prod_v1_c1_p24t8`、`crossover_only=false`、`xi_values=""`、`xi_min=-0.5`、`xi_max=0.5`、`xi_step=0.05`、`T_min=60`、`T_max=240`、`T_step=5`、`rho_min=0.0`、`rho_max=4.0`、`rho_step=0.05`、`p_num=24`、`t_num=8`、`iterations=80`、`crossover_n_mu=16`、`crossover_mu_max=450`。
- action 结果：成功；`Run dense reference builder`、`Validate dense reference artifact`、`Upload dense reference artifact` 均成功。
- artifact 已下载到 result-side convergence staging：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/convergence/action_artifacts/c1_p24t8/`。
- C1 行数摘要：boundary 304、CEP 21、spinodal 304、crossover 336；crossover converged 336/336。
- C1 覆盖：`xi=-0.5:0.05:0.5` 共 21 个点；四类 CSV 均无 NaN/Inf，除 `xi` 外未发现负值。
- 下载后直接运行旧 validator 会出现 manifest path mismatch，因为 manifest 记录的是 action runner staging path，下载位置不同；该现象属于 artifact relocation，不代表 action 内部校验失败。

### C2 refined anchors

- GitHub Actions run：`https://github.com/w5851/Julia_RelaxTime/actions/runs/28741709519`
- head SHA：`aca739f6bac4565dddef888abf9588b8c0fa583f`
- workflow：`PNJL Dense Reference`
- 输入：`tag=figure4_phase_diagram_prod_v1_c2_p32t12_anchors`、`crossover_only=false`、`xi_values="-0.5,-0.25,0.0,0.25,0.5"`、`T_min=60`、`T_max=240`、`T_step=5`、`rho_min=0.0`、`rho_max=4.0`、`rho_step=0.05`、`p_num=32`、`t_num=12`、`iterations=100`、`crossover_n_mu=16`、`crossover_mu_max=450`。
- action 结果：成功；`Run dense reference builder`、`Validate dense reference artifact`、`Upload dense reference artifact` 均成功。
- artifact 已下载到 result-side convergence staging：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/convergence/action_artifacts/c2_p32t12_anchors/`。
- C2 行数摘要：boundary 72、CEP 5、spinodal 72、crossover 80；crossover converged 80/80。
- C2 覆盖：`xi=-0.5,-0.25,0.0,0.25,0.5`；四类 CSV 均无 NaN/Inf，除 `xi` 外未发现负值。

### C1 vs C2 比较与 verdict

- 机器可读比较产物：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/convergence/c1_vs_c2_anchor_comparison/`。
- 人工摘要：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/convergence/P4-C_CONVERGENCE_SUMMARY.md`。
- 核心差异：boundary `mu_transition_MeV` 最大差约 `0.2614 MeV`；CEP `T_CEP_MeV` 最大差约 `0.0293 MeV`，`muB_CEP_MeV` 最大差约 `0.0760 MeV`；crossover `T_crossover_MeV` 最大差约 `0.0304 MeV`；spinodal quark 分支 `mu` 最大差约 `1.4341 MeV`。
- 较大相对差主要出现在 crossover `rho` / `derivative` 诊断量和部分 spinodal density；这些记录为 residual risk，但不改变核心 `T-mu_B` 相线和 CEP 收敛判断。
- P4-C verdict：`production-grade`。C1 full-grid artifact 可作为 P4-D 正式 full-grid phase reference source，需在后续 result-side audit 中保留上述 residual risk。

## 18. P4-D 正式 phase reference source lock

日期：2026-07-05。

- 正式 reference root：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/reference/`。
- 正式 reference tag：`figure4_phase_diagram_prod_v1_c1_p24t8`。
- source manifest：`data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/phase_reference_source_manifest.json`。
- source kind：GitHub Actions artifact promoted copy。
- source run：`https://github.com/w5851/Julia_RelaxTime/actions/runs/28736257287`。
- source head SHA：`aca739f6bac4565dddef888abf9588b8c0fa583f`。
- 数值参数：`p_num=24`、`t_num=8`、`iterations=80`、`T_min=60`、`T_max=240`、`T_step=5`、`rho_min=0.0`、`rho_max=4.0`、`rho_step=0.05`。
- 覆盖：`xi=-0.5:0.05:0.5`，共 21 个点。
- 文件：boundary、CEP、crossover、crossover meta、phase reference manifest、spinodals CSV 已复制到正式 reference root，并在 source manifest 中记录 SHA256。
- P4-D 决策：后续 P4-E/P4-F 的图像与 audit 应以该 result-side reference root/tag 为正式输入；`convergence/action_artifacts/` 保留为原始证据，不作为默认消费入口。

## 19. P4-E/P4-F/P4-G/P4-H 正式图与审计记录

日期：2026-07-05。

### P4-E 可复用绘图脚本审计与实现结论

- 已审计 `scripts/pnjl/plot_phase_diagram.py`：可复用其 T-mu / T-rho 双联布局、按 `xi` 分组和图层开关思路，但该脚本仍以旧 `data/reference/pnjl` 和旧字段口径为主，CEP loader 期待 `mu_CEP_MeV`，crossover loader 期待旧 `T_crossover_chiral_MeV` / `rho_chiral`，且 T-mu 面板使用 quark `mu` 而非 `mu_B`，不适合作为正式 Figure 4 入口。
- 已审计 `scripts/relaxtime/build_paper_p1_figure_assets.py`：其 result-side `figure_assets/`、figure-side `plot_manifest.json`、`phase_reference_root/tag` 与 `mu_B=3*mu_q` overlay 逻辑适合作为正式资产组织模板。
- P4-E 决策：不 patch 旧预览脚本；新增正式后处理入口 `scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py`。
- 新脚本默认读取 result-side formal reference root/tag，输出 `T-mu_B + T-rho` 双联 PNG/PDF，支持 `--xi-values` 和 `--include-spinodal`，并写出 result-side figure assets 与 figure-side `plot_manifest.json`。

### P4-F 正式图与 manifest 输出

- result-side figure assets：
  - `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure_assets/figure4_phase_lines_TmuB.csv`
  - `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure_assets/figure4_phase_lines_Trho.csv`
  - `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure_assets/figure4_phase_plot_inputs_summary.json`
- figure-side outputs：
  - `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure4_phase_diagram_TmuB_Trho.png`
  - `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure4_phase_diagram_TmuB_Trho.pdf`
  - `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/plot_manifest.json`
- result-side audit/index：
  - `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/README.md`
  - `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/PRODUCTION_AUDIT.md`
  - `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/manifest.json`
- 默认绘图参数：`xi=-0.5,-0.25,0.0,0.25,0.5`，`include_spinodal=false`，`formats=png,pdf`，`dpi=300`，`figsize=7.6,3.8`。
- `plot_manifest.json` 记录 source CSV hashes、axis conversion、row counts、plotted/source `xi`、style、script、git commit、dpi、figure size 和 residual risks。

### P4-G 验证与目视审核

- `python -m py_compile scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py scripts/analysis/pnjl/compare_phase_reference_convergence.py scripts/pnjl/validate_dense_reference_artifact.py`：通过。
- `python scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py`：通过。
- `--include-spinodal` 临时低 DPI smoke：通过；临时目录已清理。
- `plot_manifest.json` hash 自检：通过，10 个 input/asset/figure 文件均匹配。
- PNG 尺寸检查：`2244x1156`，`RGBA`。
- 目视审核：默认 no-spinodal 双联图可读；坐标轴、图例、panel labels、line-style 语义无明显遮挡。
- `julia --project=. scripts/dev/check_docs_consistency.jl`：通过。
- `julia --project=. scripts/dev/check_script_entrypoints.jl`：通过。
- `julia --project=. scripts/dev/check_data_output_path_guard.jl`：通过。
- `git diff --check`：通过。
- `julia --project=. scripts/dev/check_active_docs_governance.jl`：未通过，原因仍是两个既有 stale active doc（`2026-04-30_介子数密度与BU工作流任务单.md`、`2026-05-01_PNJL可选功能盘点与优先级任务单.md`），与本任务产物无直接关系。

### P4-H 外部消费说明

- 本项目正式消费入口为 result-side `README.md`、`manifest.json`、`PRODUCTION_AUDIT.md` 和 figure-side `plot_manifest.json`。
- 本项目不修改 paper 项目文件；paper 侧如需采用默认图或派生裁剪版本，应从上述正式产物索引读取。
- 仍需人工判断的是论文版面选择：默认 `xi` subset 是否合适、是否另出 spinodal variant、是否需要为正文宽度再派生更稀疏版本。

## 20. P4-I T-rho CEP 视觉连接调整

日期：2026-07-05。

- 用户确认默认绘图 subset `xi=-0.5,-0.25,0.0,0.25,0.5` 合适。
- 用户要求删除桌面 PRD 镜像；`D:\Desktop\PRD_Julia_RelaxTime_figure4_phase_diagram.md` 已删除，后续只保留本项目内 canonical active doc。
- 用户指出 `T-rho` 面板中 CEP 两侧曲线因附近数据间隔而视觉上未连接；本轮按图像资产层处理，不改正式 reference CSV。
- `scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py` 已加入 `T-rho` plotting-only `CEP_visual_connector` 行：每个 `xi` 增加两条 coexistence connector 和一条 crossover connector，使共存区两支与 crossover 曲线视觉上穿过 CEP。
- connector 的 `rho` 坐标沿用 CEP 在 `T-rho` 面板中的绘图坐标：同 `xi` 下最接近 CEP 温度的 boundary 两侧密度均值；该口径已写入 result-side README、PRODUCTION_AUDIT、figure asset summary 与 figure-side `plot_manifest.json`。
- 已重新生成 `figure4_phase_diagram_TmuB_Trho.png`、`figure4_phase_diagram_TmuB_Trho.pdf`、`figure_assets/` 和 `plot_manifest.json`。
- connector 检查：`figure4_phase_lines_Trho.csv` 中有 63 条 `CEP_visual_connector` 行；默认五个 plotted `xi` 对应 15 条 connector 被绘制。
