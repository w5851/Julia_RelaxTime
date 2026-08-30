# <工作流名称>

> 本文件是 SOP 模板。复制后必须在 `config/governance/docs_authority_map.toml` 登记；不要直接把模板登记为 active SOP。

## 1. 目的与适用范围

说明这份 SOP 负责什么计算、面向什么使用者，以及完成后得到什么。

## 2. 非适用范围

列出不由本 SOP 负责的模型、参数区间、研究结论或产物类型。

## 3. 权威入口

- 稳定 CLI：`scripts/...`
- 统一 API：`Models...`
- 默认配置：`config/...`
- 脚本白名单：链接到 `docs/guides/scripts/README.md`

不得把内部 include 文件或历史兼容路径写成当前入口。

## 4. 物理口径、单位与参数约束

至少说明：

- 外部输入单位与内部自然单位边界；
- 化学势是 `mu_q` 还是 `mu_B`；
- 无量纲量；
- 参数值域和不可静默接受的非法输入；
- 必须保持不变的物理/求解语义。

## 4A. 公式路线闭合门禁（适用时）

当本任务引入新物理领域、新可观测量、新生产候选路线，或从外部文献引入会
改变计算语义的公式时，必须先登记公式路线闭合包。闭合包至少包含：

- 微观模型起点、近似层次和符号/单位表；
- 从模型到输出的逐步推导，及每个外部公式的来源、方程号和转换因子；
- 公式 → 代码 → 测试映射、守恒量/极限/维度不变量；
- 当前状态（`draft`、`candidate` 或 `production_authorized`）与未决项；
- 适用的 `config/governance/formula_route_closure.toml` 登记和检查命令。

公式闭合通过不等于数值收敛或 production 授权；两者必须分别在第 8、11、13
节给出证据。未解决的符号、归一化、正则化或边界问题必须保留为未决项，不能
由代码默认值静默决定。

## 5. 输入配置及优先级

说明默认配置、命名 profile、环境变量和 CLI 覆盖顺序，并指明最终有效配置如何落盘。

## 6. 环境与版本冻结

记录 Julia 版本、根项目环境、git commit、sysimage 适配策略和外部数据版本。

## 7. Smoke 预检

给出可复制的低成本命令，并明确它验证什么、不验证什么。

## 8. 收敛性验证

列出必须扫描的离散化/积分/求解参数、比较观测量和进入正式计算的判据。不得通过随意放宽容差绕过失败。

## 9. 正式计算命令

给出稳定 wrapper 命令或配置驱动入口。正式参数必须来自收敛性证据，不得直接复用 smoke 参数。

## 10. 输出目录与产物合同

区分：

- `data/outputs/results/`：CSV、JSON、README、日志、审计与 manifest；
- `data/outputs/figures/`：图像和 `plot_manifest.json`。

列出必须产物、字段单位、缺失值语义和排序假设。

## 11. Regression / Validation 验收

明确对应的 unit、integration、regression、validation 层及最小命令。

## 12. 失败点、断点续算与重跑

说明失败点 sidecar、resume/overwrite、配置不兼容、部分产物和重跑策略。

## 13. Diagnostic 与 Formal Production 的边界

列出 diagnostic-only 参数或输出，以及升格 formal production 必需的收敛、回归、溯源和人工审阅证据。

## 14. 后处理与作图

说明后处理是否只消费冻结产物、图像输出目录和 plot manifest 要求。

## 15. 关联公式、API 和测试

链接权威公式、API、测试和稳定脚本指南，不复制大段事实源。

## 16. 最后验证记录

- 验证日期：`YYYY-MM-DD`
- git commit：由实际执行者填写
- 执行命令：与 authority map 一致
- 结果：通过 / 失败 / diagnostic-only
- 备注：记录环境差异或未覆盖项
