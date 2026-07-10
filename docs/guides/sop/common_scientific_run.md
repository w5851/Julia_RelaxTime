# 公共科学计算生命周期 SOP

## 1. 目的与适用范围

本 SOP 规定 Julia_RelaxTime 中正式科学计算共有的执行阶段：范围冻结、环境记录、smoke、收敛性、正式运行、产物审计、回归/验证和结果升格。

适用于相结构、输运、介子热力学、介子数密度及后续新增的正式数值工作流。专题参数与观测量由对应工作流 SOP 补充。

## 2. 非适用范围

- 不定义具体物理公式或替代公式文档。
- 不决定某个专题的积分节点、网格步长或容差。
- 不把探索性脚本自动变成稳定入口。
- 不授权改变 solver、mixed-meson 或 non-fixedmu joint-solve 语义。

## 3. 权威入口

- SOP 索引：`docs/guides/sop/README.md`
- 稳定脚本白名单：`docs/guides/scripts/README.md`
- 测试分层：`docs/dev/testing_governance.md`
- baseline 治理：`docs/guides/BASELINE_VERSION_MANAGEMENT.md`
- 输出路径治理：`scripts/dev/check_data_output_path_guard.jl`

具体计算必须使用专题 SOP 登记的稳定 CLI 或 `Models` 统一入口。

## 4. 物理口径、单位与参数约束

执行前必须写清：

- 温度、化学势、质量等外部输入是 MeV 还是 `fm^-1`；
- 使用 `mu_q`、`mu_B` 还是 flavor-resolved chemical potentials；
- `xi`、Polyakov loop 等无量纲量的定义；
- 模型配置来源和共享物理常数来源；
- 本次计算允许改变和必须保持不变的语义。

CLI 面向用户的 MeV 参数应显式命名；内部量遵循仓库自然单位约定。非法物理输入不得静默修正。

## 5. 输入配置及优先级

通用原则：

1. 先确定仓库内默认配置或命名 profile；
2. 再加载专题配置文件；
3. 最后应用显式 CLI 覆盖；
4. 运行产物必须保存最终有效配置，而不能只保存原始模板。

当工作流实际优先级不同，以专题 SOP 和代码实现为准。正式运行不应依赖未记录的 shell 历史或本地临时修改。

## 6. 环境与版本冻结

至少记录：

```powershell
julia --version
git rev-parse HEAD
git status --short
```

要求：

- 从仓库根目录使用 `--project=.`；
- CI 当前固定 Julia `1.12.5`，根项目兼容下限为 Julia `1.10`；
- 稳定 CLI 优先使用 `scripts/dev/run_with_sysimage.ps1` 或 `.sh`；
- sysimage 必须通过 wrapper 的版本、平台和 git commit 校验；
- 正式结果必须能定位到 git commit、argv、有效配置和输入数据版本。

工作树不干净不必自动阻断探索运行，但正式产物必须记录差异，并避免把未提交代码生成的结果误认为可长期复现基线。

## 7. Smoke 预检

Smoke 只验证：

- 环境能加载；
- CLI 参数能解析；
- 主链能运行最小输入；
- 约定产物和 manifest 能生成；
- 错误能够被显式报告。

Smoke 不验证：

- 网格或积分收敛；
- 全参数区间稳定性；
- 与文献或正式 baseline 一致；
- 论文结论可成立。

优先执行专题 SOP 登记的 smoke 命令，并保存输出到独立目录，禁止覆盖正式结果。

## 8. 收敛性验证

正式计算前必须识别所有可能改变数值结果的离散参数，例如：

- 动量、角度与频率积分节点；
- 截断、积分上限和奇点处理策略；
- T、化学势、密度或 `xi` 网格步长；
- solver 迭代上限、seed/continuation 策略；
- 插值、缓存或 threshold-subtraction 策略。

至少使用两个逐级加密配置比较核心观测量、失败率和分支连续性。收敛证据应写入独立分析/产物目录，并记录比较指标。

不得仅为让测试通过而放宽容差。若调整容差，必须说明误差来源、物理意义和验证证据。

## 9. 正式计算命令

正式命令必须满足：

- 使用专题稳定入口和 wrapper；
- 使用由收敛性证据支持的参数；
- 输出到新的 case 目录；
- 记录失败点和 provenance sidecar；
- 默认不覆盖已有正式 case；
- 运行期间的诊断开关不得静默变成正式口径。

正式运行前应保存最终命令或配置；运行后以 manifest 中的 effective config 为准核对实际生效参数。

## 10. 输出目录与产物合同

- 数值、配置、日志和审计：`data/outputs/results/`。
- 图像和图像 manifest：`data/outputs/figures/`。
- `outputs/` 根目录仅作历史兼容，不作为新正式产物路径。

正式数值 case 至少应包含：

- 主结果 CSV/JSON；
- `effective_config.json` 或 manifest 内等价快照；
- `run_manifest.json`；
- 失败点或失败计数；
- README/summary，说明用途、状态和已完成验证；
- 图像存在时的 `plot_manifest.json`。

## 11. Regression / Validation 验收

- Unit：公式局部关系、参数校验和纯函数行为。
- Integration：跨模块主链、CLI 和产物合同。
- Regression：当前代码与仓库内部 baseline 的数值漂移。
- Validation：当前代码与外部实现、文献或独立参考数据的一致性。
- Benchmark：性能变化，不承担正确性验收。

专题 SOP 必须指明最小适用层。若数值结果可能漂移而未运行 regression，必须明确原因和风险。

## 12. 失败点、断点续算与重跑

- 保留失败点坐标、异常类型和必要 solver diagnostics；
- resume 前核对现有输出 header 与有效配置是否兼容；
- 参数或口径改变时创建新 case，不在旧 CSV 上追加；
- `overwrite` 只用于明确的重建，不用于掩盖失败；
- 部分成功产物默认保持 diagnostic-only，直到失败点处置和完整审计完成。

## 13. Diagnostic 与 Formal Production 的边界

以下任一情况存在时，产物只能标为 diagnostic-only：

- 使用 smoke/低精度参数；
- 尚无收敛性证据；
- 非默认 solver、propagator、cache 或插值策略仍在诊断；
- 存在未解释失败点、分支跳变或非有限值；
- provenance 不完整；
- 对应 regression/validation 未执行或失败。

升格 formal production 至少需要：

1. 收敛性证据；
2. 完整 provenance；
3. 失败点审计；
4. 对应 regression，必要时 validation；
5. 人工确认物理口径和图像/表格含义。

## 14. 后处理与作图

后处理应只消费冻结的正式输入，不在绘图脚本中隐式重跑物理求解。图像写入 `data/outputs/figures/`，并通过 manifest 反向定位数值来源和绘图参数。

如果后处理发现异常，应回到 diagnostic 分支产生新证据，不直接手改正式 CSV 或图像。

## 15. 关联公式、API 和测试

- [测试组织与入口规范](../../dev/testing_governance.md)
- [Baseline Version Management](../BASELINE_VERSION_MANAGEMENT.md)
- [脚本入口清单](../scripts/README.md)
- [API 文档总入口](../../api/README.md)
- `scripts/dev/check_docs_consistency.jl`
- `scripts/dev/check_data_output_path_guard.jl`

## 16. 最后验证记录

- 验证日期：2026-07-10
- 验证范围：目录职责、公共门禁和当前稳定入口规则
- 执行命令：`julia --project=. scripts/dev/check_sop_governance.jl`
- 状态：通过；authority map、必需章节、权威唯一性和复核周期检查为绿色
- 备注：公共 SOP 不替代专题数值 smoke 或 convergence evidence
