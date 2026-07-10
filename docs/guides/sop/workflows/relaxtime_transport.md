# Relaxtime transport 计算 SOP

## 1. 目的与适用范围

本 SOP 用于执行 PNJL 平衡态、散射截面/平均散射率、弛豫时间和 RTA 输运系数的联合扫描，并生成可恢复、可审计的结果与 provenance sidecar。

主链是：

```text
equilibrium -> A/B0 -> polarization -> propagator -> |M|^2
-> sigma(s) -> average rate -> tau -> eta/zeta/sigma
```

## 2. 非适用范围

- 不负责 phase-guided transport 的专题 canonical production 口径；该流程仍以现有脚本指南和正式产物 README 为准。
- 不负责 Beth-Uhlenbeck 介子热力学或介子数密度。
- 不把低节点单点 pilot 作为回归 baseline。
- 不把 `isotropic` propagator 或 `validated_anchored` cache 诊断策略自动升格为默认修复。

## 3. 权威入口

- 稳定扫描 CLI：`scripts/relaxtime/run_gap_transport_scan.jl`
- CLI 参数合同：`scripts/relaxtime/gap_transport_scan_cli.jl`
- 单点统一入口：`Models.solve_gap_and_transport`
- 稳定脚本白名单：[脚本入口清单](../../scripts/README.md)
- 端到端语义：[输运主链路](../../../reference/formula/relaxtime/transport/Transport_EndToEnd_Pipeline.md)

正式 CLI 默认优先通过 `run_with_sysimage` wrapper 启动。

## 4. 物理口径、单位与参数约束

- CLI `tmin/tmax/tstep`：MeV。
- CLI `mubmin/mubmax/mubstep`：重子化学势 `mu_B`，MeV。
- 兼容参数 `mumin/mumax/mustep` 表示 `mu_q`，脚本按 `mu_B = 3 mu_q` 转换；新命令优先使用 `mub*`。
- `xi`：无量纲。
- `tr-p-max`：`fm^-1`。
- 输出中 `tau` 为 fm，`eta/zeta` 为 `fm^-3`，`sigma` 为 `fm^-1`；派生无量纲比值按列合同解释。
- 默认 `propagator_xi_policy=match_thermo`、`sigma_cache_policy=default`，保持当前主线语义。

所有节点数和步长必须为有效正值；非有限 tau 或上游失败不得静默进入输运积分。

## 5. 输入配置及优先级

`run_gap_transport_scan.jl` 当前以 CLI 参数为主要扫描合同。未显式传入时使用 `gap_transport_scan_cli.jl` 中的默认值；有效参数写入 `effective_config.json` 和 `run_manifest.json`。

与单点 workflow 相关的物理默认值仍可能来自 `config/physics/<profile>.toml`。正式 case 必须同时记录 physics profile、CLI argv 和 effective config，不能只依赖默认值的文字描述。

## 6. 环境与版本冻结

在仓库根目录记录：

```powershell
julia --version
git rev-parse HEAD
git status --short
```

优先使用 Windows/POSIX wrapper。长扫描前应确认：

- sysimage 与当前 commit 匹配；
- 输出 case 不与旧配置混用；
- 目标磁盘空间足够保存主 CSV、diagnostics、failed points 和日志；
- 正式运行使用已完成收敛性评估的参数。

## 7. Smoke 预检

先验证 CLI 与入口合同：

```powershell
julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --help
julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'
```

再执行单点低成本 diagnostic pilot。Windows 示例：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_gap_transport_scan.jl --tmin 150 --tmax 150 --tstep 1 --mubmin 0 --mubmax 0 --mubstep 1 --xi-list 0.0 --no-compute-bulk --p-num 8 --t-num 4 --max-iter 10 --tau-p-nodes 8 --tau-angle-nodes 2 --tau-phi-nodes 4 --tau-n-sigma 4 --sigma-grid-n 16 --tr-p-nodes 8 --tr-p-max 6.0 --output data/outputs/results/sop_smoke/relaxtime_transport/scan.csv --failed-points-output data/outputs/results/sop_smoke/relaxtime_transport/failed_points.csv --provenance-dir data/outputs/results/sop_smoke/relaxtime_transport --overwrite
```

该 pilot 只验证联合主链和 provenance；节点明显低于正式精度，不得用于 baseline、论文图或物理趋势结论。

## 8. 收敛性验证

正式参数至少应覆盖以下加密轴：

- 平衡层 `p_num/t_num/max_iter`；
- tau 层 `tau_p_nodes/tau_angle_nodes/tau_phi_nodes/tau_n_sigma`；
- `sigma_grid_n` 与插值模式；
- transport 层 `tr_p_nodes/tr_p_max`；
- T、`mu_B`、`xi` 网格步长；
- threshold-subtraction 参数；
- 失败点与 continuation/seed 稳定性。

比较量至少包括质量/序参量、代表通道 sigma、tau 六分量、`eta/s`、`zeta/s`、`sigma/T`、失败率和局部跳变。局部共振或阈值区域应额外补点，不能只看全局最大相对误差。

`propagator_xi_policy=isotropic` 和 `sigma_cache_policy=validated_anchored` 当前属于显式诊断分支；在完成复算、收敛和回归前，不得写成正式修复口径。

## 9. 正式计算命令

正式命令使用 case 目录，并显式保存 channel diagnostics、failed points 和 provenance：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_gap_transport_scan.jl <由收敛性证据确定的网格与积分参数> --propagator-xi-policy match_thermo --sigma-cache-policy default --output data/outputs/results/relaxtime/transport/<case_id>/scan.csv --channel-diagnostics-output data/outputs/results/relaxtime/transport/<case_id>/channel_diagnostics.csv --failed-points-output data/outputs/results/relaxtime/transport/<case_id>/failed_points.csv --provenance-dir data/outputs/results/relaxtime/transport/<case_id>
```

这里不提供通用“正式节点数”，因为正式参数必须来自该 case 的收敛性证据。实际生效值以 `effective_config.json` 为准。

默认 resume 会跳过已有点；若改变任何数值口径，应创建新 case 或显式 `--overwrite` 重建，不能向旧配置 CSV 追加。

## 10. 输出目录与产物合同

case 目录至少包含：

- `scan.csv`；
- `effective_config.json`；
- `run_manifest.json`；
- 建议 `channel_diagnostics.csv`；
- 建议 `failed_points.csv`；
- 正式 case 的 README/summary 和收敛性证据链接。

主 CSV 应保持行级 T、`mu_B`、`xi`、平衡求解状态、tau 与输运字段。缺点和失败点不得仅通过行数减少表达，应在 summary/sidecar 中显式记录。

## 11. Regression / Validation 验收

最小层级：

- Unit：`tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`；
- Integration：`tests/integration/relaxtime/test_transport_workflow_smoke.jl`；
- Regression：`tests/regression/relaxtime/test_transport_fixedpoint_regression.jl`；
- Workflow consistency：`tests/regression/relaxtime/test_workflow_vs_direct_consistency.jl`；
- 外部接受性：按目标选择 `tests/validation/relaxtime/` 中的 tau、sigma 或 transport ratio targets。

影响截面、散射率、tau 或输运系数时，原则上必须运行对应 regression。跳过时要记录原因和数值漂移风险。

## 12. 失败点、断点续算与重跑

- 默认 resume 只适用于相同 header 和相同数值口径；
- header 不兼容时脚本会拒绝继续，应选择新输出或显式重建；
- 使用 `failed_points.csv` 保留坐标、seed/phase diagnostics 和异常；
- 单点失败后可在独立 diagnostic case 加密，不直接手填主 CSV；
- 输出中出现非有限 tau、异常负值或局部大跳变时，停止 formal promotion，先沿上游截面/传播子链诊断。

## 13. Diagnostic 与 Formal Production 的边界

以下产物只能标为 diagnostic-only：

- 第 7 节低节点 pilot；
- 使用非默认 propagator/cache 策略；
- 无收敛性证据；
- failed points 未解释；
- tau 或输运存在未审计非有限值/突变；
- 未运行对应 transport regression。

Formal production 需要：

1. 参数加密与局部阈值窗口收敛证据；
2. 完整 effective config、manifest 和 diagnostics；
3. 固定点/路径 regression；
4. 必要的文献/legacy validation；
5. 人工确认输出单位、物理口径和图像映射。

## 14. 后处理与作图

绘图脚本应只读取冻结 `scan.csv`，图像写入 `data/outputs/figures/relaxtime/transport/<case_id>/`，并生成 `plot_manifest.json`。不得在绘图阶段静默替换失败点、插值跨越缺口或重新计算物理量。

## 15. 关联公式、API 和测试

- [Transport workflow API](../../../api/relaxtime/workflow/TransportWorkflow.md)
- [Transport API](../../../api/relaxtime/transport/README.md)
- [输运端到端链路](../../../reference/formula/relaxtime/transport/Transport_EndToEnd_Pipeline.md)
- [测试治理](../../../dev/testing_governance.md)
- [脚本入口清单](../../scripts/README.md)
- `tests/regression/relaxtime/test_transport_fixedpoint_regression.jl`

## 16. 最后验证记录

- 验证日期：2026-07-10
- 验证范围：CLI 参数合同、Models 入口、provenance sidecar 和现有分层测试映射
- 执行命令：见第 7 节及 authority map
- 状态：通过；入口单测全部通过，单点 pilot `success_count=1`、`error_count=0`、`failed_points=0`
- 备注：pilot 使用 `compute_bulk=false`，因此 zeta 相关字段为 `NaN` 是预期合同；低节点结果保持 diagnostic-only
