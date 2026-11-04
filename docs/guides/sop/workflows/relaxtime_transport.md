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
- `isotropic` propagator 仍是显式诊断选项。旧 `validated_anchored` 缓存已退役，当前生产只使用 thermal-weighted。

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
- profile 中的动量上限和映射尺度单位为 `fm^-1`；末端默认采用半无穷积分。
- 输出中 `tau` 为 fm，`eta/zeta` 为 `fm^-3`，`sigma` 为 `fm^-1`；派生无量纲比值按列合同解释。
- 默认 `propagator_xi_policy=match_thermo`、`sampling_profile=thermal_weighted_v1`。

所有节点数和步长必须为有效正值；非有限 tau 或上游失败不得静默进入输运积分。

## 5. 输入配置及优先级

`run_gap_transport_scan.jl` 以 CLI 选择扫描坐标，数值参数默认由 `config/models/pnjl/transport_sampling/thermal_weighted_v1.toml` 统一配置。与 profile 冲突的参数会被拒绝；有效参数写入 `effective_config.json` 和 `run_manifest.json`，数值源码和逐道预算另有 sampling sidecar。

与单点 workflow 相关的物理默认值仍可能来自 `config/physics/<profile>.toml`。正式 case 必须同时记录 physics profile、CLI argv 和 effective config，不能只依赖默认值的文字描述。

Issue #130 runtime 使用 accepted reference，显式 strict 模式使用 strict layer。
旧参考仅提供起始提示，稳定分支和 mode-A 温度锚点由当前 profile 的背景模型重算。
legacy runtime 已退役，诊断视图需显式提供 reference root；render 层不作为正式 runtime 输入。

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

需要检查实际数值主链时，执行同一 profile 的单点 pilot。Windows 示例：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_gap_transport_scan.jl --tmin 150 --tmax 150 --tstep 1 --mubmin 0 --mubmax 0 --mubstep 1 --xi-list 0.0 --no-compute-bulk --sampling-profile thermal_weighted_v1 --output data/outputs/results/sop_smoke/relaxtime_transport/scan.csv --failed-points-output data/outputs/results/sop_smoke/relaxtime_transport/failed_points.csv --provenance-dir data/outputs/results/sop_smoke/relaxtime_transport --overwrite
```

该 pilot 使用默认采样精度但关闭 bulk，只检查这一点的 η/σ/τ 主链与 provenance。完整资格还需 ζ 和代表点独立收敛检查；只测试参数和输出合同可改用 phase-guided mode-B 的 `--dry-run`。

## 8. 收敛性验证

正式参数至少应覆盖以下加密轴：

- 平衡层 `p_num/t_num/max_iter`；
- tau 层 `tau_p_nodes/tau_angle_nodes/tau_phi_nodes/tau_n_sigma`；
- 热权重缓存的误差预算、阈值/共振发现和独立直算参考；
- A 与末端输运的半无穷映射尺度、径向和角度节点；
- T、`mu_B`、`xi` 网格步长；
- 失败点与 continuation/seed 稳定性。

比较量至少包括质量/序参量、代表通道 sigma、tau 六分量、`eta/s`、`zeta/s`、`sigma/T`、失败率和局部跳变。局部共振或阈值区域应额外补点，不能只看全局最大相对误差。

缓存插值预算不等于完整输运误差。默认 profile 的代表集以最终六种 τ、η/s、ζ/s、σ/T 的 1% 为目标；当前资格及适用范围见[输运采样指南](../../scripts/transport_sampling_qualification.md)。

## 9. 正式计算命令

正式命令使用 case 目录，并显式保存 channel diagnostics、failed points 和 provenance：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_gap_transport_scan.jl <由收敛性证据确定的网格与积分参数> --propagator-xi-policy match_thermo --sampling-profile thermal_weighted_v1 --output data/outputs/results/relaxtime/transport/<case_id>/scan.csv --channel-diagnostics-output data/outputs/results/relaxtime/transport/<case_id>/channel_diagnostics.csv --failed-points-output data/outputs/results/relaxtime/transport/<case_id>/failed_points.csv --provenance-dir data/outputs/results/relaxtime/transport/<case_id>
```

这里不提供通用“正式节点数”，因为正式参数必须来自该 case 的收敛性证据。实际生效值以 `effective_config.json` 为准。

显式 `--resume` 会跳过已有点并校验数值源码/profile 身份；若改变数值口径，应创建新 case，不能向旧配置 CSV 追加。

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
