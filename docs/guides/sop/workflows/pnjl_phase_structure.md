# PNJL 相结构计算 SOP

## 1. 目的与适用范围

本 SOP 用于通过 `Models` 统一相结构产线执行 PNJL/RPNJL 的 T-ρ 扫描、一阶相变边界、spinodal、CEP 和可选 crossover 计算，并生成带 provenance 的成套产物。

适用于 smoke、研究扫描和收敛性证据支持后的正式相结构 case。

## 2. 非适用范围

- 不负责 transport、介子质量或介子数密度计算。
- 不把单次 smoke 结果作为物理基线。
- 不使用已移除的历史 PNJL 内部主线实现。
- 不自动批准 `--promote_reference`；reference 晋升需要独立审阅。

## 3. 权威入口

- 稳定 CLI：`scripts/pnjl/calculate_phase_structure.jl`
- 统一实现：`Models.run_phase_pipeline`
- CLI 支持：`scripts/pnjl/phase_cli_support.jl`
- 默认配置：`config/models/pnjl/phase_pipeline_default.toml`
- 稳定脚本白名单：[脚本入口清单](../../scripts/README.md)

Windows 和 POSIX 均优先通过 `run_with_sysimage` wrapper 启动。

## 4. 物理口径、单位与参数约束

- `T_min/T_max/T_step`：MeV。
- `rho_min/rho_max/rho_step`：`rho/rho0`。
- dense-reference builder 始终把声明的 `T_max`、`rho_max` 和 ranged `xi_max` 作为末端锚点；当步长不能整除区间时，最后一个间隔可以短于名义步长。
- GitHub `workflow_dispatch` 的顶层输入受 25-property 限制。dense-reference workflow 保留常用轴、节点和开关作为独立输入；积分容差、几何量门限、refine level、CEP width 与显式 crossover 温度上限通过 `advanced_config_json` 传入，并由 `resolve_dense_reference_action_config.py` 白名单校验后转换为 CLI 参数。
- `xi`：无量纲各向异性参数。
- PNJL phase 标量热项中，角度只通过 RS 分布自变量 `E_xi` 进入；chi、Polyakov 势、vacuum、Maxwell、spinodal 和 CEP 不额外引入角核。`E_xi` 不作为物理色散关系。
- `thermo_quadrature_policy=tensor_gauss` 保留固定有限区间兼容路径；`rs_reduced_adaptive` 使用 RS 角约化和无穷径向自适应积分，仅适用于上述标量热核，不得外推到 magnetic 或 transport。
- 输出中的化学势字段必须根据列名判断 `mu_q` 或 `mu_B`，两者满足 `mu_B = 3 mu_q` 的项目约定时才能转换。
- 模型参数必须来自 `config/models/<model>/`；共享物理常数来自 `config/physics/`。
- 默认 solver 主线为 `models`，不得把历史兼容实现重新写成当前入口。

本 SOP 不改变 mixed-meson 治理和 non-fixedmu unified joint-solve 语义。

## 5. 输入配置及优先级

实际顺序是：

1. 显式 `--config=...`；若未提供，则按 `--model_kind` 自动加载 `config/models/<model>/phase_pipeline_default.toml`；
2. `--preset=smoke` 覆盖为轻量研究配置；
3. 其后的显式 CLI 参数具有最高优先级。

最终采用的参数写入 `run_manifest.json.effective_config`。正式 case 应保留 `config_path`、`config_hash`、argv 和 git commit。

## 6. 环境与版本冻结

在仓库根目录记录：

```powershell
julia --version
git rev-parse HEAD
git status --short
```

稳定运行优先使用：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 <script> <args>
```

sysimage 不匹配时遵循 wrapper 的 `rebuild` 默认策略；正式 case 不应绕过版本与 commit 校验后继续复用旧 sysimage。

## 7. Smoke 预检

Windows：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/sop_smoke/pnjl_phase
```

Linux/macOS：

```bash
sh scripts/dev/run_with_sysimage.sh scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/sop_smoke/pnjl_phase
```

该命令验证 CLI、Models 相产线、最小扫描和 manifest，不验证 CEP 精度、全温区相边界或网格收敛。

对应 integration smoke：

```powershell
julia --project=. -e 'ENV["INTEGRATION_FILES"]="models/test_phase_cli_smoke.jl"; include("tests/integration/runtests.jl")'
```

## 8. 收敛性验证

正式计算前至少比较：

- 固定网格路径的 `p_num` 与 `t_num`，以及自适应路径的 `thermo_quadrature_rtol/atol/maxevals`；
- 低温费米面邻域的固定网格 oracle 与 RS-reduced adaptive 结果；
- `T_step`；
- `rho_step`，尤其 CEP 和 spinodal 邻域；
- solver `iterations` 与失败/unknown 比例；
- CEP `cep_tol`、refine level 和 adaptive-rho 参数；
- `rho_geometry_convergence` 下 Maxwell/spinodal 的位置量、密度量和面积残差粗细网格误差；
- `adaptive_temperature` 下相线中点相对端点插值的误差与最大细化层；
- dense reference 启用 adaptive xi 时的区间中点误差；
- seed/continuation 与 reverse-rho 方向的分支稳定性。

核心比较量包括 CEP 坐标、first-order boundary、spinodal、crossover、失败率、unknown rate 和 Maxwell area residual。只有这些量在证据支持的精度内稳定后，才能确定正式参数。

严格 `T=0` 只在固定状态热核/数密度层定义；五变量 PNJL phase solve 要求 `T>0`。全温区 reference 的最低正温必须单独做积分、求解和相线三级收敛，不能把接近零温结果标成严格零温。

## 9. 正式计算命令

正式运行使用命名 case 目录，不复用 smoke 输出：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --config=config/models/pnjl/phase_pipeline_default.toml --model_kind=PNJL --mode=production --output_dir=data/outputs/results/pnjl/phase/<case_id>
```

命令中的默认配置只表示当前仓库模板。若收敛性证据要求不同节点或步长，应使用已审阅的命名配置或显式 CLI 覆盖，并以 manifest 的 effective config 为实际口径。

production CLI 可显式配置 `crossover_T_max_MeV`、rho 几何量收敛门限和温度中点自适应门限。
`crossover_T_max_MeV=NaN` 表示继承 `T_max`；不得依赖历史隐藏的 `220 MeV` 截断。启用
`crossover_mu0_only=true` 时，实际 crossover 求解与产物都只能包含 `mu_q=0`，不能只修改 manifest 标记。
`rho_geometry_convergence=true` 时，`cep_max_refine_level` 至少为 1。正式 dense reference 默认启用
rho 几何量与温度中点收敛；adaptive xi 为显式 opt-in，且每个 xi 区间至少求一次中点后才能估计误差。

默认不要传 `--promote_reference`。晋升 reference 必须在产物审计、regression/validation 和人工审阅后单独执行。

## 10. 输出目录与产物合同

目标目录至少包含：

- `trho_scan.csv`；
- `first_order_boundary.csv`；
- `spinodal.csv`；
- 可选 `crossover_line.csv`；
- `phase_grid_convergence.csv`；
- `phase_summary.json`；
- `phase_report.md`；
- `run_manifest.json`。

`run_manifest.json` 至少记录 preset、argv、config path/hash、git commit、effective config、run ID 和 artifact paths。产物路径必须使用仓库相对正斜杠表示，避免机器相关绝对路径泄漏。

production 的 `production_eval/` 保存每个温度的全部 rho 粗细层；聚合 `trho_scan.csv` 只包含各温度最终采用层。
需要审计网格演化时联合读取 `production_eval/` 与 `phase_grid_convergence.csv`，下游相线消费方使用聚合 CSV，
不得把不同层的重复采样自行拼接成一条曲线。

GitHub Actions 的 dense reference workflow 使用“一 xi 一 shard”，`fail-fast=false`。启用 adaptive xi 时，初始锚点和
第一层必需中点并行计算；每层 merge 后复用 `Models._phase_result_midpoint_error` 审计区间，只为未收敛子区间创建
下一层 midpoint shard。Action 入口当前最多支持三层 xi refinement；更深层级必须先扩展 workflow DAG，不能把多个
xi 重新塞回同一 job。每层 matrix 显式限制 `max-parallel=20`，收敛档位按 C0、C1、C2 顺序触发，避免多个完整
run 同时占满 runner 配额；失败后只重跑失败的单 xi shard。

每个 shard 的正式计算 step 具有显式五小时上限，并以一分钟 heartbeat 报告 elapsed time 和已落盘文件数；job 失败时
单独上传 diagnostic artifact。最终 merge 要求所有 shard 绑定同一 commit 与同一非-xi 配置，按物理键确定性排序、
去重，并把各层 xi convergence record 合并进 `phase_grid_convergence_<tag>.csv`。
merge 产物通过 validator 后仍只是候选 artifact，不自动写入或晋升 canonical reference。

dense-reference CSV 必须遵循 RFC 4180 字段转义：字段含逗号、双引号或换行时必须加双引号，字段内双引号写成
两个双引号。shard validator 在上传前逐行拒绝列数不匹配，并报告文件、物理行号、期望/实际列数；merge 在消费
任何 shard 时执行同一结构检查，不能把 malformed row 推迟到最终 `DictWriter` 才失败。

数值计算与后处理使用双 provenance。聚合 manifest 的 `provenance.calculation_git_commit` 绑定生成 shard 的计算
提交，`provenance.postprocess_git_commit` 绑定 merge/validator 提交，并可记录 `source_workflow_run_id`。只有 solver、
数值内核、grid/refinement 语义或 effective config 改变时才必须重算对应 C0/C1/C2 档位。完整 shard artifacts 仍在
保留期内且变更只涉及 CSV、merge、validator 或 Action 后处理时，可通过
`.github/workflows/pnjl-dense-reference-replay.yml` 重放该档位；replay 必须显式给出原 run ID、原 calculation SHA、
tag 与必需 xi anchors，SHA 不匹配即失败。replay 产物始终是 diagnostic candidate，不自动晋升 reference。

主 dense-reference workflow 的可选 `calculation_ref` 只接受不可变的小写 40 字符 Git SHA；留空时数值 shard 与
workflow head 相同。设置该输入时，所有 initial/level-2/level-3 shard checkout 该计算提交，assessment 和 final merge
显式验证同一 `calculation_git_commit`，而调度、merge 与 validator 继续记录当前 workflow head 为
`postprocess_git_commit`。这允许 C0/C1/C2 在仅 CI 后处理修复后仍严格比较同一计算实现，禁止传入 branch、tag 或短 SHA。

跨 GitHub Actions rerun attempt 下载 artifact 时，不能依赖仅对当前 attempt 有效的内部 runtime artifact token。
所有跨 job 或跨 run 的 pattern 下载均显式授予 `actions: read`，设置 `GH_TOKEN`/`GH_REPO`，并通过
`gh run download <run-id> --pattern <pattern>` 获取 artifact。该 CLI 路径会分页枚举整个 run，不能退回只覆盖前/后
100 个 artifact 的单页下载实现。failed-only rerun 后若数值 shard 已完整、但 staged assessment 尚未完成，可使用
`.github/workflows/pnjl-dense-reference-resume.yml` 从源 run 恢复后续 xi refinement；该入口只接受 source run ID、
source calculation SHA、tag 和 initial xi anchors。resume planner 必须从 source manifest 重建 effective CLI config，
验证源 grid 精确包含 anchors 与 level-1 midpoints，并拒绝 calculation SHA、tag、配置或 grid 不一致。

resume 新生成的 level-2/level-3 shard 必须 checkout 原 calculation SHA；修复提交只提供 workflow 调度与后处理代码。
最终 manifest 同时记录原 calculation SHA、新 postprocess SHA 和 source workflow run ID。resume 输出始终是
diagnostic candidate，不自动晋升 reference；若 solver、数值内核、grid/refinement 语义或 effective config 已改变，
禁止使用 resume，必须从受影响的最早档位重新计算。

若完整 source run 已生成所有数值 shard、仅最后一级 assessment/merge/validator 未完成，优先使用 postprocess replay
的可选 final-assessment 输入：显式给出 assessment level、intervals JSON 和原档位四项 xi 容差。replay 的 assessment
只读取 source run shard，随后 final merge 同时读取 source run 已有 convergence records 与 replay run 新生成的最终
record；不得在 replay 中生成新数值 shard。输出继续标记为 diagnostic candidate，并保留 calculation/postprocess 双 SHA。

## 11. Regression / Validation 验收

最小层级：

- Integration：`tests/integration/models/test_phase_cli_smoke.jl`；
- Regression：`tests/regression/phase/test_phase_pipeline_regression.jl`；
- Regression consistency：`tests/regression/models/test_phase_pipeline_consistency.jl`；
- 外部/文献接受性验证按具体物理目标选择 `tests/validation/pnjl/`。

改变网格、CEP 策略、solver 或 seed 语义时，除目标变化外还要验证原有统一入口和非变更约束仍成立。

## 12. 失败点、断点续算与重跑

- smoke 失败时先检查配置路径、sysimage 和 CLI 参数；
- solver unknown 或失败比例异常时保留 `phase_summary.json` 和 report，不直接删点；
- 改变网格或求解策略时使用新 case 目录；
- 后处理修复优先重放已完成且 SHA 已核验的 shard artifacts；不得用 replay 掩盖数值代码或 effective config 变化；
- failed-only rerun 后 assessment 因跨 attempt 下载失败时，先用 tokenized REST 路径核验 source artifacts，再按上述
  source-run resume 合同恢复 refinement；不得为纯下载故障机械重算完整 initial grid；
- run 内 artifact 数量超过 100 时，检查下载日志中的实际 artifact 目录数；若与 GitHub run artifact 总数不符，先按
  分页 `gh run download` 合同修复/重放后处理，不得把缺页误判为物理网格缺失；
- 不在已有正式目录上用不同 effective config 覆盖运行；
- CEP 未找到不自动等价于“物理上不存在 CEP”，必须结合扫描覆盖和诊断字段判断。

## 13. Diagnostic 与 Formal Production 的边界

以下产物只能视为 diagnostic-only：

- `--preset=smoke` 输出；
- `mode=research` 且使用低精度网格；
- CEP direct/interpolate 策略仍在比较；
- unknown rate、失败率或分支连续性未通过审阅；
- 使用 `--promote_reference` 之前未完成回归和验证。

Formal production 需要网格/积分/solver 收敛证据、完整 manifest、失败审计和对应 regression；涉及文献结论时还需要 validation 或明确的证据边界。

## 14. 后处理与作图

绘图应消费冻结的 CSV/JSON，不在绘图脚本中隐式重跑相产线。图像写入 `data/outputs/figures/pnjl/` 的对应 case，并记录数值来源和绘图参数。

## 15. 关联公式、API 和测试

- [Phase API](../../../api/models/phase/README.md)
- [Phase Overview](../../../api/models/phase/Overview.md)
- [扫描 API](../../../api/models/scans/README.md)
- [脚本入口清单](../../scripts/README.md)
- `tests/integration/models/test_phase_cli_smoke.jl`
- `tests/regression/phase/test_phase_pipeline_regression.jl`

## 16. 最后验证记录

- 验证日期：2026-07-18
- 验证范围：CLI preset、默认配置加载、rho/T 网格误差参数、显式 crossover 温区、manifest 和 integration smoke 合同
- 执行命令：见第 7 节及 authority map
- 状态：本 PR 聚焦验证通过；完整 CI 结果以对应 PR checks 为准
- 备注：冷启动验证耗时约 12 分钟；smoke 输出不得升格为正式相结构基线
