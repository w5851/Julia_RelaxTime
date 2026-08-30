# 科学计算 SOP 总览

本目录是 Julia_RelaxTime 面向计算执行、验证、复现与产物升格的权威 SOP 入口。

SOP 回答“如何可靠地运行”；公式、API、研究结论和任务历史仍由既有目录负责。权威身份与复核周期以 `config/governance/docs_authority_map.toml` 为机读单一来源。

## 阅读顺序

1. 先读 [公共科学计算生命周期](common_scientific_run.md)。
2. 再选择具体工作流 SOP。
3. 需要理解公式时转到 `docs/reference/formula/`。
4. 需要调用 Julia API 时转到 `docs/api/`。
5. 需要判断某次结果支持什么结论时转到 `docs/analysis/`。

如果任务属于新领域、新可观测量、新生产候选或引入外部公式，先执行
[公式路线闭合 SOP](formula_route_closure.md)，再进入专题工作流 SOP。该门禁只
保证模型、来源、代码和测试可追踪，不替代专题的数值收敛和 production 审核。

## 当前 active SOP

| SOP | 权威范围 | 执行入口 |
|---|---|---|
| [公共科学计算生命周期](common_scientific_run.md) | 所有正式数值计算的公共阶段与产物升格 | 各专题入口 |
| [PNJL 相结构](workflows/pnjl_phase_structure.md) | PNJL/RPNJL 相结构扫描、CEP、spinodal 与 crossover 产线 | `scripts/pnjl/calculate_phase_structure.jl` |
| [Relaxtime transport](workflows/relaxtime_transport.md) | 平衡态 → 截面/散射率 → τ → RTA 输运扫描 | `scripts/relaxtime/run_gap_transport_scan.jl` |
| [介子热力学](workflows/meson_thermodynamics.md) | 介子 pressure、QP/LD 与 canonical EOS | `scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl`（专题入口） |
| [介子数密度](workflows/meson_density.md) | stable/BW/BU 数密度与 path × regime production | `scripts/relaxtime/run_combined_meson_density_scan.jl`（domain-candidate） |
| [论文级绘图资产与生产](workflows/figure_production.md) | 四层 figure mode、style profile、plot manifest、strict layout gate 与 Origin 边界 | `scripts/plotting/validate_plot_artifact.py` |
| [公式路线闭合](formula_route_closure.md) | 新领域/新公式从微观模型到代码、测试和生产状态的闭合门禁 | `scripts/dev/check_formula_route_closure.jl` |

## 状态语义

- `draft`：结构已建立，但尚未完成最小执行验证，不能作为正式计算依据。
- `active`：入口、产物和最小验证已经核对，可作为当前执行依据。
- `deprecated`：保留迁移说明，但必须指向替代 SOP。

状态不由本页人工重复维护，以 `config/governance/docs_authority_map.toml` 为准。

## 完成语义

- **written**：SOP 文件和索引存在。
- **configured**：入口、配置、单位、输出和验证命令已登记。
- **effective**：最小验证真实可执行并产生约定产物。
- **usable**：SOP 足以指导失败恢复、收敛判断和产物升格。

只有同时达到四层，才能把 SOP 视为“有效可用”。

## 新增或修改 SOP

1. 从 [_template.md](_template.md) 开始。
2. 在 `config/governance/docs_authority_map.toml` 登记唯一 ID 和 `authoritative_for`。
3. 确认稳定脚本已进入 [脚本白名单](../scripts/README.md)。
4. 运行：

```powershell
julia --project=. scripts/dev/check_sop_governance.jl
```

5. 实际执行登记的最小验证命令，再更新 `last_verified`。

## 边界

- 本目录不复制自动生成的 API 全集。
- 本目录不保存一次性运行日志或结果结论。
- smoke 成功不自动等价于数值收敛。
- diagnostic 产物不自动进入论文或正式 baseline。
- 改变物理语义、默认求解策略或容差时，必须另立任务并执行对应 regression/validation 治理。
