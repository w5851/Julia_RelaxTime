# BQS 单点：T = 50 MeV，μB = 620 MeV

## 实算结果

沿用项目当前 charged GBU 无限热积分链路，背景为 **quark-only BQS**：

- ρQ/ρB = 0.4；ρS = 0；ξ = 0。
- 直接指定 T、μB，不应用冻结线映射。
- 没有介子反馈、额外 pion fugacity 或末态衰变修正。

**正电部分产额比：**

\[
\frac{n_{K^+}}{n_{\pi^+}}
=\frac{3.7183733172004564\times10^{-7}}{5.8481395635567366\times10^{-6}}
=0.06358215765526304\simeq\mathbf{0.0636}.
\]

密度单位均为 fm⁻³。它比 0.06 高约 **5.97%**，没有进入 0.05–0.06，但已相当接近该区间。此点不是重新拟合后的“3 GeV 冻结点”，也不能单凭这一点断言整个可信 T–μB 区域是否存在匹配解。

## BQS 背景

| 量 | 结果 |
|---|---:|
| μQ / MeV | −2.39578390849 |
| μS / MeV | 103.631411979 |
| μu / MeV | 205.069477394 |
| μd / MeV | 207.465261303 |
| μs / MeV | 103.833849324 |
| μπ⁺ = μu−μd / MeV | −2.39578390849 |
| μK⁺ = μu−μs / MeV | 101.235628070 |
| ρB / fm⁻³ | 1.27531906174×10⁻⁶ |
| ρQ/ρB | 0.400000000009782 |
| ρu/ρd | 0.875000000011465 |
| ρS / fm⁻³ | −3.94137260617×10⁻²¹ |
| 八维解残差范数 | 4.59254818221×10⁻¹⁵ |

μS（守恒奇异化学势）与 μs（s 夸克化学势）不是同一个量。BQS 映射采用 μB=μu+2μd、μQ=μu−μd、μS=μd−μs；没有假设三个夸克化学势都等于 μB/3。

48→96 动量节点重新求解的背景解完全一致；对该背景以 48、96、192 节点计算净密度也一致到浮点精度。详细值见 `bqs_diagnostics.csv`。

## 积分验收与失败边界

| 通道 | 密度 / fm⁻³ | 采用的 q 阶数 | 相邻 q 加阶相对差 | 状态 |
|---|---:|---:|---:|---|
| π⁺ | 5.84813956356×10⁻⁶ | 32 | 0.0173151% | accepted |
| K⁺ | 3.71837331720×10⁻⁷ | 16 | 0.100084% | accepted |
| π⁻ | 6.43901165330×10⁻⁶ | 32 | 0.0176901% | accepted |
| K⁻ | NaN | — | — | local_gate_failed |

正电两通道通过原有局部核、围道、η 序列、内层加阶、Mott 分段、q 加阶和尾带门槛。上述相邻加阶差是数值诊断，不是完整误差条或模型误差。独立 Python 重算壳层加权和，与保存的三个合格密度之差均不超过 8.5×10⁻²² fm⁻³。

**不得把本次整体运行标为全通道通过：**

1. K⁻ 的 `eta_passed=false`，因此 K⁻/π⁻ 按原合同保留 NaN，不用估计值替代。
2. 新增的冷启动一致性检查首次失败：原求解器接受的 cold 解残差约 9.60×10⁻¹⁰，s 夸克化学势与高精度续算解相差 0.0283832 fm⁻¹（约 5.60 MeV）。这发生在非常低的净密度下；续算背景本身的 BQS 残差与加阶检查见上表。更严格的同方程冷启动复核单独保存，不改写原始失败记录。

因此原始 `manifest_status` 为 `fixed_point_with_failed_checks`，Actions 计算步骤退出码为 1；`ratios.csv` 的 `plus_passed=true` 仅是正电通道积分合同，不能读成所有额外背景/负电检查都通过。

### 冷启动问题的后续闭合

[复核 run 35589706534](https://github.com/w5851/Julia_RelaxTime/actions/runs/35589706534) **通过**。复核脚本 `check_charged_gbu_bqs_cold_start.jl` 直接复用现有 `Models.build_residual!` 八维联立方程，仅将这次诊断的 `ftol/xtol` 收紧至 10⁻¹⁴，不修改 production 求解器或物理模型。

- 原 cold 解在 192 节点下 ρQ/ρB=0.4001204113、ρS=1.4961×10⁻¹⁰ fm⁻³，说明原绝对停止标准在如此低的净密度下不足以保证化学势精度。
- 收紧后 cold 解与原续算解最大分量差为 **1.13454×10⁻¹²**；复核残差约 1.07×10⁻¹⁵，ρQ/ρB=0.4000000000000044。
- 用于密度积分的原续算解一开始就满足更严格的标准，复核迭代数为 0。无需替换背景或重算一个不同物理点。

因此本报告的 **正电比值 0.0636 有经过复核的 BQS 背景与合格的两通道积分支持**。冷启动差异已在同方程精度层面解释并闭合；这不构成全局解支唯一性证明，也不替代 K⁻ 的失败门槛。原始失败文件仍保持原字节，后续证据另存 `cold_start_verification.json`。该文件来自提交 `219ee25` 的 GitHub check `106302708540`。

## 来源、保存与复现

- 数值源码提交：`4ff25878932c1997748ad994e918bc2283396bc3`。
- Julia 1.12.5，根 Project/Manifest。
- [原始计算 run 35587522965](https://github.com/w5851/Julia_RelaxTime/actions/runs/35587522965)。
- [完整 artifact 10633463464](https://github.com/w5851/Julia_RelaxTime/actions/runs/35587522965/artifacts/10633463464)：包含原始日志、全部源快照、壳层及校验侧文件，保留期 30 天。
- artifact ZIP SHA256：`0f9ca0a580c2fc06ae676390cf8d207cf11978ba5545878609d3425c8c296c38`。
- [独立只读结果审计 run 35589536934](https://github.com/w5851/Julia_RelaxTime/actions/runs/35589536934)：原 manifest 列出的所有输出哈希验证通过。
- `remote_audit.json` 保存审计 API 返回的背景、通道门槛摘要和独立求和。此目录只保存小型精选输出，不伪装为完整 artifact；三个 CSV、三个背景 JSON 和 `background_checks.json` 与原始文件 SHA256 一致。
- 前置验证：BQS unit 39/39、GBU contract unit 24/24、adapter unit 14/14、entrypoint integration 12/12，共 89/89。

计算路径：

```text
Models → ChargedGBUResearchWorkflow.background
       → FixedMuBConservedCharges 联立八维求解
       → ChargedGBUResearchWorkflow.channel_density / run_job
       → causal_gbu_infinite_* 原有核和门槛
       → ratio_row
```

单点 adapter：`scripts/analysis/relaxtime/compute_charged_gbu_bqs_point.jl`。复现（输出目录须不存在）：

```sh
bash scripts/dev/run_with_sysimage.sh --mismatch-policy=fallback \
  scripts/analysis/relaxtime/compute_charged_gbu_bqs_point.jl \
  --T-MeV 50 --muB-MeV 620 --workers 2 \
  --output data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/my_T50_muB620
```

所有远程操作只在本会话分支运行，没有合并或修改 main；src、物理配置、旧入口、数值容差和 baseline 均未修改。
