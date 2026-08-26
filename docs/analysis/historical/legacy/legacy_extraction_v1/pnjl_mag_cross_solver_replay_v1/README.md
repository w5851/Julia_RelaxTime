# `pnjl_mag` / Julia 跨求解器平衡态对照 v1

本目录记录 `pnjl_mag@e1fc81d3c3c9d220c49972e54307b66a604cb9db` 的 9 个外部代表点与当前 Julia MFIR 磁场求解器的轻量跨求解器诊断。它不是新的生产扫描，也不创建 validation acceptance target。

## 口径

- 输入点来自 [`pnjl_mag_equilibrium_replay_v1`](../pnjl_mag_equilibrium_replay_v1/)，固定 `muB=0`，`T={50,150,240} MeV`，`eB={0.2,0.4,0.8} GeV^2`；
- Julia 求解使用与外部 kernel oracle 相同的参数映射：`hc=197.33`、`a=0.0108805`、MFIR/Hurwitz-zeta、`n=0..79`；
- 每个点使用 3 个 seed：外部 continuation 状态、手征方向扰动、项目热 seed；记录 trust-region primary、Newton fallback、候选去重和 residual；
- `full_screening` 使用 `p_num=64, zeta_num=64, pz_max=25 fm^-1` 覆盖全部 9 点；
- `high_matched` 只复核 `T=240 MeV` 的 3 点，使用 `p_num=128, zeta_num=256, pz_max=40 fm^-1` 与外部 replay 口径匹配。

## 观察结果

9 点筛查中所有点都找到有限且物理的 Julia 候选，候选 residual 最大为 `4.67e-12`。筛查节点下，外部状态到 Julia 候选的最大状态距离为 `2.34e-3`；该组的 Omega 差最多 `1.84e-4`，不能脱离节点差异解释为模型差异。

高温匹配节点复核中，3 个外部状态在 Julia 路径上的 residual 最大为 `3.05e-12`，Omega 差最大为 `1.18e-12`。这表明在明确固定参数和积分节点后，MFIR kernel 与外部 equilibrium 分支在这些点上数值一致。

`T=240 MeV, eB=0.8 GeV^2` 的 3-seed 复核找到两个物理候选：外部 continuation 分支的
`Omega=-24.1229642024 fm^-4`，以及另一个候选的 `Omega=-23.5233009589 fm^-4`。这只证明候选集合中存在额外局部根，不证明已经枚举全部分支；按当前磁场约定，convenience state 仍选择较低 Omega 候选，branch-aware 输出保留两者。

## 结论边界

本诊断支持以下有限结论：外部 9 点可被当前 Julia MFIR 路径重放；节点匹配后外部状态和 Omega 在所测点通过数值对照；高温高场处确实需要保留多候选信息。

它不支持以下结论：全域 branch completeness、外部 continuation 与 Julia seed governance 的完全等价、默认生产节点的收敛证明、或正式 validation target admission。当前仍需独立的积分收敛矩阵和明确的分支覆盖合同。

生成命令：

```powershell
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_julia_equilibrium.jl --scope full --nodes screening
julia --startup-file=no --project=. scripts/analysis/pnjl/compare_pnjl_mag_julia_equilibrium.jl --scope high --nodes matched
```

机器可读结果、输入 hash 和脚本 hash 见 [`manifest.json`](manifest.json) 与 [`provenance.json`](provenance.json)。
