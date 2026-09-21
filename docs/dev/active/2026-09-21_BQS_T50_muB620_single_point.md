# T=50 MeV、muB=620 MeV 的 BQS charged GBU 单点

## 范围

- 分类：research；独立单点核验，不替换现有 primary track。
- 用户请求：沿用现有计算方法，实际求解极端温度/化学势下的介子数密度比，而非灵敏度外推。
- 用户授权远程计算，要求 main 最终不改变。本次采用更严格边界：main 全程不修改；仅推送会话分支 `arena/01a0c36e-julia-relaxtime`，由 branch-only push workflow 运行。
- 计算前 main SHA：`0ef8ba0918c94a4b3d3bf7dbdeb2278793447c57`。

## 已定位的路径与量的定义

1. `scripts/relaxtime/run_charged_gbu_freezeout_scan.jl` → `Models.run_charged_gbu_freezeout_scan`。
2. `src/models/workflow_apps/ChargedGBUResearchWorkflow.jl`：`background`、`channel_density`、`run_job`、`ratio_row`。
3. `config/models/pnjl/charged_gbu_infinite_v1.toml`：quark-only BQS，rhoQ/rhoB=0.4、rhoS=0、xi=0、无限热积分、无介子反馈。
4. `scripts/analysis/relaxtime/causal_gbu_infinite_{qgate,yield,profile,thermal}.jl` 等现有单一数值核。

现有 `freezeout_20260907_v2/ratios.csv` 的前三点温度/化学势与用户报告一致：3 GeV 为 T=79.95694382094246 MeV、muB=719.0764156129742 MeV；正电比值 0.5058610078337089。5 和 7.7 GeV 的 T、muB 亦一致。

本次直接指定 T、muB，不修改冻结线参数，也不把此点标成新的 3 GeV 冻结线预测。主观察量是 `n(K+)/n(pi+)`，同时计算负电对照；是模型固定背景上的 GBU 部分产额，不是含衰变、全反馈的最终实验产额。

## 实施与复现

新增 analysis-only adapter：`scripts/analysis/relaxtime/compute_charged_gbu_bqs_point.jl`。它调用现有 workflow 内部函数，不复制/修改物理核或容差；稳定 API 和旧入口不变。

```sh
bash scripts/dev/run_with_sysimage.sh --mismatch-policy=fallback \
  scripts/analysis/relaxtime/compute_charged_gbu_bqs_point.jl \
  --T-MeV 50 --muB-MeV 620 --workers 2 \
  --output data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/point_T50_muB620_bqs_20260921
```

输出目录必须不存在。使用历史 3 GeV 背景仅作为初始种子，重新联立求解八维背景；另做 cold seed 及 48→96 节点重新求解，固定背景以 48/96/192 节点评估净 BQS 密度。四通道复用局部门槛、Mott 分段、q 加阶和尾带门槛；失败原样保存，不能报告为通过的 ratio。

本地无 Julia，安装站与二进制下载 TLS 被阻断；验证和计算由 `.github/workflows/arena-bqs-single-point.yml` 在 Julia 1.12.5、根 Manifest 环境运行，artifact 保留日志、来源快照、背景、约束、密度、比值和门槛。

验证层：现有 BQS unit、GBU production contract unit、新 adapter unit、现有 charged GBU entrypoint integration，以及单点的背景/积分数值收敛。未修改 src/config/baseline，未以全仓回归代替此点的数值验收。

## 执行结果

实际计算 run **35587522965**（源码 `4ff2587`），不是灵敏度外推：

- n(pi+) = 5.8481395635567366e-6 fm^-3。
- n(K+) = 3.7183733172004564e-7 fm^-3。
- **K+/pi+ = 0.06358215765526304，建议报告 0.0636。** 比 0.06 高 5.97%，未进入 0.05–0.06，但已接近上沿。
- muQ = -2.39578390849 MeV，muS = 103.631411979 MeV；rhoQ/rhoB = 0.400000000009782、rhoS = -3.94e-21 fm^-3、背景残差 4.59e-15。
- 48→96 节点背景重解一致，48/96/192 净密度一致到浮点精度。pi+ 最终 q32、加阶差 0.0173%；K+ 最终 q16、加阶差 0.1001%。两通道均通过既有全部积分门槛。

### 失败记录及复核

原运行整体退出 1，保留 `fixed_point_with_failed_checks`，不可声称四通道全通过：

1. K- `eta_passed=false`，不输出其比值，仍保留 NaN。
2. 新增 cold-vs-continuation 检查失败：默认 cold 解虽被原绝对容差接受，mu_s 仍相差约 5.6 MeV。原续算背景的净荷残差与节点检查没有失败。

更严格的同方程冷启动复核 **35589706534**（源码 `219ee25`）通过：使用原 `Models.build_residual!` 的八维联立方程，诊断 `ftol/xtol` 收紧至 1e-14。cold 解继续收敛后，与原续算解最大分量差 1.13454e-12；续算解原本即满足更严格标准，迭代 0 次。该差异是低密度下的停止精度问题，不是这两次求解的不同解支；不做全局唯一性宣称，不改 production solver，不改原失败记录。正电比值由此获得背景复核支持。

### 结果保存与审计

- 目录：`data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/point_T50_muB620_bqs_20260921/`。
- 主报告：该目录 `README.md`，含原始 CSV、精选背景、门槛摘要、冷启动复核和复现命令。
- 完整计算 artifact ID `10633463464`，ZIP SHA256 `0f9ca0a580c2fc06ae676390cf8d207cf11978ba5545878609d3425c8c296c38`，远端保留 30 天。
- 沙箱无法直接访问 artifact blob。因此增加只读结果审计 workflow（run **35589536934**），在 runner 内验证原 manifest 的全部输出 SHA256 并独立求和，通过 GitHub check API 返回精选证据。该流程不修改源码或 main。
- 三个合格通道独立壳层求和与原密度最大差 8.47e-22 fm^-3；本地再验证正电比值等于原密度相除，并验证所保存原始文件哈希及 nominal BQS 约束。
- 前置验证：BQS unit 39/39、GBU contract unit 24/24、新 adapter unit 14/14、入口 integration 12/12，共 **89/89**。没有修改 baseline，也未声称运行全仓回归。
- `src/`、`config/`、Project/Manifest 和 baselines 相对初始提交均无改动。最终 main SHA 再次核验仍为 `0ef8ba0918c94a4b3d3bf7dbdeb2278793447c57`。

结论边界：本次完成用户指定极端点的正电模型部分产额实算。该点不在原默认冻结线上，不能单独证明整个统计热模型可信区间内有解或无解；quark-only BQS 不能与包含介子反馈的总守恒荷背景混同。
