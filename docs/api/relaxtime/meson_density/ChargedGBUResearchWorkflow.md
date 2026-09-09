# 无限热 charged GBU smoke production 默认入口

## 授权与默认

`Models.run_charged_gbu_freezeout_scan` 是 charged 介子密度的 smoke production
默认入口。它固定使用 quark-only BQS 背景、无限 PNJL 热差和无介子反馈；通用
legacy `MesonDensity` API 仍保留为显式兼容入口。
它不修改 `Models.run_freezeout_meson_density_scan`、MesonDensity regime 或
PNJLCore；不更新正式 baseline，也不把 smoke default 晋升为 formal baseline。
观察量是固定 quark-only BQS 背景上的 GBU 部分产额，非最终强子产额，
不含介子反馈、额外 pion fugacity 或负谱裁零。

计算核保持已验收的单一源码，由默认 workflow 惰性加载；没有复制一个不同的
生产物理核。该入口的稳定范围是 charged GBU smoke production，不代表全部
analysis API 已获稳定公共 API 授权。运行 manifest 保存实际依赖快照。

作者已审核2026-09-07的v2冻结线图像；10个配置能量点、40个通道全部通过。
这是固定方法/网格的研究产额验收，不是实验定量拟合或连续背景全域认证。

## 调用

```julia
include("src/models/Models.jl")
Models.run_charged_gbu_freezeout_scan(
    output="data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/my_run",
    energies=[3.,5.,7.7,11.5,19.6,27.,39.,62.4,130.,200.],
    workers=3,
)
```

CLI：`scripts/relaxtime/run_charged_gbu_freezeout_scan.jl`；参数包括
`--output DIR`、`--figure-output DIR`、`--config TOML`、`--sqrts-list ...`、
`--workers 1..4`、`--resume` 和 `--no-plot`。默认方法配置为
`config/models/pnjl/charged_gbu_infinite_v1.toml`。
使用项目 sysimage wrapper；元数据不匹配时用 `-MismatchPolicy fallback`，
不自动构建或静默加载旧物理源码的 sysimage。

## 数值与失败合同

- 固定 `FixedMuBConservedCharges`、rhoQ/rhoB=0.4、rhoS=0、xi=0。
  从高能向低能延续背景种子；背景求解48个动量节点、残差<=1e-7。
- 真空双球Lambda，PNJL热差通过半无限映射积分到无穷；外部频率为k0。
- 四通道逐个q1.4代表探针：原cut/独立径向/快速核、两种近轴围道、
  eta序列0.001/0.0003/0.0001及两个IR端点、内层加阶。
- 正常Mott动量处分割q积分；bulk每个节点独立根/cut检查。
  q采用8/16阶，不通过则32阶；总数和bound/pair/Landau分账均需收敛。
  q尾带8--12、12--16、16--24、24--32；尾部属于数值衰减证据，非解析全频界。
- 全频正谱、全复平面区间证明和所有可能拓扑不在已通过声明内。
  额外拓扑、Bose不安全或验收失败会拒绝合格ratio。
- 宽坐标求积使用临时BigFloat精度，因此使用独立进程而非共享线程。

目录只可新建或显式续算；续算核对配置、能量列表、源码和checkpoint哈希。
已结束的失败通道也保留，不因resume变成成功；需修复时开新版本目录。

## 产物

`ratios.csv`始终保留全部请求能量：`sqrt_s_NN_GeV`、`T_MeV`、`muB_MeV`、
`Kplus_over_pi_plus`、`Kminus_over_pi_minus`及独立passed标记。
失败ratio为NaN，图中断线，不插值跨失败点。`densities.csv`保留通道状态/原因。
密度单位fm^-3，ratio无量纲，内部q/omega为fm^-1。

每个energy子目录保存背景、通道结果、积分shell checkpoint。
JSON中的非有限值使用null，CSV中使用NaN；两者不代表零。
`run.json`固定配置/源码身份；`manifest.json`汇总最终状态及结果输出哈希，
图像目录另有 `plot_manifest.json` 记录图像与 `ratios.csv` 的输入哈希。
完整遍历但有失败点的状态为`complete_scan_with_failed_points`，不能冒充
`complete_research_curve_accepted`。PNG/PDF仅连通过点，线段用于引导视线。

公式及既有验收见
[无限热目标](../../../reference/formula/relaxtime/ChargedGBU_InfiniteThermal.md)。
