# Magnetic IMC 参数 provenance v1

## 结论

当前 Julia 磁场实现中的 `a=0.108805` 没有发现隐藏的单位补偿或额外十倍缩放。
`eB` 的生产入口以 `MeV^2` 接收，扫描入口只执行
`eB_fm2=eB_MeV2/hbarc^2`；`coupling_GB` 再以同一内部单位计算
`zeta=eB_fm2/(Lambda_QCD/hbarc)^2`。因此转换因子在 `zeta` 中抵消，不能解释
`0.108805` 与 `0.0108805` 的差异。

两篇 Ferreira 来源都明确给出 `a=0.0108805`：

| 来源 | 固定标识 | 原文位置 | `a` |
| --- | --- | --- | ---: |
| Ferreira, Costa, Lourenco, Frederico, Providencia | DOI [`10.1103/PhysRevD.89.116011`](https://doi.org/10.1103/PhysRevD.89.116011) | 式(12)参数段 | `0.0108805` |
| Ferreira, Costa, Providencia | DOI [`10.1103/PhysRevD.97.014014`](https://doi.org/10.1103/PhysRevD.97.014014) | 模型参数段（正文第2页） | `0.0108805` |

当前项目公式审核表引用的高雪艳博士论文表5-1给出 `a=0.108805`。本地证据文件为
`D:\Desktop\PNJL论文\同课题组毕业论文\高雪艳博士论文.pdf`，SHA-256 为
`421109FF06E8000D4BDD904EB27AAA1081F8692338480E1DB0DB92887086A904`。因此存在真实的
文献/版本分歧，而不是 Julia、Fortran 和 Python 之间的单位转换问题。

## 对当前实现的判断

- 如果生产目标是 Ferreira 的 IMC 参数化，当前默认值很可能是十倍抄录错误，属于
  高风险实现/配置问题；它会改变 `G(eB)` 的数值，不能通过放宽容差或重新积分来修复。
- 如果生产目标有意采用博士论文中的参数变体，当前数值不应直接删除，但必须把它命名为
  独立 model/profile，并在文档中同时列出 Ferreira 版本；不能继续把两个值描述成同一来源。
- 在作者确认 profile 语义前，本项目保留当前默认值，不生成新的 acceptance baseline，也不
  重跑 C0/C1/C2。

## 可复核证据

当前 Julia 调用链：

1. `src/models/scans/MagneticScan.jl` 将外部 `eB_MeV2` 转为 `eB_fm2`。
2. `src/models/pnjl_physics/core/MagneticThermodynamics.jl` 的
   `MagneticIMCParams` 和 `coupling_GB` 使用该内部值及 `Lambda_QCD_MeV=300`。
3. `config/models/pnjl/magnetic_default.toml` 与公式审核表重复声明 `a=0.108805`。

本轮下载的 Ferreira PDF 快照 hash 为：

- DOI `10.1103/PhysRevD.89.116011`：`3B1B2611FE4AA76EC8852FAF9E6CC7DDF687B752D2EF50FB7F36341F63EDC9AE`
- DOI `10.1103/PhysRevD.97.014014`：`DBB8773A711DB306BC52C0DD87005EFD2CE261520691AE29172AFE96BB5AA8EA`

外部快照的文件 hash、来源定位和 acceptance 状态见：
[`legacy_magnetic_external_crosscheck_v1.csv`](tables/legacy_magnetic_external_crosscheck_v1.csv)。
本文件不把外部 PDF 或外部仓库复制进主项目；只保留可追溯的 DOI、公式位置和结论边界。

## 物理影响的量级提示

在其他参数不变时，`a` 的十倍差异并非小扰动。使用当前 `b,c,d` 和
`zeta=eB/Lambda_QCD^2` 的公式，`eB=0.05 GeV^2`、`0.155 GeV^2`、`0.4 GeV^2`
处的 `G(eB)/G0` 会分别落在约 `1.0265/0.9965`、`1.2384/0.9663`、
`2.0768/0.7976`（前者为 `a=0.108805`，后者为 `a=0.0108805`）。这些值只是
参数敏感性诊断，不是生产结果或 acceptance target。
