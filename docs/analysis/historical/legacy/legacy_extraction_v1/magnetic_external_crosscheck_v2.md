# 磁场外部 FixedMu 交叉验证 v2

## 范围

本记录是在 v1 provenance 基础上的短、定向 replay。它只覆盖五维
`(phi_u, phi_d, phi_s, Phi, PhiBar)` FixedMu 局部驻点，不触发磁场生产扫描，
也不改变 `D:\Desktop\legacy` 原始源文件。

本轮新增的机器可读证据是
[`magnetic_external_fixedmu_replay_v2.csv`](tables/magnetic_external_fixedmu_replay_v2.csv)。
Julia/Fortran 同路线状态差异见
[`magnetic_julia_fortran_same_route_v2.csv`](tables/magnetic_julia_fortran_same_route_v2.csv)。

## 1. IMC 参数结论

`a=0.108805` 与 Fortran/`pnjl_mag` 的 `a=0.0108805` 不是单位转换造成的差异。
Julia 入口把 `eB_MeV2` 转为 `eB_fm2`，而 `zeta=eB/Lambda_QCD^2` 的分子和分母
使用同一内部单位，转换因子抵消。高雪艳博士论文表 5-1 明确写出
`0.108805`，Ferreira 2014/2018、旧 Fortran 和 `pnjl_mag` 明确写出
`0.0108805`，因此这是两个 model profile 的真实冲突。

源码审计还发现：`config/models/pnjl/magnetic_default.toml` 声明了 `a=0.108805`，
但 `PNJLMagneticModel` 当前直接调用硬编码的 `default_imc_params()`；没有看到该
TOML 被运行时构造器读取。这是配置权威性缺口，不是十倍差异的隐藏补偿。

因此本轮不擅自切换生产默认值，也不重写 C0/C1/C2。若项目目标是 Ferreira profile，
当前默认值应视为高风险实现/配置问题；若目标是论文 profile，则必须把它命名为
独立 profile，不能继续把两者当作同一模型。

## 2. Fortran 截断探针

在 `D:\Temp\pnjl_fortran_fixedmu_probe` 临时副本中，把 `n_max` 和 `p_z` 上限变为
运行参数，并将 Landau 数组扩展到 `0:400`。原始 `D:\Desktop\legacy` 未修改。
使用 `p_z` Gauss 128 节点、平滑 cutoff `N=10`、`a=0.0108805`。

主要观察：

- `n_max=100, p_z=20 -> 25 fm^-1` 对点 1/2/4/5 的状态变化很小，说明
  `p_z=20` 不是本组点的主误差源；
- `n_max=100 -> 300/400` 在点 1/2 的状态仍有约 `1e-4 fm^-1` 量级漂移，旧的
  `n=0..100` 不能当作收敛证明；
- 点 3 `(250,0,0.04)` 的 residual norm 约 `3.7--3.9`，排除；
- 点 1、2、4、5 在 `n_max=400,p_z=25` 下 residual norm 分别约为
  `4.14e-8, 7.88e-7, 4.08e-7, 1.26e-8`，可作为候选状态证据；
- 老程序的 `check=.false.` 不能单独视为收敛，所有资格判断都以显式 residual gate
  为准。`BB quantity.f90` 的熵/密度派生量和含 `NaN` 的固定密度导出仍然排除。

## 3. 同路线 Julia 对照

为隔离 IMC profile 影响，临时 Julia probe 显式传入
`MagneticIMCParams(a=0.0108805)`，并使用 `n_max=300,p_num=128,p_z=25`。
点 1、2 的 Julia residual 分别为 `4.16e-9`、`4.41e-9`；与 Fortran
`n_max=400` 状态的最大绝对分量差约为 `1.1e-4 fm^-1`。这支持如下较强但有边界的
结论：当前 Julia 的平滑 Landau residual/state kernel，在 profile、截断和积分节点对齐后，
与 Fortran 路线数值一致到本轮探针的误差范围内。Fortran 临时 adapter 没有输出完整
`Omega`，因此本轮不宣称两者的 `Omega` 已完成对照；同路线表中的 `omega_fm4` 对
Fortran 行保持为空。

这不是当前生产默认 profile 的 acceptance：对照使用了显式 override
`a=0.0108805`，而当前默认仍为 `0.108805`。

## 4. `pnjl_mag` 的位置

`pnjl_mag@e1fc81d3c3c9d220c49972e54307b66a604cb9db` 使用五维 FixedMu、JAX/Optimistix
和 MFIR/Hurwitz-zeta 真空项。低密度点 1/2 得到 residual 约 `1e-15` 的局部根；
其状态与 Fortran/Julia 平滑 Landau 根明显不同。点 4/5 在喂入 Fortran 状态后仍未找到
可接受根；`throw=false` 得到的部分状态还越出 `Phi, PhiBar in [0,1]` 的物理范围。

这不是“某一方已经被证明错误”：两条路线的真空正则化、Landau 离散化和 branch policy
不同。仅把 `n_max` 改成相同值不能消除 MFIR 与平滑 Landau 的模型差异。

因此当前 admission 分层是：

| source family | 可复用范围 | 当前状态 |
| --- | --- | --- |
| Fortran smooth-Landau | 五维状态、显式 residual；与 Julia 同路线对照 | `candidate_external_reference`，待 profile 冻结后才能晋升 acceptance |
| `pnjl_mag` MFIR | 独立公式/局部根/分支诊断 | `diagnostic_only`，不能直接作为当前平滑 Landau 路线的数值 target |
| Fortran `quantity.f90` 派生量 | 熵、密度、固定密度导出 | `excluded`，存在味电荷复用、截断不一致和 `NaN` |

## 5. 对主项目的影响

本轮确实暴露了两个实现层问题：

1. IMC profile 的论文来源与外部实现不一致，且运行时硬编码没有兑现磁场 TOML 的
   配置权威性；
2. 磁场外部验证缺少“profile、ensemble、正则化、截断、branch、residual”六元合同。

本轮没有证明需要修改 Landau 能谱、净密度定义或求解器算法；也没有触发完整磁场扫描。

## 下一步与 acceptance gate

1. 作者先决定生产 profile：Ferreira `a=0.0108805`，或显式保留论文 profile
   `a=0.108805`；随后修复/明确 TOML 与硬编码的权威关系。
2. 若采用 Ferreira profile，把 Fortran 点 1/2（必要时 4/5）作为小型外部 state/residual
   validation candidate；在同一 profile 下更新 Julia expected values，并保留 source/tree
   hash 和截断设置。
3. `pnjl_mag` 只有在主项目增加同一 MFIR 路线，或完成可审计的正则化转换后，才能进入
   独立 acceptance family；在此之前保持诊断证据。
4. 点 3、高 `mu_B` 的 `pnjl_mag` 失败点和所有 Fortran entropy/density 派生量不进入
   acceptance CSV。

本轮没有必要重跑 C0/C1/C2。

## Replay provenance

- Julia project: `codex/legacy-extraction@b9a2441e`
- Fortran temporary adapter: `fixedmu_probe.f90`, SHA-256
  `3015CC8196432E3973017A15E29CE7C7AA019644BD9179DD29D664F72B8541B9`
- Fortran temporary constants: `BB constants.f90`, SHA-256
  `9C5BA625EC3F8053ACF8DB97B4F02E6EF3D35BC3D91776B601E89005C47556E0`
- Fortran temporary quadrature setup: `BB setting constants.f90`, SHA-256
  `86A09341CBE6CCEDAD65596E41CBD8A385A0B420BE6D54EFB37DE3B4808E5D30`
- `pnjl_mag` strict replay script, SHA-256
  `8F1426A61C4CE62F9E172D18B23256D9756B44ACED6E9A3A2695E3F3A3AB9CBE`
- `pnjl_mag` Fortran-seed replay script, SHA-256
  `ED40186BA1E6882F129B91EFEC7C0C3C6AF2979DB6415F3DE0A0B5EF10217069`
- `pnjl_mag` continuation diagnostic script, SHA-256
  `0CB828D3949468203B687F219DBFDFB56644EBA0C4E08410DC70015E1D35FB64`
- Julia same-route comparison script, SHA-256
  `7276389C0E437EACEFFBA39440FA95CEDC1D32063C9F7E824D786BD9FFA44D29`
