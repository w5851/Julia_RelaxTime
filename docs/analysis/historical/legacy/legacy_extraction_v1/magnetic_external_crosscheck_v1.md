# 磁场外部数值交叉核对 v1

## 目的和结论边界

本记录把旧 Fortran 磁场路线和
[`ZhouRui-xzit/pnjl_mag`](https://github.com/ZhouRui-xzit/pnjl_mag) 作为
独立参考源登记下来。它们目前都是 `diagnostic_only` 或
`candidate_not_aligned`，还不是默认 validation acceptance。机器可读的逐源
登记见
[`legacy_magnetic_external_crosscheck_v1.csv`](tables/legacy_magnetic_external_crosscheck_v1.csv)。

当前不能把两个外部项目直接拼成“可信测试集”，原因不是缺少输出，而是比较
合同还没有闭合：

1. 当前 Julia 磁场参数使用 `a=0.108805`，旧 Fortran 和 `pnjl_mag` 使用
   `a=0.0108805`。这不是单位转换造成的十倍补偿：Julia 只把外部 `eB` 从
   `MeV^2` 转为内部 `fm^-2`，随后在 `zeta=eB/Lambda_QCD^2` 中使用同一单位，
   分子分母的单位转换会抵消。Ferreira 2014（DOI `10.1103/PhysRevD.89.116011`）
   的式(12)和 Ferreira 2018（DOI `10.1103/PhysRevD.97.014014`）正文都明确给出
   `a=0.0108805`；高雪艳博士论文表5-1和当前 Julia 给出 `0.108805`。因此，若当前
   目标是 Ferreira 参数化，`0.108805` 很可能是十倍抄录/版本选择错误；若项目有意采用
   学位论文变体，则它是一个需要显式命名的模型版本，不能静默当作同一模型。该分歧会
   改变 `G(eB)` 的数量级和 IMC 方向，不能用数值容差吸收。
2. Fortran 的 `BB equation.f90` 实际是六维固定密度 residual：
   `(phi_u, phi_d, phi_s, Phi, PhiBar, mu_B)`，前五个分量才是给定 `T,rho_B,eB` 下的
   五维平衡态字段。`BB quantity.f90` 在熵/密度后处理把 `E_d_n`、`E_s_n` 的横向动量
   构造误用 `q_u`；固定密度导出还包含 `NaN`，因此这些派生量不具备 acceptance 资格。
   修正后的 residual 根和五维状态仍可作为候选外部参考，但必须保留固定密度合同标签。
3. Fortran residual 与 postprocess 使用的 Landau 截断并不一致，且目录中有多套 source/build
   变体；`n=0..100`、128 点 `p_z in [0,20]` 只能作为历史配置，不能视为收敛证明。用户已
   允许在临时诊断副本中统一截断以控制变量；这不会把旧截断自动提升为项目标准。
4. `pnjl_mag` 是 Python/JAX/Optimistix 研究实现，PNJL 主模块采用 `fm` 单位，独立
   NJL-MFIR 模块采用 `GeV` 单位；它提供五维 FixedMu 局部驻点，但没有 FixedRho 或全局
   极小值搜索。其 PNJL 真空项采用 MFIR/Hurwitz-zeta 路线，不能仅靠同步 `n_max` 就与
   Fortran 平滑截断路线视为同一离散化。本审计没有执行完整平衡态扫描。

## 当前可复用的比较合同

未来若要把外部点升级为可信测试集，顺序应固定为：

1. 先冻结模型版本和参数，特别是 `G(eB)` 的 `a` 参数。当前 Julia 的
   `0.108805` 与 Ferreira/外部实现的 `0.0108805` 必须先由作者决定是修正默认值，
   还是显式保留为独立参数 profile；本轮不擅自修改默认值。
2. 为 `(T,mu_B,eB)` 建立 FixedMu 对照点。`pnjl_mag` 可直接提供五维局部驻点；Fortran
   需要从其六维 FixedRho 主程序中抽取一个固定 `mu_B` 的临时诊断入口，或明确把比较
   标记为不同 ensemble。两者都先比较五维状态和 stationarity residual，再比较 `Omega`。
3. 对 `p_z` 节点、Landau `n_max`、真空/热积分节点分别做增量收敛；用户允许在 Fortran
   临时副本中把截断改到与对照实现一致，但只有收敛证据才可进入候选表。
4. Fortran 与 `pnjl_mag` 的真空正则化并不相同，不能在没有 MFIR/平滑截断转换说明时把
   两者的 `Omega` 差异解释成 Julia 错误。Fortran 后处理缺陷覆盖到的 entropy/density
   行永不作为 gate 输入；可用范围限于修正后的 residual、五维状态及明确口径的 `Omega`。
5. 只有参数、单位、ensemble、截断、分支和 observable 定义全部一致，且多个代表点通过
   预先声明的容差，才可以新增
   `tests/validation/data/targets/pnjl/reference/` 下的轻量 acceptance CSV
   和对应测试。当前没有新增 acceptance target。

## 来源快照

| 来源 | 固定标识 | 关键实现 | 当前状态 |
| --- | --- | --- | --- |
| 旧 Fortran | `external://legacy/fortran-magnetic-field-validation`, tree SHA-256 `a51a488d...` | `BB equation.f90`、`BB quantity.f90`、`BB setting constants.f90` | `diagnostic_only`；保留公式/缺陷证据，六维 FixedRho 合同 |
| `pnjl_mag` | Git commit `e1fc81d3c3c9d220c49972e54307b66a604cb9db` | `src/pnjl_mag.py`、`src/gap.py`、`src/magnetic_jet.py` | `candidate_not_aligned`；尚未执行独立 replay |
| Ferreira 2014 | DOI `10.1103/PhysRevD.89.116011`，式(12) | `G_s(eB)` IMC 参数化 | `parameter_source`；明确给出 `a=0.0108805` |
| Ferreira 2018 | DOI `10.1103/PhysRevD.97.014014`，式(3)附近参数段 | 磁场 PNJL 参数化 | `parameter_source`；明确给出 `a=0.0108805` |
| 高雪艳博士论文 | 外部论文第5章表5-1 | 当前公式审核表/配置来源 | `project_source_variant`；给出 `a=0.108805`，与 Ferreira 来源冲突 |
| 当前 Julia | project commit `1ccf29310fb20c30bcd154f0b4966e25a7565225` | `PNJLMagneticModel` 与 magnetic fixed-point baselines | 当前内部合同，不是外部证据 |

完整文件 SHA、单位、截断、已知问题和 gate 状态以 CSV 为准；这里的短 SHA 仅用于阅读。

## 对 legacy 清理的影响

本记录不要求把 Fortran 源码、二进制或 1.5 亿字节级扫描复制进仓库。清理
`D:\Desktop\legacy` 前必须保留本表、提取包 manifest 及当前 Mott provenance；
若后续真的执行匹配 replay，应在删除前完成并把新结果以独立 hash-bound artifact
登记。否则删除后只能保留本轮的诊断性结论，不能声称外部 acceptance 已完成。
