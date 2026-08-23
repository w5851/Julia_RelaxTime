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
   `a=0.0108805`，相差 10 倍；这可能是版本差异或旧代码笔误，不能未经确认
   当成相同模型。
2. Fortran 的 `BB quantity.f90` 在熵/密度后处理把 `E_d_n`、`E_s_n` 的
   横向动量构造误用 `q_u`；固定密度导出还包含 `NaN`，因此相关派生量不具备
   acceptance 资格。
3. Fortran residual 与 postprocess 使用的 Landau 截断并不一致，且目录中有
   多套 source/build 变体；`n=0..100`、128 点 `p_z in [0,20]` 只能作为历史
   配置，不能视为收敛证明。
4. `pnjl_mag` 是 Python/JAX/Optimistix 研究实现，PNJL 主模块采用 `fm` 单位，
   独立 NJL-MFIR 模块采用 `GeV` 单位；它没有 FixedRho 或全局极小值搜索。
   本审计只核对源码、README、commit 和已有测试声明，没有在本轮执行其平衡态
   扫描。

## 当前可复用的比较合同

未来若要把外部点升级为可信测试集，顺序应固定为：

1. 先冻结模型版本和参数，特别是 `G(eB)` 的 `a` 参数；若确认外部的
   `0.0108805` 是旧版，应为旧版建立显式版本标签，不能静默改当前 Julia。
2. 在 `(T,mu,eB)` 的同一单位和同一五维状态顺序下比较 residual，再比较
   `Omega`、压力和净密度。固定密度或等熵结果不能代替 FixedMu 点对点核对。
3. 对 `p_z` 节点、Landau `n_max`、真空/热积分节点分别做增量收敛；外部输出
   只有在自身收敛且与当前 Julia 的 `n_max` 语义一致时才可进入候选表。
4. 先做 Fortran 与 `pnjl_mag` 的 source-consistency，再与 Julia 比较。已知
   Fortran 后处理缺陷覆盖到的 entropy/density 行不得作为 gate 输入。
5. 只有参数、单位、截断、分支和 observable 定义全部一致，且多个代表点通过
   预先声明的容差，才可以新增
   `tests/validation/data/targets/pnjl/reference/` 下的轻量 acceptance CSV
   和对应测试。当前没有新增 acceptance target。

## 来源快照

| 来源 | 固定标识 | 关键实现 | 当前状态 |
| --- | --- | --- | --- |
| 旧 Fortran | `external://legacy/fortran-magnetic-field-validation`, tree SHA-256 `a51a488d...` | `BB equation.f90`、`BB quantity.f90`、`BB setting constants.f90` | `diagnostic_only`；保留公式/缺陷证据 |
| `pnjl_mag` | Git commit `e1fc81d3c3c9d220c49972e54307b66a604cb9db` | `src/pnjl_mag.py`、`src/gap.py`、`src/magnetic_jet.py` | `candidate_not_aligned`；尚未执行独立 replay |
| 当前 Julia | project commit `1ccf29310fb20c30bcd154f0b4966e25a7565225` | `PNJLMagneticModel` 与 magnetic fixed-point baselines | 当前内部合同，不是外部证据 |

完整文件 SHA、单位、截断、已知问题和 gate 状态以 CSV 为准；这里的短 SHA 仅用于阅读。

## 对 legacy 清理的影响

本记录不要求把 Fortran 源码、二进制或 1.5 亿字节级扫描复制进仓库。清理
`D:\Desktop\legacy` 前必须保留本表、提取包 manifest 及当前 Mott provenance；
若后续真的执行匹配 replay，应在删除前完成并把新结果以独立 hash-bound artifact
登记。否则删除后只能保留本轮的诊断性结论，不能声称外部 acceptance 已完成。
