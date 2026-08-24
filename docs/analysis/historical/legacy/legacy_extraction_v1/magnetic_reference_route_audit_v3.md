# 磁场参考路线审计 v3：贺伟博论文与 `pnjl_cep` 转向

## 状态

`blocked_pending_pnjl_cep_source`。本页完成了参考论文的公式审计，但尚未拿到
`pnjl_cep` 的仓库 URL、commit/SHA、生成脚本或机器可读输出，因此还不能生成
外部 acceptance target，也没有运行任何 equilibrium solver。

## 参考文献与可复核来源

- 文件：`D:\Desktop\PNJL论文\同课题组毕业论文\贺伟博 毕业论文.pdf`
- 文件 SHA-256：`4C2B3A51C8B3B8A3D976E2D12CE221A679E25336A1A0FD97BD36808980A524B5`
- PDF 总页数：57
- 磁场章节：第 4 章，印刷页 23--29，对应 PDF 页 33--39
- 相关小节：4.1“有磁场时的 PNJL 模型”、4.2“相关特殊函数的计算”、4.3 的
  `rho_B-mu_B`、Landau 能级和 Maxwell 等面积法则讨论

本页只记录论文明确写出的公式和由公式直接得到的路线判断；没有把论文中未说明
的数值节点、根选择或截断策略补成隐含合同。

## 论文明确给出的路线

| 论文位置 | 内容 | 对当前实现的含义 |
| --- | --- | --- |
| 式 (4-1) | 外部恒定磁场通过协变导数进入，并单独写出 Maxwell 项 | 当前物质巨势接口仍不包含 Maxwell 自能；这是既有边界，不是本轮修复目标 |
| 式 (4-2) | `Omega = sum_i (Omega_vac + Omega_med + Omega_mag + G phi_i^2) + 4 K phi_u phi_d phi_s + U` | 论文路线把磁场真空修正单独拆出，不是完整 Landau 真空项的一种记号改写 |
| 式 (4-3) | `Omega_vac` 是零场三动量截断真空项 | 真空基项使用零场 cutoff |
| 式 (4-4)-(4-5) | 热项按 Landau 能级求和并使用 PNJL `Z_i^+ + Z_i^-` | Landau 离散化保留在热项 |
| 式 (4-6)-(4-8) | `Omega_mag` 使用 Hurwitz-zeta 导数的 MFIR 磁场真空修正 | 这是论文给出的磁场真空正则化路线 |
| 式 (4-9)、表 4-1 | `G(eB)` 的 IMC 参数化 | 表中 `a=0.108805` 已由用户确认是笔误；正式 profile 采用 `a=0.0108805` |
| 式 (4-15)-(4-16) | 逐 Landau 能级的夸克密度与低温占据结构 | 低温 Landau onset、多个小 `S` 结构和高场反 `S` 结构可作为后续测试集的物理现象目标 |
| 4.3.4 | Maxwell 等面积法则 | 这是论文的相变后处理结论，不代表当前磁场 phase/Maxwell workflow 已实现 |

因此，论文路线是：

```text
zero-field vacuum cutoff
  + MFIR/Hurwitz-zeta magnetic vacuum correction
  + Landau-level PNJL thermal term
  -> five-field FixedMu stationarity
  -> density / phase-structure diagnostics
```

它不是当前 `MagneticIntegrals.omega0_flavor_landau` 所使用的“有限 `n_max` +
`p_z` Gauss 积分 + `smooth_cutoff(N=10)` 的完整 Landau 真空求和”。作者关于
“光滑截断有问题”的反馈与这一公式差异一致；因此不能再把 smooth-Landau
replay 当成磁场物理正确性的外部证明。

## 既有证据的重新分层

`magnetic_external_crosscheck_v2.md` 中的 Fortran/Julia 同路线比较仍可保留，但
它只证明在相同 profile、节点和 smooth-Landau kernel 下两个实现相互接近。它不证明
该正则化与论文或生产物理路线一致。Fortran 的 `quantity.f90` 派生量仍因味电荷
复用和 `NaN` 排除。

`pnjl_mag` 的 MFIR/Hurwitz-zeta 实现是有价值的独立路线诊断，但不是当前 Julia
smooth-Landau 路径的直接数值 target；除非主项目实现并冻结同一 MFIR 合同，否则不
把它提升为 acceptance。

## `pnjl_cep` 测试集转向

后续外部测试集优先采用作者指定的 `pnjl_cep`，流程如下：

1. 获取仓库 URL/路径、commit/SHA、运行依赖和输出文件或生成脚本；
2. solver-free 审计其磁场公式、真空正则化、`a` profile、单位、`n_max`/特殊函数
   约定、五维变量顺序、分支选择和 residual/status 字段；
3. 先抽取少量代表点，覆盖低温/高温、低/中/高 `eB`、Landau onset 与相结构变化，
   同时保留所有候选分支而不是只保留最低 `Omega`；
4. 只有在 source consistency、branch clarity、有限性和收敛门槛通过后，才把轻量
   fixed-point CSV 放入 `tests/validation/data/targets/`；原始输出只留在外部 provenance
   或分析目录；
5. 若 `pnjl_cep` 与论文 MFIR 路线不一致，建立独立 validation family，不把不同
   正则化的数值混成一个 target。

在第 1 步完成前，`pnjl_cep` 只能标记为 `pending_external_source`。不重跑 C0/C1/C2，
不扩展 Fortran 的 `n_max` 收敛扫描，也不修改当前 solver 或生产扫描路径。

## 主项目待办边界

- **必须后续处理**：把正式 profile `a=0.0108805` 接入唯一运行时配置，并修复
  `magnetic_default.toml` 未实际接入的问题；这与正则化路线确认分开验收。
- **路线决策后处理**：决定实现 MFIR/Hurwitz-zeta 路线，或明确保留 smooth-Landau
  仅作历史/诊断路线；在决定前不做性能或精度优化。
- **不在本轮**：Maxwell 自能、磁化率、各向异性压力、FixedRho/phase magnetic
  workflow、全域磁场扫描和 C0/C1/C2 重跑。
