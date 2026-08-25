# 磁场参考路线审计 v3：贺伟博论文与 `pnjl_mag` 外部参考

## 状态

`pnjl_mag_cross_solver_diagnostic_passed_pending_convergence_target_gate`。本页已完成参考论文
和 `pnjl_mag` 源码的静态 source-gate，并已固定仓库、commit、关键文件 hash 和
合同差异；当前已在本机重放外部单 seed equilibrium continuation，并保留 9 个代表点，
同时仅接纳一个 fixed-state `Omega` kernel-only oracle。外部 equilibrium 重放仍为
`diagnostic_only`；已形成小规模 Julia/外部跨求解器对照，但尚未形成 acceptance target。

机器可读 gate 表：
[`tables/pnjl_mag_source_gate_v1.csv`](tables/pnjl_mag_source_gate_v1.csv)。
结论边界见
[`tables/pnjl_mag_source_gate_claim_ledger_v1.csv`](tables/pnjl_mag_source_gate_claim_ledger_v1.csv)。
外部仓库已有数据文件的范围和 hash 见
[`tables/pnjl_mag_external_data_inventory_v1.csv`](tables/pnjl_mag_external_data_inventory_v1.csv)。

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

`pnjl_mag` 的 MFIR/Hurwitz-zeta 实现与当前 Julia 默认 MFIR 路线属于同一公式家族，
因此是首选外部 source candidate。但历史 replay 使用的求解设置和输出合同尚未完成
同口径审计，不能把既有 replay 行直接提升为 acceptance。

## `pnjl_mag` 测试集接入

外部源码已经固定：

- remote：`https://github.com/ZhouRui-xzit/pnjl_mag.git`
- commit：`e1fc81d3c3c9d220c49972e54307b66a604cb9db`
- 本机只读审计副本：`D:\Temp\pnjl_mag_audit_20260823`
- 已确认的高层结构：五维 FixedMu、JAX 自动微分、Optimistix 求根、
  MFIR/Hurwitz-zeta 真空项和 Landau 热项

静态 source-gate 结论如下：

- **已通过**：MFIR/Hurwitz-zeta 公式家族、五维变量顺序、声明的 PNJL/IMC 参数。
- **需 adapter**：外部 `muB` 是 baryon chemical potential，源码内部除以 3；Julia
  入口接受共同的 flavor `mu_vec`。
- **条件通过**：两边都使用自然单位，但外部 `hc=197.33` 与 Julia
  `hbarc=197.3269804` 不完全相同，接纳前必须固定换算常数。
- **不一致**：外部允许 `eB=0` 并对负值取绝对值；Julia 生产契约只接受
  `eB >= 100 MeV^2` 的正标量。
- **已完成诊断、未完成接纳**：9 点筛查和 3 个高温匹配节点已完成 Julia 多 seed
  对照；匹配节点下 Omega 差最大约 `1.18e-12`。`T=240 MeV,eB=0.8 GeV^2`
  发现第二个物理候选，说明 branch-aware 输出确有必要。外部实现仍是一 seed
  continuation，不是 Julia 的多 seed candidate 集合。
- **待完成**：独立 `p_z`、Landau level、zeta 节点收敛矩阵、正式输出字段映射和
  branch-policy/target admission；筛查节点下的差异不能替代匹配节点收敛证据。
- **已有数据但不接纳**：`orders_muB0.csv` 和 `reduced_magnetic_susceptibilities_muB0.csv`
  的行数、范围、hash 和生成口径已登记；它们只有 `muB=0`，输出没有完整 solver
  status/候选分支/独立收敛矩阵，因此目前只能作为诊断输入清单。
- **本机 equilibrium 重放通过但不接纳**：外部 `.venv` 已通过锁定的 `uv.lock` 建立，
  作者测试 `9/9` 通过；按 `T=300..50 MeV`、步长 `1 MeV` 的下降 continuation，在
  `eB={0.2,0.4,0.8} GeV^2` 提取 `T={50,150,240} MeV` 共 9 点。五维状态与作者提交的
  `orders_muB0.csv` 逐字段一致，最大 gap residual 为 `1.65e-13`。证据见
  [`pnjl_mag_equilibrium_replay_v1/`](pnjl_mag_equilibrium_replay_v1/)。该路线仍只有一个
  continuation 分支，不能证明 Julia 多 seed candidate 集合等价。
- **跨求解器诊断通过但不接纳**：[`pnjl_mag_cross_solver_replay_v1/`](pnjl_mag_cross_solver_replay_v1/)
  在全部 9 点做筛查，并对 `T=240 MeV` 的 3 点使用匹配节点；匹配节点下外部状态
  residual 最大约 `3.1e-12`、Omega 差最大约 `1.2e-12`。`eB=0.8 GeV^2` 的高温点
  找到两个物理候选，外部 continuation 分支为较低 Omega 候选。该证据尚不包含
  全分支完备性或独立收敛矩阵。
- **已接纳的最小外部 oracle**：`tests/validation/data/targets/pnjl/reference/`
  中的单点 `Omega` 只验证 MFIR kernel；测试显式记录外部 `hc=197.33` 的单位映射，
  不调用 root solver，也不约束 equilibrium branch 或扫描。

后续接入流程如下：

1. 已固定可重放的依赖、生成脚本、参数清单，并建立外部单分支诊断输出 schema；
2. 已对代表点执行 Julia 多 seed candidate 重放，显式记录 continuation 分支和额外
   候选；
3. 下一步完成 source consistency、branch clarity、有限性和收敛门槛后，才把轻量
   fixed-point CSV 放入 `tests/validation/data/targets/`；原始输出只留在外部 provenance
   或分析目录；
4. 若 `pnjl_mag` 的具体公式或参数合同与当前 Julia MFIR 路线不一致，建立独立
   validation family，不把不同
   正则化的数值混成一个 target。

源码和外部 equilibrium 可重放性、有限范围的跨求解器 adapter 已不再阻塞；当前阻塞点是
numeric convergence、正式 output schema、branch policy 和 target admission。不重跑
C0/C1/C2，不扩展 Fortran 的 `n_max` 收敛扫描，也不把历史 replay 直接复制成 target。

## 主项目待办边界

- **已完成的主项目路线决策**：正式 profile `a=0.0108805` 已接入
  `magnetic_default.toml`；MFIR/Hurwitz-zeta 为默认生产真空路线，完整
  smooth-Landau 只通过 `route=:landau_legacy` 保留为历史/诊断路线。
- **必须后续处理**：在不改变生产路线的前提下，补充多节点收敛矩阵、固定 branch
  policy 和正式输出映射，再决定是否生成轻量 equilibrium acceptance targets。
- **不在本轮**：Maxwell 自能、磁化率、各向异性压力、FixedRho/phase magnetic
  workflow、全域磁场扫描和 C0/C1/C2 重跑。
