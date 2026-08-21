# PNJL Magnetic 语义与稳定性审计任务单

## 1. 背景与审计边界

本任务从 `origin/main` 冻结基线 `6697feddd1c6382955727b93881b68569c165779` 的独立
worktree `D:\Temp\julia_relaxtime_core_algorithm_perf_ab` 建立。当前 `origin/main`
已前进到 `51c93cf0111415b35bb199376c358782c0f5a2f4`，该提交只增加 Maxwell CEP-local
diagnostic workflow，与 magnetic 文件无关；本审计继续以冻结 SHA 为代码证据，并在结论中
保留该分支差异。

目标是把 About 指向的 magnetic 路线公式、实现入口、脚本、测试和 fixed-point baseline
逐项对齐，先判断“当前输出代表什么物理量”，再决定是否值得做底层性能优化。

本任务只读：不修改 magnetic solver、不运行 `run_magnetic_*`、不生成新 baseline、不改变
默认精度/容差，也不触发 C0/C1/C2 或任何生产 workflow。若仓库材料不足以决定公式或
workflow 语义，状态保持为 `blocked_pending_author_confirmation`，不得用数值试跑代替物理确认。

## 2. 权威来源与路线

### About/能力清单

- `README.md` 将 magnetic 单点/扫描列为专题入口；完整能力和边界见
  `docs/reference/implemented_capabilities.md` 第 14 节。
- 该节给出的流程是：
  `eB + MagneticConfig -> Landau levels/smooth cutoff -> Omega_vac + Omega_T -> 5D gap solve -> pressure/densities + n_max report`。

### 公式来源

- `docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md`
  - Eq. (5-2)--(5-4)：磁场巨势、真空项和含 Polyakov 的热项；
  - Eq. (5-6)--(5-7)：Landau 能谱与平滑截断；
  - Eq. (5-9)：`G(B)` IMC 参数化；
  - Eq. (5-10)--(5-11)：完整 Polyakov 修正的密度；
  - Eq. (5-13)：低温近似密度；
  - 文档另列 Eq. (2-63) 的五个稳定平衡条件。
- 公式文档声称来源为高雪艳博士论文第 5.1 节；原文 PDF 已在本任务中逐页视觉核验，
  但论文自身存在符号/排版冲突，公式文档仍不是无条件的外部真值。

### 代码入口

- `src/models/Models.jl` 实际加载
  `src/models/pnjl_physics/core/MagneticIntegrals.jl` 与
  `src/models/pnjl_physics/core/MagneticThermodynamics.jl`。
- `src/models/pnjl_physics/PNJLMagneticModel.jl` 提供 `PNJLMagneticModel` 聚合接口。
- `scripts/pnjl/run_magnetic_point.jl`、`run_magnetic_eb_scan.jl`、
  `run_magnetic_stability_scan.jl` 和 baseline exporter 当前都使用固定的 5D `x_state`，
  没有在脚本内调用 `PNJLMagneticModel.solve_gap`。

### 新增原文证据（2026-08-19）

用户提供并核验了本地论文 `D:\Desktop\PNJL论文\同课题组毕业论文\高雪艳博士论文.pdf`，SHA-256 为
`421109FF06E8000D4BDD904EB27AAA1081F8692338480E1DB0DB92887086A904`。已使用原始 PDF 页面视觉核对而不只依赖文本抽取：印刷页 65--66 给出式(5-1)--(5-9)，印刷页 67--68 给出式(5-10)--(5-13)；印刷页 21--24 给出式(2-50)--(2-70)。

原文因此能够确认：磁场路线的完整 mean-field \(\Omega\) 由 Landau 真空项、未截断的 Polyakov 热项、手征项和对数 Polyakov 势组成；动力学质量使用式(2-55)并在磁场路线中采用 `G(eB)`；磁场计算固定 `mu_u=mu_d=mu_s`；一般温度密度应来自式(5-10)--(5-12)，式(5-13)只是 `Phi≈PhiBar≈0` 的低温近似；低温占据能级估计为 `Floor((mu_f^2-M_f^2)/(2|q_f|B))`。

原文不能无条件闭合所有实现合同。论文自身存在：式(2-63)重复 `partial Omega/partial phi_u` 的排版错误；式(5-10)导数符号与式(2-65)冲突；式(5-11)使用带符号 `q_f B`，而式(5-3)/(5-4)使用 `|q_f|eB`。因此密度符号/电荷绝对值、外部 Maxwell 项是否计入输出、固定态导数还是 on-shell 导数，以及数值 `n_max`/quadrature 参数，仍需项目物理合同确认。

原文核验后，M-1/M-2/M-4 的风险判断不降级：论文确认了“应对磁场 Omega 做五变量驻点求解”这一物理意图；实现已按用户合同建立磁场专用五维入口，仍需短 probe 检查数值分支覆盖。

## 3. 公式到实现对照

### 已可确认的一致部分

- `energy_landau` 实现 `sqrt(pz^2 + M^2 + 2n|q_f eB|)`，与 Eq. (5-6) 一致。
- `alpha_n(0)=1`、`alpha_n(n>0)=2`，与 Landau 简并度一致。
- `omega0_flavor_landau` 使用 `-N_c |q_f eB|/(2pi)`、Landau 求和和 `f_Lambda^2 E`；
  `smooth_cutoff` 的平方对应文档平滑截断形式。
- `omegat_flavor_landau` 使用 Eq. (5-4) 的 Polyakov 双对数结构。
- 非零 `eB` 的 `calculate_magnetic_omega_components` 将手征项、Polyakov 势、Landau 真空项
  和 Landau 热项组装为 `omega`，并在同一固定状态内使用 `G(B)` 生成质量和手征项。
- `|eB| <= 1e-14` 的零场退化路径与 API 文档声明的 `eB -> 0` 兼容合同一致；现有测试
  覆盖了零场 Omega 与普通 PNJL 热项的数值接近性。
- `magnetic_nmax_convergence_report` 确实比较 `n_base` 与 `n_base + delta_n`，没有把
  `resolve_nmax_from_cutoff` 单独当作收敛证明。

### 阻塞性不一致/未决项

#### M-1：gap solve 与 magnetic Omega 的驻点语义未闭合（已修复，待短 solver 验证）

`PNJLMagneticModel.solve_magnetic_gap` 现在在非零 `eB` 下以磁场 Omega 为目标，对
`(phi_u, phi_d, phi_s, Phi, PhiBar)` 做五维有限差分驻点 residual，并对多 seed 收集、
去重和按 Omega 选择。`solve_gap` 只返回选择后的 `MeanFieldState`；
`MagneticGapResult.candidates` 保留分支、方法、残差、物理性和可选 Hessian 稳定性标签。
零场仍显式退化到 PNJL 兼容路径。当前尚缺一次短的非零 `eB` solver stationarity probe，
所以本项实现已落地但不宣称数值分支覆盖完备。

#### M-2：`rho` 与 `number_densities` 使用了不同的分布语义（已修复）

- `omegat_flavor_landau` 使用含 `Phi/PhiBar` 的 Polyakov 对数；
- `calculate_magnetic_rho` 与 `density_flavor_landau` 现在共用含 `Phi/PhiBar` 的
  Polyakov 净密度 integrand；
- `calculate_magnetic_number_densities(...).quark` 保留历史字段名但正式表示 `q-qbar`，
  `baryon` 由其和除以 `3 rho0`。

固定点 smoke/nightly baseline 的 baryon 列已按该合同更新，unit test 也覆盖逐味 parity
和 `mu -> -mu` 奇宇称。

论文原页已经确认 Eq. (5-12) 的反夸克分布应使用 `E+mu`，且分子/分母的
`Phi/PhiBar` 交换结构与普通 PNJL 分布一致；公式文档已修正文本抽取造成的错位。
但 Eq. (5-10) 的导数符号和 Eq. (5-11) 的带符号 `q_f B` 仍与第 2 章热力学约定冲突；
本实现按用户确认的净密度和 `|q_f eB|` 相空间约定落地。

#### M-3：调用方的积分参数没有完整传入 magnetic Omega（已修复）

`omega_components`、`number_densities`、`model_rho/model_thermo` 现在显式转发
`p_num/t_num/pz_max/n_max/cutoff_N/xi`；零场退化也使用调用方的 `p_num/t_num`。模型
构造器可直接配置 magnetic `p_num/pz_max/n_max/cutoff_N`。`thermal_nodes` 仅作为 PNJL
兼容入参接收并忽略，因为 Landau 路线有独立的 pz 节点缓存。


#### M-4：`n_max` 估计公式的意图需要确认（已落实为显式数值合同）

公式文档的数值建议使用类似 `(μ^2-M^2)/(2|q|B)` 的占据能级估计，而实现
`resolve_nmax_from_cutoff` 使用 `max(Lambda^2, mu^2) + mass^2` 的 cutoff-based 估计。
当前实现继续把 `resolve_nmax_from_cutoff` 定义为数值求和起点，而不是低温物理占据上限；
调用方可显式传 `n_max`，并用 `magnetic_nmax_convergence_report` 对 `n_base` 与增量截断
比较。该边界已写入 API/公式文档，尚未把低温占据估计替换进真空项求和。

#### M-5：磁场文件存在结构性重复

`src/models/pnjl_physics/core/MagneticThermodynamics.jl`（SHA-256
`B138F18F478B507EB35907EE034E0EEDCC502188654D51377F09AD5F0A73031F`）与
`src/models/pnjl_magnetic/core/MagneticThermodynamics.jl`（SHA-256
`190EAAC73E2F72729F6559C5EF13F5A79179AFBCEAE572E00907DB53EE127A06`）并不相同；当前入口只加载
前者，且不再在初始化时向 `Main` 注入额外副本。后者是未接入当前入口的旧重复实现，未被证明造成运行时成本，但会增加维护和审计风险。
删除/迁移旧副本需要单独确认引用图，不在本语义审计中直接修改。

## 4. 现有测试与证据缺口

已有证据：

- `tests/unit/pnjl/test_magnetic_integrals.jl`：能谱、简并度、`n_max` 起点、平滑截断；
- `tests/unit/pnjl/test_magnetic_thermodynamics.jl`：有限性、`G(B)` 正性、零场 Omega parity；
- `tests/unit/pnjl/test_pnjl_magnetic_model.jl`：构造、委托接口、零场 Omega 可调用；
- `tests/regression/pnjl/test_magnetic_fixedpoint_regression.jl`：固定 `x_state` 的 smoke/nightly
  baseline、密度有限性和 `n_max` 报告。

已补充或仍缺：

- [x] 非零 `eB` 下五维 magnetic residual contract 与低节点短 solver stationarity/branch probe；
- [x] `calculate_magnetic_rho` 与 `calculate_magnetic_number_densities` 的逐味 parity；
- [x] `T>0`、`xi=0` 与一般 Polyakov 净密度边界测试；
- [x] `p_num/pz_max/n_max` 从模型调用方到 Omega/密度内核的参数传递测试；
- [x] `n_max` cutoff-based 起点与文档边界已固定；
- [x] magnetic workflow 正式要求 equilibrium solve，多分支候选由 branch-aware 结果保留。

## 5. 当前决定：语义实现完成，等待统一调用边界收尾

用户确认已解除物理合同阻塞；当前不建议：

- 运行新的生产 magnetic scan、C0/C1/C2 或默认精度变更；
- 在统一调用边界和分支稳定性验证完成前做 Landau pass 融合或磁场性能 A/B。

已完成的最小语义实现包括：非零 `eB` 五维 residual/多 seed branch-aware solver、
Polyakov 净密度统一、积分参数转发、`T>0`/`xi=0` 边界和 fixed-point baseline 迁移。

## 6. 作者确认后的实现合同

- `PNJLMagneticModel.solve_gap` 对非零 `eB` 必须求解磁场 Omega 的完整五维驻点；
  `solve_magnetic_gap` 暴露 branch-aware 候选，普通入口按 Omega 选择 equilibrium。
- `calculate_magnetic_number_densities` 使用一般温度的完整 Polyakov 净密度；Eq. (5-13)
  仅作为论文低温近似来源，不再作为默认实现分布。
- `rho_u/d/s` 与 `nd.baryon` 可以写入同一产物，但都必须来自同一净密度语义；固定点
  baseline 已迁移并保留 provenance。
- `resolve_nmax_from_cutoff` 是 cutoff-based 数值求和起点；显式 `n_max` 和增量收敛报告
  用于正式数值边界，不能把低温占据上限直接替代真空项截断。
- magnetic 路线只接受 `T_fm > 0`、`xi=0`；`p_num/pz_max/n_max/cutoff_N` 由模型配置或
  调用方显式传入，`thermal_nodes` 仅保留 PNJL 兼容入参。

## 7. 后续门槛

上述问题已得到明确回答并补充 unit/regression contract；低节点 stationarity/branch
probe 已完成。当前仍不进入磁场性能优化 PR，原因是共享 ProblemSpec 约束链不承载
磁场专用 Omega；该入口现在显式拒绝磁场约束模式，避免返回静默的 NaN 失败候选。

## 7A. 作者确认（2026-08-19）

用户已明确确认以下正式物理合同：

1. 磁场路线要求完整五维平衡态解：
   `∂Ω_B/∂(φu, φd, φs, Φ, Φbar) = 0`。
2. 磁场下可能存在多个稳定、亚稳和不稳定分支；实现必须保留多 seed 候选，不能只
   依赖单一 seed 或只暴露一个数值根。标准 equilibrium API 可以按明确的 branch policy
   返回一个 convenience state，但 stability/phase 路线必须能够访问候选集合及其分支/稳定性元数据。
3. 密度按 PNJL 家族约定以净密度为正式语义；`rho` 与 `number_densities` 不得使用
   互相矛盾的分布定义。保留既有 quark/antiquark 兼容字段时，应显式区分其与净密度
   字段的含义。
4. 不要求额外的 `T=0` 兼容路径；磁场五维 solver 可以明确要求有限且 `T>0`。
5. 当前不考虑联合 RS；磁场路线只允许 `xi=0`，非零 `xi` 应显式拒绝。
6. `n_max`、quadrature、参数传递等按审计建议落实为显式实现合同，并通过 focused
   测试固定边界行为。

上述确认解除 M-1/M-2/M-4 的作者确认阻塞，已进入最小语义闭合实现；当前仅允许一次
受控的短 solver stationarity/branch probe，仍不触发生产磁场扫描、C0/C1/C2 或性能 A/B。

### 短 stationarity/branch probe（2026-08-19）

- 参数：`T=150 MeV`、`mu_u=mu_d=mu_s=60 MeV`、`eB=20000 MeV^2`、`xi=0`；
  `p_num=4`、`t_num=4`、`pz_max=5`、显式 `n_max=1`；两枚显式 seed；关闭默认 seed
  catalog；`iterations=8`、`residual_norm_max=1e-3`；不计算 Hessian。
- 结果：`converged=true`，`attempt_count=2`，`failed_attempts=0`，去重后
  `candidate_count=1`，选中 `Omega=-1.3192642899314526`、`residual_norm=1.5701e-11`
  的候选，方法为 `trust_region`。
- 证据边界：证明非零磁场 Omega 的五维 finite-difference residual、NLsolve 调用、
  候选去重和 equilibrium 选择在低节点短样本可运行；不证明默认高节点精度、分支全集
  覆盖、Hessian 稳定性或生产扫描资格。

### Hessian stability probe（2026-08-19）

- 同一低节点、单 seed 参数开启 `classify_stability=true`；根 residual 为
  `4.1541e-11`，但有限差分 Hessian 分类为 `saddle_or_maximum`。
- 该结果保留在 `candidates`；Hessian 标签作为诊断信息，不提升为 PNJL 默认生产过滤条件，
  也不因没有 `local_minimum` 就否定已经收敛的驻点。普通入口继续按显式 branch policy
  选择 convenience state，stability/phase 路线保留全部候选；显式分类模式也不因
  Hessian 标签把 state 置空，调用方应自行应用 branch policy。
- 证据边界：仅验证稳定性分类标签可以工作；不证明 Hessian 在生产节点的物理完备性或全分支覆盖。

### 受限复核（2026-08-20）

- 在同一隔离 worktree、`T=0.7 fm^-1`、`mu=(0.4,0.4,0.4) fm^-1`、`eB=0.1 fm^-2`、
  `p_num=4`、`t_num=4`、`pz_max=5`、`n_max=1`、单显式 seed、`iterations=8` 下复核：
  未启用分类时得到 `residual_norm=1.1102e-11` 的单一候选，`Omega=-0.2472633958`，
  `trust_region` 4 次迭代；启用分类后同一根的 Hessian 特征值含约 `-4.9151`，
  标记为 `saddle_or_maximum`，但仍保留该候选和 convenience state。
- 该结果确认“驻点 residual 很小”不等于“有限差分 Hessian 局部最小”；结合 PNJL
  平均场的既有物理约定，Hessian 不提升为默认生产合同。标准 equilibrium 入口按明确的
  branch policy 选择 convenience state，stability/phase 路线必须访问并保留全部候选；
  Hessian 分类仅作诊断，不宣称已经完成生产分支全集或物理稳定性验收。
- `calculate_magnetic_number_densities` 增加显式 `net` 别名，保留历史 `quark` 字段，
  并将没有独立反夸克输出的事实写为 `antiquark=nothing`；这避免把净密度误接到普通
  PNJL 输运的独立 `quark/antiquark` 输入。
- 零场兼容复核：`PNJLMagneticModel.number_densities` 在 `eB≈0` 时委托普通
  `PNJLModel.number_densities`，仍返回独立 `quark/antiquark`；非零 `eB` 才使用上述
  magnetic 净密度合同。磁场 core 的 `calculate_magnetic_number_densities` 不作为普通
  PNJL 输运接口的替代品。

### capability 边界复核（2026-08-20）

- 通用 `number_densities` 合同明确要求独立 `quark`/`antiquark` 三味向量；非零
  `eB` 的磁场适配器只提供 `q-qbar` 净密度，因此不能把
  `supports_number_densities=true` 解读为普通输运兼容。
- `model_capabilities(::PNJLMagneticModel)` 现按配置动态报告：`eB≈0` 为 `true`，
  非零 `eB` 为 `false`；专用 `calculate_magnetic_number_densities(...).net` 仍是可用的
  磁场密度入口。该变化已同步 `implemented_capabilities.md`、magnetic API 文档和单测。
- 当前验证：磁场热力学 `26/26`、磁场模型 `27/27`、fixed-point smoke `30/30`、
  模型 API 同构 `21/21`、结构同构 `42/42`、反向依赖治理 `19/19`、文档一致性通过。
- 该 capability 修正不改变磁场 Omega、五维 solver 或既有 fixed-point baseline 的数值；
  仍未运行生产 magnetic scan、C0/C1/C2 或高节点全分支验证。
- `calculate_magnetic_pressure=-Omega` 的实现只给出固定外部磁场背景下的标量物质压力；
  当前没有磁化强度、横向/纵向压力拆分或 Maxwell 自能，不能直接宣称完整磁化介质 EOS。
- 公式到实现的逐项审核表已加入
  [`PNJL_magnetic_core.md`](../../reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md#开发者审核表公式到实现的对应关系)。
  表中同时区分当前 `main` 的固定态内核和未合并的五维候选分支，供后续 PR review 使用。

### 7B. 磁场 production adapter（当前隔离分支）

在作者确认的限定范围内，已新增专用 `Models.run_magnetic_scan`：

- 只接受 `model_kind=:PNJLMagnetic`、`solver_mode=:fixed_mu`、`T>0` 和 `xi=0`；
  外部输入为 `T_MeV`、共同 `mu_MeV` 和 `eB_MeV2`，内部显式转换到自然单位；
- 每个 `(T,mu,eB)` 点调用完整五维 `solve_magnetic_gap`，沿 `(T,eB)` 的 `mu` 顺序
  传递 continuation seed，同时保留多 seed 去重后的全部候选；
- selected CSV 写入 convenience state、热力学和净密度；candidates CSV 写入每个候选的
  seed、Omega、残差、方法、迭代数、分支/稳定性标签与 `n_max`；
- `run_unified_scan.jl scan magnetic` 转发该入口。普通 `TmuScan`、`TrhoScan` 和
  phase pipeline 对 `PNJLMagnetic` 显式报错，拒绝错误的 FixedRho/phase 语义；
- README、能力总账、API 细账和公式审核表已互相链接；旧 `run_magnetic_*` 脚本明确标为
  固定 `x_state` 诊断，不再描述为 equilibrium 扫描。

该实现阶段尚未声称默认高节点全域分支完备、外部 validation 或磁化介质 EOS；生产运行前
仍需代表点 `n_max/p_z/quadrature` 收敛审计和候选覆盖检查。

## 8. Definition of Done

- [x] About、公式来源、API、实现、脚本、测试和 baseline 已建立交叉索引；
- [x] 已确认无需运行新 magnetic solver 即可发现的高风险语义缺口；
- [x] 已明确停止条件和作者确认问题；
- [x] 作者确认 M-1/M-2/M-4 的正式物理语义；
- [x] 补充 parity/stationarity-contract/parameter propagation 测试；
- [x] 完成一次短的非零 `eB` solver stationarity/branch probe（低节点 diagnostic-only）；
- [x] 明确共享 `solve_constraint` 不接管 magnetic route，并对其返回显式边界错误；
- [x] 明确非零 `eB` 的通用独立密度 capability 边界，并同步能力清单/API 文档；
- [x] 建立论文公式、代码入口、主线现状、候选实现状态和证据边界的开发者审核表；
- [x] 在限定范围内接入 `(T,mu,eB)` 完整五维 FixedMu production adapter，并保留 selected/candidates 双输出；
- [x] 对普通 Tmu/Trho/phase magnetic 路径增加显式拒绝，并补充 adapter contract unit test；
- [x] 语义闭合后决定 magnetic 性能 A/B；已完成 A/B-M1 的低节点 diagnostic-only 部分测量，
  不改变生产 solver 语义；
- [ ] 任务完成后归档本文件。
