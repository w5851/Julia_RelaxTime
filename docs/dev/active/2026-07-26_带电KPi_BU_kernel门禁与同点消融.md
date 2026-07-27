# 带电 K/π BU kernel 门禁与同点消融任务单

更新日期：2026-07-26

当前状态：completed（verdict：`diagnostic-only`），承接 `docs/dev/backlog/2026-07-25_带电KPi介子数密度物理约束与验证中期路线.md` 的 Phase 0、Phase 1 和 Phase 2 首批工作。

## 1. 背景与目标

现有 `trho_asymmetric_kplus_piplus_scan_v1` 在既定数值策略下已经通过 production convergence gate，但部分 BU/GBU `K^+/π^+` 结果远高于实验常见的 `O(10^-1)`。在修改平衡态约束前，必须先回答两个更窄的问题：

1. current/GBU 导数形式与分部积分形式在同一个有限 `ω` 窗口内，加入上下边界项后能否闭合；
2. 在固定外部点和固定 kernel 下，把 `μ_s/μ_q` 从 `0.2` 改为 `0.55` 时，`K^+/π^+`、`K^-/π^-` 及组成密度如何变化。

本任务的交付不是修改正式 BU 数值语义，而是建立一个能决定后续是否允许进入 `μ_Q、μ_S` solver 开发的 analysis gate。

## 2. 本轮 harness

### 2.1 允许修改

- `scripts/analysis/relaxtime/` 下的 BU kernel gate 与文献对齐诊断脚本；
- 支撑上述脚本的窄 analysis helper；
- 对应 unit/integration smoke；
- `scripts/analysis/README.md`、本任务单和必要的分析说明；
- 由本任务命令生成的独立 diagnostic output。

### 2.2 不要修改

- `phase_shift_meson_number_density` 的正式默认公式、相移 convention、display policy 或 Bose-domain policy；
- `FixedAsymmetricRho`、`μ_Q/μ_S` solver、冻结线默认值和正式 production artifact；
- mixed-meson governance、non-fixedmu unified joint-solve 语义；
- transport 数值口径和现有 regression baseline。

### 2.3 必须保持

- 内部量继续使用自然单位，MeV/fm^-1 字段显式标注；
- `K/π=n_K/n_π` 方向不变；
- `K^+/π^+` 与 `K^-/π^-` 必须来自同一外部点、同一 profile 对应的同一个 upstream state；
- `μ_π^neq` 与 flavor difference 产生的 signed `μ_K` 分栏输出；
- analysis 结论不得升级为正式 kernel 修复或实验拟合结论。

## 3. 当前能力与缺口

- [x] `scripts/analysis/relaxtime/compare_bu_derivative_vs_byparts_e5.jl` 已在同一 GL 网格比较 current/GBU 的 derivative 与 by-parts bulk。
- [x] `scripts/analysis/relaxtime/check_bu_levinson_e5.jl` 已输出低能端、高能端、相位跨度及最近 `π` 倍数。
- [x] `scripts/analysis/relaxtime/audit_bu_meson_density_literature_alignment.jl` 已覆盖 `μ_s/μ_q=0.2` 下 `μ_K=0` 与 signed `μ_K` 的正负通道诊断。
- [x] derivative-vs-byparts 诊断已加入有限窗口边界项 `[gF]_{ω_min}^{ω_max}`，旧相对差保留为 bulk-only 历史映射。
- [x] charged audit 已用同一字段合同覆盖 `μ_s/μ_q=0.2/0.55`，并输出 flavor chemical potentials、profile、Bose status 和 `m_M-μ_M`。
- [x] 已新增不依赖 PNJL 重求解的纯数值 unit gate。

## 4. 任务分解

### A. 有限窗口恒等式 helper

- [x] 实现纯数值 helper，输入同一 `ω` 节点/权重上的 `F(δ)` 与 `dF/dω`，输出 derivative integral、by-parts bulk、上下边界项、by-parts total 和 closure residual。
- [x] 对 `T<=0`、`ω_max<=ω_min`、`ω_min<=μ`、长度不一致和非有限输入给出明确 `ArgumentError`。
- [x] unit 覆盖 constant/linear phase weight、current/GBU 两种 `F`、边界项符号和非均匀三点导数。

### B. 现有 E5 对比脚本接入边界项

- [x] 在每个 `(T, meson, q)` 上把 `ω_min`、GL 内点和 `ω_max` 作为同一序列 unwrap，避免端点与内点使用不同相移分支。
- [x] 同时输出 current/GBU 的 bulk、boundary、total、derivative 和 closure residual。
- [x] 保留旧 density 列的可解释映射，避免历史 CSV 无法对读。
- [x] 输出端点 weighted phase、Bose occupation 和累计边界贡献，支持判断低能端还是高能端主导。

### C. `μ_s=0.2/0.55` 双通道同点消融

- [x] 从 `FlavorChemicalProfiles` 配置读取 `bu2020_mu_s_0p2` 与 `friesen2019_mu_s_0p55`，不在脚本中复制 profile 系数。
- [x] 每个 profile 只求一次 upstream state，再在其上计算 `K^+/π^+`、`K^-/π^-` 的 `μ_K=0` 与 signed `μ_K`。
- [x] 输出 `flavor_profile`、`μ_u/d/s`、`μ_π`、`μ_K`、四个组成密度、ratio、channel status、`m_M-μ_M` 和数值窗口。
- [x] 明确 `0.2/0.55` 是 sensitivity/literature scenario，不是守恒荷平衡解。

### D. 测试与诊断产物

- [x] unit helper 测试通过。
- [x] analysis script contract/smoke 覆盖新增列与两个 flavor profiles、两个 charge channels。
- [x] 运行一个低成本代表点，确认输出有限性、状态和 closure 字段可解释。
- [x] 根据代表点结果给出 `diagnostic-only` gate；判据来自 closure/Bose 域，不来自 ratio 是否接近 `0.3`。

### E. 文档与治理

- [x] 更新 `scripts/analysis/README.md` 的脚本说明和输出口径。
- [x] 本轮只新增 analysis helper，未新增稳定公共入口，因此无需更新 `docs/api/`。
- [x] 运行 `git diff --check`、docs consistency、脚本治理和与改动风险匹配的聚焦测试。

## 5. 测试与验收标准

### Unit

- 有限窗口恒等式对解析平滑样例随节点加密收敛；
- boundary 使用 `[g(ω)F(ω)]_{ω_min}^{ω_max}/(2π)` 的符号；
- `by_parts_total = by_parts_bulk + boundary`；
- 非法 Bose 支撑域不静默用 `nextfloat(μ)` 继续。

### Integration / analysis smoke

- E5 输出同时包含 current/GBU 的 boundary 与 closure；
- charged 输出恰好覆盖两个 flavor profiles、两个 charge channels 和两个 `μ_K` 规则；
- 同一 `equilibrium_group` 内四种后处理共享相同 flavor chemical potentials 和 upstream state identifier。

### Regression / validation

- 本任务不更新正式数值 baseline，因为正式 kernel 和 workflow 默认语义不变；
- 既有 meson-density unit/integration smoke 必须继续通过；
- 文献和实验只用于机制、趋势与数量级解释，不作为本任务的逐点 pass/fail target。

## 6. 里程碑

1. [x] M1：纯数值 finite-window identity gate 与 unit 完成。
2. [x] M2：E5 boundary-aware CSV 和 `0.2/0.55` charged ablation CSV 完成。
3. [x] M3：代表点 verdict 为 `diagnostic-only`，允许拆出后续 solver 任务，但不允许修改正式 kernel。

只有 M3 verdict 不是 `blocked`，才拆出固定 `(T,μ_B)` 下 `μ_Q、μ_S` 联合求解的下一张 active 任务单。

## 7. DoD

- [x] 所有 A-E 项完成或明确记录未完成原因。
- [x] 代表点输出可重复且包含完整输入、边界、状态和 flavor profile 字段。
- [x] derivative/by-parts 的差异已分解为 bulk、boundary 和 residual，不能再把漏边界的旧差直接称为 kernel mismatch。
- [x] `μ_s=0.2/0.55` 的影响能拆成 upstream-state 变化与同 state 上 signed `μ_K` 后处理变化。
- [x] 未改变正式 BU 默认、solver、mixed-meson、transport 或 baseline 语义。
- [x] 聚焦测试和治理检查通过，任务单状态与证据一致。

## 8. 代表点证据与 verdict

### 8.1 finite-window identity

代表点使用 `T=208 MeV`、`qmax=4 fm^-1`、`q_nodes=4`、`omega=[0.05,3] fm^-1`，对同一上游状态比较 `omega_nodes=12/24`：

| meson | current closure rel (12→24) | GBU closure rel (12→24) |
|---|---:|---:|
| `pi` | `0.257% → 0.063%` | `39.8% → 16.6%` |
| `K` | `0.435% → 0.118%` | `6.44% → 3.02%` |

旧 bulk-only current relative difference 在该点约为 `pi 17%`、`K 26%`；加入 boundary 后 current closure 降至千分量级，证明旧差异的主要缺项确为有限窗口边界。GBU closure 随节点加密下降但 `pi` 仍偏大，因此不把它升级为正式 kernel pass，也不把它解释为 kernel bug。

### 8.2 `mu_s/mu_q=0.2/0.55` 同点消融

代表点使用 `T=170 MeV`、`mu_q=80 MeV`、`mu_pi=100 MeV`、`qmax=4 fm^-1`、`q_nodes=4`、`omega=[0.05,3] fm^-1`、`omega_nodes=6`。每个 flavor profile 只求一个 upstream state；`mu_s/mu_q=0.2/0.55` 分别对应 `(mu_u,mu_d,mu_s)=(80,80,16)/(80,80,44) MeV`。

| channel/rule | `mu_s/mu_q=0.2` | `mu_s/mu_q=0.55` |
|---|---:|---:|
| `K+/pi+`, `mu_K=0` | `0.1159` | `0.0579` |
| `K-/pi-`, `mu_K=0` | `0.1684` | `0.0877` |
| `K+/pi+`, signed `mu_K` | `0.3514` (`mu_K=64 MeV`) | `0.2188` (`mu_K=36 MeV`) |
| `K-/pi-`, signed `mu_K` | `0.1050` (`mu_K=-64 MeV`) | `0.0999` (`mu_K=-36 MeV`) |

该结果把影响拆成两部分：`mu_K=0` 行显示 upstream-state 变化，signed 与 zero 的同 profile 差显示后处理化学势变化。原请求窗口对 `mu_pi=100 MeV` 属于 `strict_unsafe_support`；表中有限数值来自显式 `density_policy=:x_min_cut`、`bose_x_min=0.05` 的低节点诊断延拓，因此只用于机制和数量级观察。

### 8.3 结论

总 verdict：`diagnostic-only`。

- current finite-window identity 已通过代表点门禁；
- GBU 尤其 `pi` 仍需更高节点/导数 backend/phase branch 诊断后才能定量认证；
- `0.2/0.55` 对正负通道的影响不对称且足以改变 `O(10^-1)` 内的 ratio，支持进入固定 `(T,mu_B)` 的 `mu_Q,mu_S` 守恒荷求解阶段；
- 本结论不授权修改正式 BU kernel 或把 `0.2/0.55` 当成守恒荷平衡解。

### 8.4 验证记录

- `UNIT_FILES='relaxtime/test_bu_kernel_gate_utils.jl'`：`20/20` pass；
- `INTEGRATION_FILES='relaxtime/test_bu_analysis_script_contract.jl'`：`34/34` pass；
- `UNIT_FILES='relaxtime/test_meson_density.jl;models/test_flavor_chemical_profiles.jl'`：`149/149` pass；
- `scripts/dev/check_script_entrypoints.jl`：检查 `179` 个 Julia 文件，pass；
- `scripts/dev/check_relaxtime_script_governance.jl`：pass；
- `scripts/dev/check_docs_consistency.jl`、`git diff --check`：pass。

## 9. 风险与回退方案

1. **相移 unwrap 在端点与内点之间换支**
   - 缓解：端点和内点一次性排序 unwrap，并输出端点 phase 与跨度；若仍不稳定，verdict 降级为 `blocked`。
2. **folded/no-anomalous 相移不可直接做平滑 AD identity**
   - 缓解：本轮数学门禁只认证 unwrapped smooth branch；folded/no-anomalous 作为单独非光滑诊断，不据此修改正式 helper。
3. **全 48x48 charged phase-shift 诊断成本过高**
   - 缓解：先用低节点 smoke 验合同，再在一个代表点做收敛补档；不提前扩大到 T-ρ 全网格。
4. **边界项加入后 residual 仍大**
   - 缓解：依次检查 derivative backend、节点收敛、phase branch 和有限窗口，而不是直接调整 kernel normalization。
