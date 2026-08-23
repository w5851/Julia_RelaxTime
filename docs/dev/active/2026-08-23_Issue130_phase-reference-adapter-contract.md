# Issue #130：phase-reference adapter contract 与引用迁移

状态：active。PR253 的 solver-free compatibility audit 已完成并接受；本任务把
versioned candidate 的 schema 语义收束为显式 adapter 合同。默认 legacy reference
和所有旧文件保持不变，当前阶段不切换 runtime、不启动 RS transport、不运行 PNJL。

## 已知边界

- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`
- layer：`strict`、`derived`；`render` 仅可用于绘图，不可作为 runtime 输入。
- `mu_MeV`/`mu_CEP_proxy_MeV` 是 `mu_q`；adapter 同时派生 `muB=3*mu_q`，不得由调用方猜测。
- Maxwell 的 `grid_unresolved`、`geometry_converged`、`finite_and_converged` 和状态字符串
  共同决定 `certified`；unresolved 行必须保留为诊断，但默认不能进入 legacy runtime view。
- derived 的 `interpolated_noncertified` 只能显式允许用于派生图或诊断，不能伪装成 strict。
- CEP 同时保留 `[T_low_MeV,T_high_MeV]`、宽度和 `T_midpoint_MeV`；使用中点必须显式选择
  `cep_mode=estimated_midpoint`，内部审计可以选择 `cep_mode=bracket`。
- crossover 仅在 `physical_region=crossover_below_CEP` 时赋予 crossover 语义；Maxwell 区的
  导数峰不能被 adapter 自动暴露为 crossover。

## 任务项

- [x] 新增 solver-free `scripts/pnjl/phase_reference_adapter.py`，校验 import/layer manifest、
  CSV 字段、finite、重复键和 runtime_consumption gate。
- [x] 归一化 boundary/crossover/CEP/spinodal，明确 `mu_q`/`mu_B` 和 strict/derived 语义。
- [x] 提供只在内存中生成的 legacy view；默认拒绝 unresolved/non-certified 行，不写旧 reference。
- [x] 增加 fixture 单元测试，覆盖三态、CEP 双模式、单位换算、重复键、NaN/Inf 和 runtime gate。
- [x] 将 plotting consumer 改为显式 candidate root/layer 入口，默认保持 legacy；candidate
  入口只生成诊断图，不改变 runtime 语义。
- [x] 为 `run_gap_transport_scan.jl` / phase-guided plan 定义 Julia 侧等价 adapter，并完成
  legacy-vs-candidate parity fixtures；这一步完成前不改变 transport 默认路径。
- [x] 审核 Paper P1 的 tagged-file contract，并提供显式 candidate overlay 与 parity fixture；
  candidate 仍只用于显式诊断/图层输入，不改变 P1 默认 legacy 路径。
- [x] 完成 transport、phase-guided、paper 和 legacy plot 的显式引用入口迁移及 solver-free
  回归；默认 runtime 仍为 legacy，candidate 必须显式加载并通过 certified gate。
- [ ] 另立 versioned runtime-switch PR，默认使用 `strict` candidate，保留 legacy
  fallback/rollback；在作者审核和 switch 合并前不得删除旧 `boundary.csv`、`cep.csv`、
  `spinodals.csv`、`crossover_dense.csv` 或 `phase_reference_dense_manifest.json`。

## 验收与非目标

本阶段验收是 adapter contract 可读、可拒绝不安全输入、消费者 parity 可追溯且默认不改变
runtime；不是 phase-reference promotion，也不是 RS transport production。验证只运行 Julia/Python
fixture、candidate schema replay、CLI/dry-run、`git diff --check`、task-ledger/docs governance；
本机不调用新的 equilibrium solver 扫描。下一步必须先创建独立 runtime-switch PR，经过 CI 和
作者审核后才能标记 candidate 为默认 runtime reference；旧 reference retirement 仍需另立 PR。
