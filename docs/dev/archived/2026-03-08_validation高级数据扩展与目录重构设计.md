---
title: Validation 高级数据扩展与目录重构设计
archived: true
original: docs/dev/active/2026-03-08_validation高级数据扩展与目录重构设计.md
archived_date: 2026-03-08
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Validation 高级数据扩展与目录重构设计

## 0. 当前收口状态

本轮工作已经先完成了一批 validation MVP 落地与稳定化，因此这份文档当前的重点不再是“要不要做 validation 扩展”，而是“如何把已经新增的数据和测试体系化收编”。当前已闭环内容如下：

- [x] legacy `usbar, mu_B=0` state guardrail 数据已落地
- [x] legacy `usbar, mu_B=0` sigma consensus curve/point guardrail 数据已落地
- [x] `kappa_BB / kappa_QQ / kappa_SS / lambda` fixed-point target 数据已落地
- [x] 对应 validation 测试已接入，`kappa/lambda` fixed-point 为 `29 / 29`
- [x] tau literature validation 的 `InterruptException` 已修复，修复方式限定在 validation 侧：降低测试专用积分节点并注入低成本 `cs_caches`，未改动生产默认参数
- [x] 全量 validation 当前通过 `187 / 187`

当前尚未完成的不是数值闭环，而是结构治理。2026-03-08 本轮推进后，首批结构收编已经落地：现有 validation 数据已按 `raw_long / targets / provenance` 三层归位，`runtests.jl` 也已升级为递归发现 `test_*.jl`。当前剩余工作主要转为后续 family 扩张、helper 继续拆分与新增 legacy/high-level dataset 准入。

## 1. 背景与目标

- 当前 `tests/validation/` 已经承载两类数据来源：
  - literature digitized / textual targets
  - legacy Fortran / C++ 共识 replacement / guardrail 数据
- 近期对老版本 Fortran 与 C++ 工程的勘查表明，这两套旧实现除了当前已使用的 `sigma(sqrt(s))` 和少量 equilibrium/state 点外，还能稳定产出一批更高层的物理结果：
  - 相变/交叉线数据（如 `crossover.dat` 或 C++ `loop_cro_*` 路径）
  - 弛豫时间谱
  - `eta`, `zeta`, `sigma`
  - `kappa_BB`, `kappa_QQ`, `kappa_SS`, `lambda`
  - Lorentz / Prandtl / `zeta_eta` / `R_etasig` / `Ri` 等比值型诊断量
  - 通道级平均散射率
  - 部分 flavor-resolved 拆分表
- 随着 validation 数据面扩张，当前 `tests/validation/data/` 的平铺结构已经开始混合以下不同性质的数据：
  - 原始长表 raw digitized 数据
  - 脚本直接消费的轻量 target 数据
  - disputed / provenance / evidence 辅助记录
  - legacy consensus guardrail 数据
- 需要一份统一开发文档，把“高级数据纳入 validation”与“validation 目录重构”放到同一执行框架里，避免后续继续增量堆文件。

## 2. 本文档目标

- 设计一条分阶段实施路径，把老版本项目可稳定生成的“高级数据”逐步纳入 Julia 项目的 validation 体系。
- 明确 `tests/validation/` 的新分层方案，区分：
  - long 原始数据
  - 脚本可用轻量数据
  - provenance / disputed / legacy evidence 数据
  - test 代码的模块化组织
- 为后续实现提供任务拆解、优先级、里程碑与验收标准。

## 3. 范围与非范围

### 3.1 范围

- 设计 validation 数据目录重构方案。
- 设计 validation 测试目录重构方案。
- 设计“老版本高级数据”纳入 validation 的阶段化路线。
- 明确首批应优先纳入的高级数据类型。
- 将 `crossover.dat` / 交叉线数据纳入正式候选范围。

### 3.2 非范围

- 本文档不直接实施目录迁移或测试改写。
- 本文档不一次性承诺把所有老版本高级数据在一轮内全部落地。
- 本文档不把旧工程中的一次性 debug / test 文件自动视为 validation 候选。
- 本文档不替代具体 dataset 生成规范、容差标定和回归阈值实验记录；这些需在后续实现文档或归档中补齐。

## 4. 现状盘点

### 4.1 当前 validation 顶层结构

- `tests/validation/runtests.jl`
- `tests/validation/pnjl/`
- `tests/validation/relaxtime/`
- `tests/validation/data/`

当前特点：

- 测试代码已经按 `pnjl` / `relaxtime` 做一级分层。
- 数据目录尚未进一步分层，所有数据文件平铺在 `tests/validation/data/`。
- `runtests.jl` 目前只按一级目录 include，不支持更细的递归测试组织。

### 4.2 当前数据目录中的数据类型混杂

当前 `tests/validation/data/` 已混合以下类型：

- PNJL literature textual / phase targets
- RelaxTime sigma / tau / eta_over_s literature targets
- sigma long raw digitized data
- tau long raw digitized data
- disputed points 记录
- legacy consensus replacement point 记录
- legacy curve guardrail 数据
- legacy equilibrium/state guardrail 数据
- `kappa/lambda` fixed-point target 数据

### 4.3 老版本项目已确认可产出的高级数据面

已确认可作为 validation 候选的高层数据包括：

- 相变 / crossover 线数据
  - Fortran：`crossover.dat` 路径存在但默认主循环中当前注释，需要专项驱动或恢复产线
  - C++：`loop_cro_tmu()` / `loop_cro_tn()` 已存在，可作为交叉线候选入口
- 平衡态与相变量
  - `solution` / `A0 solution`
- 热力学状态量
  - `pse` / `A1 quantity`
  - `sound_speed`
- 弛豫时间
  - `rex_time`
- 输运系数
  - `eta`
  - `zeta`
  - `sigma`
  - `kappa_BB`, `kappa_QQ`, `kappa_SS`
  - `lambda`
- 无量纲比值与诊断量
  - Lorentz
  - `R_etasig`
  - Prandtl
  - `zeta_eta`
  - `Ri`
- 通道级平均散射率
  - Fortran：`w_ij*`
  - C++：`sca_rate.txt`
- 拆分型输运贡献表
  - C++：`qk_etat3.txt`, `qk_sigmat.txt`

## 5. 设计原则

- 单一真相源：每类 validation 数据在仓库内只保留一个权威位置。
- 数据语义优先于历史文件名：目录结构先反映“数据性质”，再反映“来源历史”。
- raw 与 target 分离：原始长表不可与脚本直接消费的小 target 表混放。
- provenance 显式化：disputed、replacement、legacy consensus、evidence 记录必须单独存放。
- 分阶段落地：先完成结构治理和高价值 MVP 数据，再扩展到全量高级数据。
- 轻量脚本优先：测试脚本默认只依赖轻量 target 数据，不直接依赖大长表，除非该测试明确属于 raw compare 类。

## 6. 目标目录结构设计

### 6.1 数据目录目标结构

建议把 `tests/validation/data/` 重构为：

```text
tests/validation/data/
  raw_long/
    pnjl/
    relaxtime/
      sigma/
      tau/
      transport/
      phase/
  targets/
    pnjl/
      literature/
      legacy/
    relaxtime/
      literature/
        sigma/
        tau/
        eta_over_s/
        phase/
        transport/
      legacy/
        sigma/
        state/
        phase/
        transport/
  provenance/
    pnjl/
    relaxtime/
      disputed/
      replacement/
      evidence/
```

解释：

- `raw_long/`
  - 存原始长表、digitized 全量表、扫描长表等。
  - 原则上不直接作为绝大多数 test 的主输入。
- `targets/`
  - 存脚本直接消费的轻量 acceptance / guardrail 数据。
  - 后续大多数 `test_*.jl` 仅从这里读取。
- `provenance/`
  - 存 disputed、replacement、legacy consensus evidence 等治理类文件。
  - 这些数据用于审计和 provenance，不直接混入 target 主路径。

### 6.2 测试代码目录目标结构

建议把 `tests/validation/pnjl/` 与 `tests/validation/relaxtime/` 继续细分：

```text
tests/validation/
  common/
    io_helpers.jl
    provenance_helpers.jl
  pnjl/
    literature/
    legacy/
    phase/
  relaxtime/
    common/
    literature/
      sigma/
      tau/
      eta_over_s/
      phase/
      transport/
    legacy/
      sigma/
      state/
      phase/
      transport/
    advanced/
      rates/
      decompositions/
```

解释：

- `common/` 放共享 helper，避免 `literature_validation_helpers.jl` 越长越难维护。
- `literature/` 与 `legacy/` 分开，避免 acceptance 基准语义混淆。
- `advanced/` 用于 flavor-resolved / rate-resolved / decomposition guardrail，和主 acceptance target 分层。

### 6.3 测试入口调整方向

- `tests/validation/runtests.jl` 需从“只 include 一级目录”升级到“递归 include 子目录中的 `test_*.jl`”。
- 递归 include 时应保持：
  - `PNJL Baselines`
  - `RelaxTime Baselines`
  两个顶层 testset 不变，避免 CI 输出面突变过大。

补充现状：当前 `tests/validation/runtests.jl` 仍停留在一级 include 模式，这也是下一轮结构治理最直接的落点之一。

## 7. 高级数据纳入 validation 的候选清单与优先级

### 7.1 P0：应最先纳入

- [ ] 相变 / crossover 线数据
  - 理由：用户已明确提出；且这是独立物理面，不与现有 sigma/tau/eta 数据重复。
  - 候选来源：Fortran `crossover.dat` 路径、C++ `loop_cro_tmu()` / `loop_cro_tn()`。
- [ ] 弛豫时间谱 guardrail
  - 理由：仓库已有 literature tau validation，补 legacy 高级数据最自然。
- [ ] `eta/s` 与 `sigma/T` 或 `sigma` guardrail
  - 理由：已有 `eta_over_s` literature validation，可自然扩到 legacy 共识 guardrail。

### 7.2 P1：第二阶段纳入

- [ ] `zeta/s` 或 `zeta/T^3`
- [x] `kappa_BB/T^2`, `kappa_QQ/T^2`, `kappa_SS/T^2` 的 fixed-point MVP 已有首版 target 与测试入口
- [x] `lambda/T^2` 的 fixed-point MVP 已有首版 target 与测试入口
- [ ] Lorentz / `R_etasig` / Pr / `zeta_eta` / `Ri`

### 7.3 P2：高级拆分与机制守卫

- [ ] 通道级平均散射率 `w_ij` / `sca_rate`
- [ ] flavor-resolved `qk_etat3`, `qk_sigmat`
- [ ] 其他明确可复现、且不属于一次性 debug 的 decomposition 数据

## 8. 推荐实施顺序

### 阶段 A：结构治理先行

- [x] 设计并落地新的 `tests/validation/data/` 目录树。
- [x] 先完成 `runtests.jl` 递归 include 改造，为后续测试目录细分清障。
- [x] 制定命名约定：raw / target / provenance / legacy / literature 的统一前缀与后缀规则，并写入 `tests/validation/README.md`。
- [x] 制定迁移兼容策略：本轮直接将 test 入口切到新路径，旧平铺路径不再作为主入口保留。

### 阶段 B：先落 crossover 与第一批 transport guardrail

- [ ] 选定 crossover validation 的最小物理口径：
  - `T(mu_B)` 线点
  - 或固定 `mu_B` / 固定 `T` 的 crossing targets
- [ ] 生成第一版 crossover legacy / literature target 数据。
- [ ] 设计 `eta/s`、`tau`、`sigma` 的 legacy target 轻量表。
- [ ] 新增对应 test 文件并完成容差标定。

补充现状：`tau` literature target 已有测试，但这一轮新增的是“稳定化闭环”而不是 taxonomy 迁移；下一步应优先补 legacy tau/transport guardrail 和目录归位，而不是继续在旧目录平铺追加文件。

### 阶段 C：扩展到 transport full family

- [ ] 纳入 `zeta`、`kappa`、`lambda`。
- [ ] 纳入 ratio diagnostics。
- [ ] 为每类数据建立 provenance / source / tolerance 说明。

### 阶段 D：扩展到 advanced guardrail

- [ ] 评估 `w_ij` / `sca_rate` 是否作为正式 guardrail 还是 analysis-only baseline。
- [ ] 评估 `qk_etat3` / `qk_sigmat` 是否作为 legacy advanced validation。
- [ ] 将不适合主 acceptance 的拆分型数据放入 `advanced/` 层，不与主 target 混合。

## 9. 文件命名与数据规范建议

### 9.1 raw long 数据

- 命名建议：
  - `*_digitized_longtable_v1.csv`
  - `*_scan_longtable_v1.csv`
  - `*_legacy_longtable_v1.csv`

### 9.2 target 轻量数据

- 命名建议：
  - `*_literature_targets_v1.csv`
  - `*_legacy_consensus_targets_v1.csv`
  - `*_legacy_state_targets_v1.csv`
  - `*_crossover_targets_v1.csv`

### 9.3 provenance 数据

- 命名建议：
  - `*_disputed_points_v1.csv`
  - `*_replacement_points_v1.csv`
  - `*_evidence_index_v1.csv`

### 9.4 每类 target 的最小字段集合

- `record_id`
- `series`
- 核心状态点坐标，如 `T_MeV`, `muB_MeV`, `sqrt_s_GeV`
- `target_value` 或结构化的字段值
- `rtol` / 必要时 `atol`
- `source`
- 对 legacy 数据，额外建议保留：
  - `fortran_value`
  - `cpp_value`
  - `consensus_value`
  - `relative_delta`

## 10. 任务分解

### 10.1 结构设计任务

- [x] 完成 `tests/validation/data/` 首批重构，现有 literature / legacy / kappa-lambda CSV 已归位到新 taxonomy。
- [ ] 起草 `tests/validation/` 测试目录迁移表。
- [~] 设计新的 helper 布局，已新增 `tests/validation/common/data_paths.jl`，但 `literature_validation_helpers.jl` 仍待继续拆分。
- [x] 完成 `runtests.jl` 的递归 include 改造。

### 10.2 高级数据纳入任务

- [ ] 定义 crossover 数据的 validation 口径与接受准则。
- [ ] 盘点 Fortran / C++ 在同一物理点上的 crossover 可交叉验证程度。
- [ ] 设计第一批 `tau / eta / sigma` legacy target 数据文件。
- [ ] 设计第二批 `zeta / kappa / lambda / ratios` target 数据文件。
- [ ] 设计 advanced `rates / decompositions` 是否纳入 validation 的准入规则。

### 10.3 文档与治理任务

- [x] 为新目录结构补一份 `tests/validation/README.md` 或等价说明文档。
- [ ] 为 data taxonomy 补一份命名与字段约定文档。
- [ ] 明确哪些文件属于 acceptance target、哪些属于 provenance only。
- [ ] 明确 raw long 数据是否允许被测试脚本直接依赖。

## 11. 测试与验收标准

### 11.1 结构验收

- [x] `tests/validation/data/` 中 raw / target / provenance 不再平铺混放。
- [ ] 新增数据类型后，目录即可直接表达其语义，不需要依赖文件名猜测。
- [x] `runtests.jl` 能稳定递归发现新的 `test_*.jl`。

### 11.2 数据验收

- [ ] 每类新增高级数据都有明确来源说明。
- [ ] 每类新增 target 都有容差依据，不以拍脑袋方式设阈值。
- [ ] legacy 数据在 Fortran / C++ 间的一致性证据被单独记录。
- [ ] disputed / replacement / evidence 数据不再与 target 主表混放。

### 11.3 维护性验收

- [ ] 新贡献者能仅通过目录结构区分 raw long 数据与脚本可用轻量数据。
- [ ] validation 代码的 helper 不再集中堆叠在单个大文件中。
- [ ] 新增一个 validation family 时，不需要再次改动大面积历史路径。

## 12. 里程碑

- M1：完成目录重构设计评审
- M2：完成 crossover + 第一批 transport legacy targets 的 MVP 落地
- M3：完成 transport full family 轻量 targets 接入
- M4：完成 advanced rates / decompositions 的准入与分层

## 13. DoD

- [x] 新的 validation 数据 taxonomy 已在仓库文档中固定。
- [ ] 至少一类高级数据（优先 crossover）已经按新结构落地。
- [x] `tests/validation/` 测试入口已支持更细目录分层。
- [x] 现有 literature / legacy 数据已按新语义路径归位。
- [ ] 新增高级数据的 provenance、容差依据、生成口径均可追溯。

补充收口判断：

- [x] 结构治理第一阶段已完成，当前未勾选项主要属于“下一批 validation family 准入”，而不是本轮目录重构残缺。

## 14. 风险与回退方案

### 风险

- 一次性迁移过多文件，导致路径回归面太大。
- 老版本高级数据的单位与口径未完全统一，直接接入会制造新的假失败。
- 通道率 / 分解类数据过于敏感，不适合直接作为 acceptance target。
- crossover 数据在 Fortran 与 C++ 侧的可复现路径可能不完全对齐。

### 回退方案

- 先只做目录设计与最小迁移，不一次性搬空全部历史数据。
- 每新增一种高级数据，先走单 family MVP，再复制模式。
- 对敏感数据先放到 `advanced/` guardrail，而不是主 acceptance target。
- 对尚未统一单位/口径的数据，先只建立 provenance / inventory，不立即写 test。

## 15. 当前推荐的第一实施切片

在当前工作树状态下，建议下一轮不要再从“新增更多目标点”开始，而是先做结构收编。更合适的第一切片是：

- [x] 当前已经存在的 state / sigma legacy / kappa-lambda / tau 文献数据已完成首轮新 taxonomy 归位。
- [x] `runtests.jl` 已完成递归 include 改造，使 `tests/validation/relaxtime/` 能继续向 `literature/legacy/advanced` 细分。
- [ ] 继续拆分现有共享 helper，把 tau 稳定化里引入的 cache/参数策略收敛到共享 validation helper，而不是继续留在单文件局部逻辑里。
- [x] 已补 `tests/validation/README.md`，明确 raw/targets/provenance 三层与文件命名规则。
- [ ] 在结构收编完成后，再以 `crossover` 或第一批 legacy transport guardrail 作为新增 family 进入点。

这样可以先把已经验证通过的成果纳入稳定目录结构，再继续扩张数据面，避免下一轮对话继续在旧路径和平铺数据上累积技术债。