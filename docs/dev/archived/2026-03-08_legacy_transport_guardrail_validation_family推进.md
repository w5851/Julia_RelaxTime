---
title: Legacy Transport Guardrail Validation Family 推进
archived: true
original: docs/dev/active/2026-03-08_legacy_transport_guardrail_validation_family推进.md
archived_date: 2026-03-09
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Legacy Transport Guardrail Validation Family 推进

## 当前进展

- 已新增第一批 legacy transport target：`tau`、`eta_over_s`、`sigma_t`（对应 legacy `sigma / T` guardrail 口径）。
- 已新增 `tests/validation/relaxtime/legacy/transport/` 下的三组测试文件。
- 已把 tau 复算主逻辑拆到共享 helper，旧 literature tau 测试已改为复用该 helper。
- 已新增第一批 provenance evidence CSV，记录本轮 fixed-point target 的来源口径。
- 已完成定向验证与 validation 总入口验证，当前 `tests/validation/runtests.jl` 可递归发现新增子目录测试并通过。
- 已纳入第二批 legacy transport family：`kappa_BB / kappa_QQ / kappa_SS / lambda`、`R_etasig`。
- 已将原有 `kappa/lambda` fixed-point 测试并入 `tests/validation/relaxtime/legacy/transport/` 目录语义下。
- 已确认 `zeta`、`lorentz_legacy`、`zeta_eta` 在当前 shared helper / workflow 口径下仍属于 secondary candidate：现有 acceptance helper 返回 `NaN`，本轮先保留 provenance 与 backlog，不强行纳入主 target。
- 已完成第二批 `kappa / lambda / R_etasig` 的定向验证与 validation 总入口验证。
- 已在 `tests/validation/README.md` 中补充 `relaxtime/legacy/transport` 与 `targets/provenance` 的字段语义说明。
- 已补上 bulk-enabled helper，并将 `zeta_over_s`、`zeta_eta` 接入 acceptance；当前 validation 总入口已提升至 `302 / 302`。
- 已将不再属于本任务单 DoD 的后续扩展项拆分到新的 follow-up 文档，不再阻塞本任务单归档。

## 0. 任务定位

本任务单作为 2026-03-08 validation taxonomy 重构工作的后续切片，目标不是继续做目录治理本身，而是在新 taxonomy 已稳定落地的前提下，推进下一批真正新增的 validation family。

本轮明确选择“第一批 legacy transport guardrail”作为切入点，而不是优先做 crossover。原因如下：

- 当前仓库刚完成 `kappa / lambda / ratio diagnostics / thermo background` 的实现与公式澄清，transport 语义边界最清楚。
- legacy Fortran / C++ 在 transport 侧已有更直接的可对照产物，生成 guardrail 的路径更短。
- 与现有 `eta_over_s` literature、`tau` literature、`kappa/lambda` fixed-point MVP 能自然衔接，便于形成第一批完整的 transport validation family。

## 1. 本轮目标

- 在 `tests/validation/` 新 taxonomy 下，新开第一批 legacy transport guardrail family。
- 先完成 transport MVP，而不是一次性纳入全部 legacy 高级数据。
- 将归档文档中仍未完成、且对当前 validation 扩张仍有价值的待办项，收编进本任务单统一推进。

## 2. 本轮范围

### 2.1 范围内

- 设计并落地第一批 legacy transport target 数据文件。
- 为第一批 transport guardrail 新增对应 validation tests。
- 补齐 transport target 的来源、容差、acceptance 与 provenance 边界。
- 继续整理 validation helper，使新增 family 不继续把逻辑堆回单个共享大文件。
- 明确 raw long 数据、targets、provenance 在 transport family 下的使用边界。

### 2.2 暂不纳入

- crossover family 的实际落地。
- `rates / decompositions` 作为正式 acceptance target 的接入。
- 一次性扩到全部 `zeta / kappa / lambda / ratios` legacy 曲线 target。

## 3. 推荐切片

本轮先做一个最小但有代表性的 legacy transport family MVP：

- `tau` legacy guardrail
- `eta_over_s` legacy guardrail
- `sigma_t` 或 legacy `sigma / T` 等价 guardrail

说明：

- `tau` 已有 literature validation，补 legacy guardrail 后可以形成“文献 + 旧实现”双来源稳定化闭环。
- `eta_over_s` 已有 literature 入口，补 legacy guardrail 的风险较低。
- `sigma` 在旧工程中有更直接的 guardrail 价值，且与已完成的 ratio / lambda 语义整理直接相关。

本轮不把 `Lorentz / R_etasig / Pr / zeta_eta / Ri` 直接作为第一批 acceptance 主 target，而是先作为第二层候选：

- 若其数值稳定且来源口径单一，可在同 family 下补充为 secondary guardrail。
- 若其数值受背景量或近零点行为影响较大，则先放 provenance 或 analysis，不强行纳入第一轮 acceptance。

## 4. 数据与测试的目标落点

### 4.1 目标数据位置

建议新增或扩展以下路径：

- `tests/validation/data/targets/relaxtime/legacy/transport/`
- `tests/validation/data/provenance/relaxtime/evidence/`
- `tests/validation/data/provenance/relaxtime/replacement/`

第一批 target 文件建议按“物理量 + 来源口径”拆分，而不是继续平铺到一个综合表：

- `relaxtime_tau_legacy_targets_v1.csv`
- `relaxtime_eta_over_s_legacy_targets_v1.csv`
- `relaxtime_sigma_t_legacy_targets_v1.csv`

若某一物理量需要保留 Fortran/C++ 双源对照证据，则额外保留：

- `relaxtime_transport_legacy_evidence_index_v1.csv`
- `relaxtime_<quantity>_legacy_consensus_points_v1.csv`

### 4.2 测试代码位置

建议在以下目录下新增第一批 test：

- `tests/validation/relaxtime/legacy/transport/`

测试文件粒度建议一物理量一文件，避免 transport family 再次堆成一个超大 test 文件：

- `test_tau_legacy_guardrail.jl`
- `test_eta_over_s_legacy_guardrail.jl`
- `test_sigma_t_legacy_guardrail.jl`

## 5. 承接自归档文档的高价值未完成项

以下待办来自已归档的 validation 重构文档，但其价值仍然直接服务于下一批 family 落地，因此保留并转入本任务单。

### 5.1 直接纳入本轮的待办

- [x] 设计第一批 `tau / eta / sigma` legacy target 数据文件。
- [x] 为新增 target 明确容差依据，避免拍脑袋设阈值。
- [x] 为新增 target 补齐 provenance、source、生成口径说明。
- [x] 明确哪些 transport 文件属于 acceptance target，哪些只属于 provenance only。
- [x] 明确 raw long 数据是否允许被 transport validation 测试直接依赖。
- [x] 记录 legacy 数据在 Fortran / C++ 间的一致性证据，并与 target 主表分离保存。
- [x] 起草 validation 测试目录迁移表，至少覆盖 relaxtime legacy transport 这一层新增路径。
- [x] 继续拆分现有共享 helper，避免把 transport family 的缓存、读取、容差逻辑继续堆入单个大文件。

### 5.2 作为本轮顺手补齐的治理项

- [x] 为 data taxonomy 补一份命名与字段约定文档，至少覆盖 transport legacy target 的字段规范。
- [x] 让目录结构本身能表达 transport 数据语义，而不是依赖文件名猜测。
- [x] 让新贡献者能直接区分 raw long、targets、provenance 三类 transport 数据。

### 5.3 已拆分到 follow-up 的后续项

- [x] 将 crossover 数据的 validation 口径与接受准则转入独立文档，并在本轮完成第一批 crossover family 落地。
- [x] 将 Fortran / C++ crossover 同点交叉验证盘点拆至 follow-up 文档继续跟踪。
- [x] 第二批 `zeta / kappa / lambda / ratios` 已完成当前稳定 acceptance 切片：`zeta_over_s`、`zeta_eta`、`kappa_BB / kappa_QQ / kappa_SS / lambda`、`R_etasig`。
- [x] 将 advanced `rates / decompositions` 准入规则拆至 follow-up 文档继续跟踪。

## 6. 实施步骤

### 阶段 A：transport family 口径冻结

- [x] 明确第一批 family 的状态点坐标体系。
  - 第一批 `tau / eta_over_s / sigma_t / zeta_over_s / zeta_eta / R_etasig` 统一采用 `T_MeV, muB_MeV` 固定点口径，当前 acceptance 点为 `muB=5 MeV` 且 `T=50,55,60 MeV`。
  - 第二批 `kappa_BB / kappa_QQ / kappa_SS / lambda` 继续采用 `T_MeV, muB_MeV` fixed-point 口径，当前 acceptance 点为 `(T, muB) = (210,0), (250,0), (210,100), (250,100)` MeV。
- [x] 明确每个 quantity 的比较口径。
  - `tau`：是点值、扫描线，还是代表性固定点。
  - `eta_over_s`：直接比无量纲值，不再额外折算。
  - `sigma_t`：显式写明与 legacy `sigma / T` 或 `sigma/T` 的对应关系。
- [x] 明确 acceptance 与 evidence 的分界。
  - 主 target 只保留测试直接消费的轻量点表。
  - Fortran/C++ 原始导出差异记录放 provenance。

### 阶段 B：数据生成与归位

- [x] 从 legacy Fortran / C++ 生成或整理第一批 transport guardrail 点。
- [x] 把生成后的轻量点表放入 `targets/relaxtime/legacy/transport/`。
- [x] 把分歧点、替换依据、双源对照索引放入 `provenance/relaxtime/`。
- [x] 如需保留长表或扫描中间结果，放入 `raw_long/relaxtime/transport/`，不直接作为默认 test 输入。
- [x] 将第二批 `kappa/lambda` fixed-point acceptance test 迁入 `relaxtime/legacy/transport/` 子目录语义。
- [x] 将稳定的 ratio diagnostic `R_etasig` 补入第二批 acceptance target。
- [x] 为 `zeta`、`zeta_eta` 明确单独的 bulk-enabled helper 并升格为 acceptance；`lorentz_legacy` 因当前低 `mu_B` legacy anchor 对应的 `lambda` 仍为 `NaN`，已拆至 follow-up 文档单独跟踪。

### 阶段 C：validation test 接入

- [x] 新增 `tests/validation/relaxtime/legacy/transport/` 下的第一批测试文件。
- [x] 将共用 CSV 读取、字段校验、容差比较逻辑下沉到新的 transport validation helper。
- [x] 验证 `tests/validation/runtests.jl` 在新增子目录下无需额外改动即可递归发现测试。
- [x] 让第二批 `kappa/lambda` 与 `R_etasig` 也沿用同一子目录递归发现路径。

### 阶段 D：治理文档补齐

- [x] 补 transport legacy target 的命名与字段约定说明。
- [x] 在 `tests/validation/README.md` 或等价文档中补 transport family 的语义示例。
- [x] 在实现完成后，把数据来源、生成命令、容差依据写入归档记录。

## 7. 字段与容差建议

第一批 legacy transport target 最小字段建议：

- `record_id`
- `series`
- `quantity`
- `T_MeV`
- `muB_MeV`
- `target_value`
- `rtol`
- 必要时 `atol`
- `source`
- `source_kind`

对于 legacy 共识 target，建议额外保留：

- `fortran_value`
- `cpp_value`
- `consensus_value`
- `relative_delta`
- `evidence_ref`

容差原则：

- 优先根据 Fortran/C++ 互相偏差、以及 Julia 对 legacy 的现状误差带来定标。
- 若某 quantity 在近零或近临界点附近波动大，则不能简单复用同一全局 `rtol`。
- 若 target 值穿零或量级跨度大，应允许按点给出 `atol` 或分段 `rtol`。

## 8. 验收标准

### 8.1 数据验收

- [x] 第一批 `tau / eta_over_s / sigma_t` legacy target 已落地到新 taxonomy。
- [x] 每类 target 都有可追溯来源和容差依据。
- [x] target 与 provenance 已明确分离，不再混表。
- [x] 第二批 `zeta_over_s / zeta_eta / kappa / lambda / R_etasig` target 已补入同一 taxonomy 语义下。

### 8.2 测试验收

- [x] 新增 test 能被 validation 总入口自动发现。
- [x] 第一批 transport guardrail 在本机口径下可稳定通过。
- [x] transport family 的新增测试不依赖 raw long 大表作为默认输入。
- [x] 第二批 `zeta_over_s / zeta_eta / kappa / lambda / R_etasig` guardrail 在本机口径下完成定向验证与总入口验证。

### 8.3 维护性验收

- [x] transport 相关 validation helper 不再继续堆入历史单大文件。
- [x] 新贡献者仅查看目录结构即可理解 transport targets 与 provenance 的区别。
- [x] 后续扩展 `zeta / kappa / lambda / ratios` 时无需再改大面积历史路径。

## 9. 里程碑

- M1：冻结第一批 legacy transport family 的 quantity 与状态点口径。
- M2：完成第一批 target / provenance 数据文件落地。
- M3：完成第一批 validation tests 接入并跑通。
- M4：完成字段规范、README 补充与归档收口。

## 10. DoD

- [x] 第一批 legacy transport guardrail family 已在新 taxonomy 下正式落地。
- [x] 至少 `tau / eta_over_s / sigma_t / zeta_over_s / zeta_eta` 五类 quantity 已形成稳定 validation 闭环。
- [x] target、raw_long、provenance 三层边界在 transport family 下是清晰且可执行的。
- [x] helper、字段规范、测试路径不再回退到旧平铺模式。
- [x] 归档材料中能回答：数据从哪里来、为什么这样设阈值、哪些是 acceptance、哪些只是 evidence。
- [x] 第二批 legacy transport guardrail family 已完成本机验证并进入稳定闭环。

## 13. 归档准备摘要

### acceptance 与 evidence 边界

- acceptance target：`tests/validation/data/targets/relaxtime/legacy/transport/` 下的轻量 fixed-point CSV，仅供 `test_*.jl` 直接读取。
- provenance / evidence：`tests/validation/data/provenance/relaxtime/evidence/` 下的来源索引，用于记录 legacy 文件路径、列语义、target 提取口径。
- raw long：当前 transport family 默认不作为测试输入；若未来保留长扫描表，仅进入 `raw_long/relaxtime/transport/` 供分析工作流使用。

### 本轮数据来源

- `tau`：`Relaxtime_fortran/results/rex_time.dat`
- `eta_over_s`、`sigma_t`、`R_etasig`：`Relaxtime_fortran/results/eta and sigma.dat` 与 `related_coe.dat`
- `zeta_over_s`、`zeta_eta`：`Relaxtime_fortran/results/zeta.dat`
- `kappa_BB / kappa_QQ / kappa_SS / lambda`：`Relaxtime_fortran/results/coe_diffusion.dat` 与 `coe_lambda.dat`

### 容差依据

- 当前 legacy transport fixed-point guardrail 统一采用宽松 `rtol=0.5` 或既有 `kappa/lambda` `rtol=0.05`，用于吸收 legacy/Julia 求解器差异与低成本 validation 参数带来的漂移。
- 对返回 `NaN` 或近零背景高度敏感的量，不强行纳入 acceptance；本轮明确将 `lorentz_legacy` 留在 follow-up backlog。

## 11. 风险与策略

### 风险

- legacy 导出量在不同程序中单位链可能不完全一致。
- 比值型量可能受背景量近零行为影响，直接做 acceptance target 风险较高。
- transport family 若一次纳入过多 quantity，容易重复上轮 taxonomy 之前的“平铺堆积”问题。

### 策略

- 先做最小 family MVP，再向 `zeta / kappa / lambda / ratios` 扩张。
- 对容易受近零背景污染的量，先做 evidence 或 secondary guardrail，不急于进入主 acceptance。
- 把单位、来源、双源偏差写进 provenance，而不是隐含在 test 常量里。

## 12. 备注

这份任务单明确承接已归档的 validation 重构文档中“仍值得继续做”的部分，但不再把“目录重构本身”当作主任务。当前主任务已经切换为：在已完成 taxonomy 治理的前提下，用第一批 legacy transport guardrail 验证这套新结构确实能持续承载后续 family 扩张。