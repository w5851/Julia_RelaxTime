---
title: Validation Transport / Phase Follow-ups
archived: true
original: docs/dev/active/2026-03-09_validation_transport_phase_followups.md
archived_date: 2026-03-09
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Validation Transport / Phase Follow-ups

## 0. 任务定位

本任务单承接 2026-03-08 两份 validation family active 文档归档后仍未完成、但不应阻塞归档的后续扩展项。

## 1. 当前挂账项

- `lorentz_legacy` legacy branch-divergence 记录：此前尝试升格的高 `mu_B` C++ fixed points 已确认不再适合作为 acceptance target；后续仅保留 provenance / evidence，不再作为默认验证项。
- crossover Fortran / C++ 双源交叉验证：旧 C++ `T-mu` crossover 导出已重生成并接入 dual-source guardrail；Fortran 侧 machine-readable crossover 导出也已补齐并正式接线到同一 dual-source guardrail。
- advanced `rates / decompositions` 准入规则：尚未定义哪些量可进入 acceptance、哪些只能作为 evidence 或 analysis。

## 当前进展

- 已确认旧 Fortran/C++ 在低 `mu_B` 仍能写出有限 `Lorentz` 的直接原因：旧工程按裸公式计算 `lambda` 与 `Lorentz`，不会像 Julia 当前 helper 那样在近零 `n_B` 处主动返回 `NaN`。
- 已确认 `lorentz_legacy` 的低 `mu_B=5 MeV` legacy anchor 不适合作为 acceptance 升级点；对应问题不是旧导出文件缺值，而是 Julia 当前 `lambda_from_kappa_BB` 在这些点触发了近零 `n_B` 保护。
- 已确认 `lorentz_legacy` 的 `muB=100 MeV, T=210/250 MeV` 候选点只保留为 provenance / evidence，不再继续作为 acceptance target 维护。
- 已从 Julia validation acceptance 中删除这两个 `lorentz_legacy` legacy 点，对应 guardrail 不再默认运行。
- 已新增 [tests/validation/data/provenance/relaxtime/evidence/relaxtime_lorentz_legacy_candidate_points_v1.csv](tests/validation/data/provenance/relaxtime/evidence/relaxtime_lorentz_legacy_candidate_points_v1.csv)，记录当前候选与拒绝状态点。
- 已完成 `lorentz_legacy` 高 `mu_B` 审计：Julia、旧 Fortran、旧 C++ 的符号公式链一致，剩余幅值差主要来自高 `mu_B` 比较点的平衡分支不一致。旧 C++ `muB=100 MeV, T=210/250 MeV` 落在低 `Phi`、高 constituent mass 的 hadron-like 分支；Julia 当前 helper 落在高 `Phi`、低质量的 crossover/quark-like 分支，因此 `kappa_BB/T^2` 与 `((e+p)/(n_B T))^2` 背景因子都会系统性偏离。
- 已形成策略判断：基于当前 `Ω` 比较与高温低密分支审计，不再把“让旧 Fortran/C++ 在 `muB=100 MeV, T=210/250 MeV` 上贴合 Julia”作为后续目标。更合理的解释是，这两个 legacy 点落在与当前 Julia 主线不同、且在 `T > 200 MeV` 低密区缺少热力学稳定性支撑的 branch 上，因此它们不再适合作为现代 acceptance target，只保留为 legacy branch-divergence evidence。
- 已重生成旧 C++ crossover `T-mu` 导出：`D:/Desktop/Cpp/2026_03_08/20250413备份/results/crossover_tmu_reference_cpp.dat`；该文件的 `muB_MeV` 列需先转换为 `mu_q = mu_B / 3`，再与 Julia `detect_crossover` 的输入口径对齐。
- 已完成 crossover 旧工程双源首轮接线：仓库内现成 reference 来源仍是 `data/reference/pnjl/crossover.csv`，同时已新增 legacy C++ dual-source guardrail fixed points。
- 已新增 [tests/validation/data/provenance/pnjl/evidence/pnjl_crossover_legacy_source_audit_v1.csv](tests/validation/data/provenance/pnjl/evidence/pnjl_crossover_legacy_source_audit_v1.csv)，记录双源盘点结论。
- 已新增旧 Fortran crossover machine-readable 导出程序 `D:/Desktop/Relaxtime_fortran/codes/main/AA_CROSSOVER_TMU_EXPORT.f90`，并成功生成 `D:/Desktop/Relaxtime_fortran/results/crossover_tmu_reference_fortran.dat`。当前 4 个 chiral fixed points 为 `muB=2/105/209/291 MeV` 对应 `T=200/199/196/193 MeV`，与已接线的旧 C++ fixed points 基本一致。
- 已把旧 Fortran crossover fixed points 正式并入 [tests/validation/data/targets/pnjl/reference/pnjl_crossover_legacy_dual_source_targets_v1.csv](tests/validation/data/targets/pnjl/reference/pnjl_crossover_legacy_dual_source_targets_v1.csv)，当前 legacy dual-source guardrail 已同时覆盖 C++ 与 Fortran 两套固定点来源。
- 已完成旧 Fortran / 旧 C++ crossover fixed-point 一致性核对：共享 `muB=2/105/209/291 MeV` 四点的 `T_c` 偏差分别为 `0.1/0.1/0.4/0.1 MeV`，现已新增独立 source-consistency target 并接入 Julia validation。
- 已对高温低密平衡解做直接 `Ω` 对比：在 `T=210 MeV, muB=0/30 MeV`，强子种子、夸克种子、弱禁闭种子和高温种子都会收敛到同一组解，`Ω` 一致到数值舍入误差内，`Phi ≈ 0.652`；在 `T=250 MeV, muB=0 MeV`，低温/强子型种子已不再收敛，只剩高温种子稳定落在 `Phi ≈ 0.762`、`phi_s ≈ -1.103` 的单一解上。这支持“`T > 200 MeV` 的低密区已不存在 competing hadron branch，高温区稳定相只能解释为 quark-like / deconfined branch”的判断。
- 已新增 [tests/validation/data/provenance/relaxtime/evidence/relaxtime_rates_decomposition_source_audit_v1.csv](tests/validation/data/provenance/relaxtime/evidence/relaxtime_rates_decomposition_source_audit_v1.csv)，记录 `rates / qk_sigmat / qk_etat3` 的当前源状态。
- 已确认 `rates / decompositions` 当前 checked-in 旧工程导出均为 0 字节：`sca_rate.txt`、`qk_sigmat.txt`、`qk_etat3.txt` 目前只能支持 schema audit，不提供新的 acceptance 信息增量。

## 2. 建议顺序

1. 先补 crossover 旧工程双源的一致性检查与 source-consistency guardrail。
2. 再把 legacy/reference 双源 family 继续收紧到 source-consistent 固定点。
3. 最后只在旧工程出现非空导出时再回看 `rates / decompositions`，当前不优先推进。

## 2a. 已确认结论

- `lorentz_legacy`：
	- 拒绝点：`muB=0 MeV, T=210/250 MeV`，因为 `lambda` 为 `NaN`。
	- 低 `mu_B=5 MeV` 旧锚点同样不适合作为 acceptance compare，因为 Julia workflow 给出的 `n_B` 只有 `1e-12 ~ 1e-11` 量级，触发了 `lambda_from_kappa_BB` 的近零 `n_B` 保护。
	- 已撤回 acceptance 升格：`muB=100 MeV, T=210/250 MeV` 这两个 legacy 点已从默认 validation target 中删除。
	- 当前策略更新：不再尝试调旧 Fortran/C++ 使其在 `muB=100 MeV, T=210/250 MeV` 上贴合 Julia，也不把这两个 legacy 点继续视作高价值 acceptance 基准。
	- 已完成首轮审计：Julia、旧 Fortran、旧 C++ 的 `lorentz_legacy = lambda / (sigma / T)` 公式链一致；剩余差异不是单一最终归一化因子，而是 `muB=100 MeV, T=210/250 MeV` 上的平衡分支不一致。旧 C++ 对应的是低 `Phi`、高 `m_u/m_s` 的 hadron-like 分支，而 Julia helper 当前落在高 `Phi`、低 `m_u/m_s` 分支。结合高温低密区 `Ω` 比较结果，这更适合解释为 legacy branch mismatch，而不是要求主线反向贴合旧分支。
	- 因而这两个点的推荐去向是：从“应继续收紧的 acceptance guardrail”转为“保留现有记录、但作为 evidence 的 legacy branch-divergence 注记”；后续验证资源优先投入 source-consistent 的 family。
- crossover 双源盘点：
	- Julia 主线已有 `data/reference/pnjl/crossover.csv`。
	- 旧 C++ `crossover_mu/crossover_T` 数组已通过专用导出器重生成 machine-readable 文件 `crossover_tmu_reference_cpp.dat`。
	- 该 legacy C++ 导出的横轴是 `mu_B`，而 Julia `detect_crossover` 当前使用的是 quark chemical potential，因此 dual-source guardrail 必须先做 `mu_q = mu_B / 3` 转换。
	- Fortran 侧现已新增 `AA_CROSSOVER_TMU_EXPORT.f90` 并产出 `crossover_tmu_reference_fortran.dat`，输出列为 `point_index, muB_MeV, T_MeV, phi_u, phi_d, phi_s, Phi, Phibar, m_u_MeV, m_d_MeV, m_s_MeV, delta_phi, check`。
	- 当前 Fortran 导出的 4 个 chiral fixed points 为 `200/199/196/193 MeV`，与旧 C++ dual-source guardrail 的 `200.1/199.1/196.4/192.9 MeV` 基本对齐，现已正式并入同一 dual-source guardrail。
	- 已新增 source-consistency validation：共享 `muB=2/105/209/291 MeV` 四点目前满足 `|ΔT| <= 0.5 MeV`。
	- 高温低密稳定分支补充结论：`T=210 MeV, muB=0/30 MeV` 上不同种子收敛到同一 `Ω` 极小点；`T=250 MeV, muB=0 MeV` 上强子型种子已无法收敛，只剩高温 branch。对当前 PNJL 参数而言，这说明 `T > 200 MeV` 的低密区没有可竞争的 hadron-like 平衡解，稳定相应视为 quark-like / deconfined branch。

## 3. 完成标准

- 每个新增 family 或派生量都要先回答：来源文件、字段语义、acceptance/evidence 边界、默认容差依据。
- 不允许把整条 raw scan 直接当默认测试输入。

## 3b. 当前固定约束

结合本轮已落地的 dual-source crossover 与 `lorentz_legacy` 降级处理，当前 validation 约束固定为：

- legacy 数据先进入 `provenance/.../evidence/`，不是先进入 `targets/`
- shared legacy fixed points 只有在 source-consistency 已验证后，才可升级为 acceptance guardrail
- 任何已知 branch-divergence comparison point 都不得再作为默认 acceptance target 继续维护
- raw scan、空文件导出、decomposition 大表只保留 analysis / evidence 角色，不直接成为默认测试输入
- `rates / decompositions` 当前维持 evidence-only，直到出现非空、可 machine-readable 读取且字段语义明确的 legacy 导出

上述约束已同步固化到 `tests/validation/README.md`，后续 family 扩展应直接遵守，而不是再次依赖活动文档口头判断。

## 3a. `rates / decompositions` 准入规则

### 可进入 acceptance 的条件

- 有稳定、非空、已检入的 legacy 或 reference 轻量导出。
- 字段语义能与 Julia 当前公共输出一一对应，而不是只在旧工程内部数组名层面存在。
- 可以抽成 fixed-point 或少量 representative points，不依赖整条大扫描文件直接入测。
- 容差可以基于旧工程对照、数值量级和低成本 validation 参数稳定性给出明确依据。

### 只进入 evidence 的条件

- 只确认了代码级输出 schema，但当前仓库中没有非空导出文件。
- 有导出文件但为空，或字段语义尚未完成 Julia 公共 API 对齐。
- 量本身更适合作为诊断输出，而不是稳定 acceptance 契约。

### 只进入 analysis 的条件

- 需要整条 raw scan、通道矩阵或大体积分解表才能解释。
- 当前 Julia 侧尚未把对应对象作为稳定公共输出公开。
- 缺少明确的 legacy/reference 锚点，强行比较会退化成实现自比较。

### 当前结论

- `rates`：目前先留在 evidence。Julia 公共 API 已有，但旧 C++ `sca_rate.txt` 在 checked-in snapshot 中为空。
- `qk_sigmat`、`qk_etat3`：目前先留在 evidence / analysis。旧 C++ 文件为空，且当前 Julia 主线也未把这些 flavor-resolved decomposition 作为正式 validation target 暴露。
- 结合当前 checked-in 文件均为 0 字节、且这些量的物理使用频率与公共 API 稳定性都弱于现有主 family，现阶段不建议优先补充这组验证。
- 若未来补到非空旧工程导出，优先先做 fixed-point target，而不是直接整表 compare。

## 4. 下一步建议

1. `lorentz_legacy` 方向已完成降级：后续不再围绕这两个 legacy 点追加对齐工作。
2. crossover 现已进入真正的 Fortran/C++ 双源 guardrail，并已补上 source-consistency 检查；下一步若继续收紧，应优先扩大 source-consistent fixed points，而不是再引入 branch 歧义点。
3. `rates / decompositions` 当前不建议优先补充验证，因为 checked-in 旧工程导出为空，且即便补齐，新增信息量也大概率低于现有主 family；除非后续出现非空且稳定的旧工程导出。
4. 若后续还要新增 transport guardrail，优先选择 source-consistent、branch-unique、不会落在高温 metastable branch 歧义区的固定点，而不是继续扩大 legacy branch mismatch family。

## 5. 收尾判断

基于当前 active / archived 文档状态，validation 这条工作流已经具备进入收尾准备的条件：

- 已完成本轮真正高价值的 family 落地：reference crossover、legacy dual-source crossover、legacy source-consistency gate
- 已把低价值且 branch-ambiguous 的 `lorentz_legacy` 点从 acceptance 中删除
- 已把 `rates / decompositions` 明确降到 evidence-only，不再作为阻塞项
- validation taxonomy、README 语义、family 目录结构和 follow-up 分工都已稳定

因此，当前更合理的动作不是继续扩大 backlog，而是准备把本 follow-up 文档转为“收尾摘要 + 残余 backlog 说明”。

当前建议的收尾口径：

- crossover mini family 视为已 formalized，可按版本化 target 维护
- `lorentz_legacy` 只保留 provenance / evidence 注记
- `rates / decompositions` 只有在未来出现非空旧工程导出时才重新激活

若本轮相关测试继续通过，则本活动文档已接近 archive-ready；后续只需在归档前补一版简短收尾摘要，不再需要新的大规模实现工作。