# PR290：公式路线闭合 SOP 与 charged-RPA/BU 试点

更新日期：2026-08-30
当前状态：in progress（等待独立审阅；不授予数值 production 资格）
独立基线：`origin/main` @ `bc9b2990bcfe3b8c32d2ec0f00066b52b4cf800b`

## 1. 背景与目标

PR289 保持 open。本任务从合并后的 `origin/main` 独立开始，落实“每次向新领域
拓展都先闭合唯一公式路线”的 SOP 约束，并以 charged-RPA/BU 作为第一个试点，
为后续 5.6Sol 审核提供固定、可逐项质询的公式包。

目标不是在本 PR 直接完成严格 charged-RPA/BU 数值生产，而是把微观拉氏量、
平均场、`K_ab`、有序泡、RPA 分母、相移、BU 权重、输出字段和未决项全部接到
同一份可审阅规范上。

## 2. 范围与非目标

### 范围

- 新增跨领域公式路线闭合 SOP、registry 和可执行 checker；
- 编写 fixed-BQS quark-only charged `pi^\pm/K^\pm` 的 candidate 公式包；
- 记录外部来源、项目符号/单位转换、代码/测试映射和 production 边界；
- 增加纯结构/代数契约测试，不改变已有数值 baseline。

### 非目标

- 不修改 `PNJLCore`、`MesonDensity` 稳定默认语义或 transport；
- 不实现 strict-support/凝聚零模、复杂极点求解或完整 `Omega_M` 反馈；
- 不把 `K03/K38` 注入当前 charged 通道；
- 不提交任何既有或新生成的 diagnostic CSV，不更新 production baseline；
- 不在本 PR 创建独立 skill；待流程在三个以上新领域复用后再评估。

## 3. 当前能力与缺口

已有能力：

- `MesonInteractionKernel` 已提供 `K12/K45/K67` 和中性 `(0,3,8)` 纯代数后端；
- `MesonRPA` 已提供中性矩阵 `2K[I-2KPi]^{-1}` 的诊断代数；
- `MesonDensity` 已能在固定背景注入 charged coupling，沿用 `2K/(1-4KPi)`；
- 现有文档/测试已锁定 `K^\pm -> K45`、旧 `K4567 -> K67` 的基本映射。

仍需严格复核：charged `Pi_us/Pi_su` 的 retarded 延拓、scalar/matrix 归一化
转换、极点/相位边界、Levinson、Bose 支撑和节点收敛。因此试点 route 保持
`candidate`。

## 4. 可勾选任务分解

- [x] 新增 `formula_route_closure` SOP，并登记到 `docs_authority_map.toml`。
- [x] 在 literature-to-implementation protocol 中加入公式闭合门禁。
- [x] 新增 ADR-0007，记录治理决策和 charged-RPA/BU 试点边界。
- [x] 新增 `ChargedRPA_BU_ProductionRoute.md`，形成模型→输出的闭合公式链。
- [x] 新增 `formula_route_closure.toml` 与 `check_formula_route_closure.jl`。
- [x] 增加 route registry 纯契约测试，并接入 unit core profile。
- [x] 将本 PR 登记到 `config/governance/task_tracks.toml`，保留独立基线和审阅状态。
- [ ] 由独立审阅者逐项复核公式、来源方程号和单位/归一化转换。
- [ ] 在后续任务中实现并验证严格 charged retarded bubble、极点和 BU gate。

## 5. 测试与验收标准

本 PR 至少通过：

```powershell
julia --project=. scripts/dev/check_formula_route_closure.jl
julia --project=. scripts/dev/check_sop_governance.jl
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. -e 'ENV["UNIT_FILES"]="config/test_formula_route_closure.jl"; include("tests/unit/runtests.jl")'
```

纯契约测试应拒绝缺失来源、缺失章节、无效状态、路径越界和
`production_authorized`/状态矛盾。上述检查通过只表示文档/登记闭合，不表示
数值收敛或物理 production 正确。

## 6. 里程碑

1. **M0 基线冻结**：确认独立于 PR289 的 `origin/main` SHA 和干净工作树。
2. **M1 治理落地**：SOP、registry、checker、ADR 和模板约束可执行。
3. **M2 公式试点**：charged-RPA/BU candidate 公式包可由独立审阅者逐式复核。
4. **M3 PR290**：选择性提交并创建 draft PR；strict 数值实现另立后续任务。

## 7. DoD

- 公式路线有唯一 `route_id`、明确状态和冻结基线；
- 每条外部公式有 DOI/arXiv、方程号/页码和与项目规范的转换说明；
- 关键公式均有代码/测试映射，未决项和非目标没有被隐藏；
- checker、SOP governance、docs consistency 和纯契约测试通过；
- PR290 不依赖 PR289，也不修改已有 production/reference 数据。

## 8. 风险与回退方案

- 若独立审阅发现符号或归一化不闭合：保持 route `candidate`，修正文档/测试，
  不修改旧数值接口；
- 若外部来源版本或方程号无法确认：保留待核验项，不把来源写成已证实事实；
- 若后续严格实现需要改变平均场、混合基底或 `Omega_M`：创建新 route id 和新
  决策记录，不在本试点静默扩大范围；
- PR289 继续保持 open，直到本 PR 和后续审阅明确其改动是否仍有必要。
