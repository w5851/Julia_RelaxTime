---
title: PNJL Crossover Validation Family 推进
archived: true
original: docs/dev/active/2026-03-08_pnjl_crossover_validation_family推进.md
archived_date: 2026-03-09
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# PNJL Crossover Validation Family 推进

## 当前进展

- 已确定第一批 crossover family 采用 reference fixed-point 口径，而不是整条扫描曲线直接入 acceptance。
- 已新增 `tests/validation/data/targets/pnjl/reference/pnjl_crossover_reference_targets_v1.csv`。
- 已新增 `tests/validation/pnjl/reference/test_crossover_reference_targets.jl`，覆盖 `chiral/deconf` 两条 crossover 线在代表性 `mu_B` 点的验证。
- 已新增 `tests/validation/data/provenance/pnjl/evidence/pnjl_crossover_reference_evidence_index_v1.csv`，记录来源文件与字段口径。
- 已完成 crossover 第一批 target 的定向验证与 validation 总入口验证。
- 已在 `tests/validation/README.md` 中补充 `pnjl/reference` 的语义与字段规范示例。
- 已抽出 PNJL phase/crossover 共用 helper，消除了 literature phase 与 crossover reference 测试中的重复 CSV 读取逻辑。
- 已将后续超出本轮 DoD 的 phase-family 深化项拆分到 follow-up 文档继续跟踪，不再阻塞本任务单归档。

## 0. 任务定位

本任务单承接 validation taxonomy 稳定化后的下一批 family 扩张，目标是把 PNJL crossover 从“脚本输出/参考文件”提升为正式 acceptance validation family。

## 1. 本轮目标

- 为 PNJL crossover family 建立第一批可维护的 acceptance target。
- 采用 fixed-point target，而不是把整条 `crossover.csv` 作为默认测试输入。
- 明确 `targets` 与 reference/provenance 的边界，避免把项目参考输出直接当作测试主输入。

## 2. 本轮范围

### 2.1 范围内

- `chiral crossover` fixed-point target。
- `deconfinement crossover` fixed-point target。
- 以 `mu_MeV` 代表点方式做 acceptance，保持测试成本可控。

### 2.2 暂不纳入

- 整条 crossover 扫描线逐点 compare。
- `rho_chiral / rho_deconf` 作为第一批 acceptance target。
- 与 CEP / first-order boundary 的统一 phase family helper 重构。

## 3. 第一批 target 设计

- 来源主表：`data/reference/pnjl/crossover.csv`
- acceptance 主表：`tests/validation/data/targets/pnjl/reference/pnjl_crossover_reference_targets_v1.csv`
- evidence 索引：`tests/validation/data/provenance/pnjl/evidence/pnjl_crossover_reference_evidence_index_v1.csv`

首批只保留四个代表性化学势点：

- `mu_B = 0`
- `mu_B ≈ 104 MeV`
- `mu_B ≈ 208 MeV`
- `mu_B ≈ 291 MeV`

每个点分别保留：

- `chiral_crossover`
- `deconf_crossover`

## 4. 验证策略

- 调用 `Models.detect_crossover(...)` 逐点复算，不直接读取整条 scan 输出。
- 当前第一批统一采用：
  - `method = inflection`
  - `solver_backend = :legacy`
  - 低成本 validation 参数：`p_num=8`, `t_num=4`, `n_scan=20`, `max_iter=20`
- 使用温度窗口而不是点对点完全相等，以吸收扫描离散度与求解器微小漂移。

## 5. 当前 checklist

- [x] 确定 crossover family 第一批采用 fixed-point acceptance。
- [x] 新增 crossover reference target CSV。
- [x] 新增 crossover validation test。
- [x] 新增 evidence 索引记录来源口径。
- [x] 完成定向验证与 validation 总入口验证。
- [x] 评估是否需要抽出 PNJL phase/crossover 共用 helper。
- [x] 在 `tests/validation/README.md` 中补充 `pnjl/reference/` 语义示例。

## 6. DoD

- [x] crossover family 已在 validation taxonomy 下拥有独立 target/test 落点。
- [x] 第一批 crossover fixed-point target 在本机验证通过。
- [x] 新贡献者可通过目录结构区分 acceptance target 与 reference/provenance 来源。

## 7. 归档准备摘要

- acceptance target 位于 `tests/validation/data/targets/pnjl/reference/`，采用从 `data/reference/pnjl/crossover.csv` 提取出的轻量 fixed-point 目标，而不是整条扫描线直接入测。
- provenance / evidence 位于 `tests/validation/data/provenance/pnjl/evidence/`，记录源文件与列语义。
- 当前第一批采用 `mu_B = 0, 104, 208, 291 MeV` 四个代表点，并分别验证 `chiral_crossover` 与 `deconf_crossover`。
- 统一验证口径：`Models.detect_crossover(...; method=:inflection, solver_backend=:legacy, p_num=8, t_num=4, n_scan=20, max_iter=20)`，使用温度窗口而非点对点相等。