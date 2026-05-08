# PNJL dense reference CI 与 artifact 固化任务单

更新日期：2026-05-08

当前状态：设计中；用于承接 `crossover_dense.csv` / `crossover_dense.meta.json` 后续扩展与 GitHub Actions 化，目标是把“长时间运行的 dense reference 生成”从本机手工执行迁移到可手动触发、可下载 artifact、可审计 provenance 的 CI 工作流。

> 目的：把 dense phase reference 的生成、校验、artifact 留存与正式落库路径从“本机长跑脚本”收束成可复用的 CI 流程，而不是继续依赖个人机器长时间运行和手工拣选文件。

---

## 1. 背景与目标

当前仓库在 PNJL phase / dense reference 方向已经具备：

- [x] 手动触发的 PNJL phase diagram workflow：
  - `.github/workflows/pnjl-phase-diagram.yml`
- [x] dense reference 生成脚本：
  - `scripts/pnjl/build_dense_phase_reference.jl`
- [x] dense crossover reference 与 sidecar meta：
  - `data/reference/pnjl/crossover_dense.csv`
  - `data/reference/pnjl/crossover_dense.meta.json`
- [x] 相图邻域 transport 当前已能消费 dense reference：
  - `scripts/relaxtime/phase_guided_transport_scan_plan.jl`

但当前仍缺一条正式的 dense reference CI 主线：

- [ ] 本机运行 `build_dense_phase_reference.jl` 时间较长，不适合频繁手工重跑
- [ ] 当前没有专门面向 dense reference 的 workflow_dispatch 入口
- [ ] 当前没有围绕 dense reference 的 artifact 命名、manifest、下载与“是否晋升到 `data/reference/`”治理口径
- [ ] 当前没有把“生成 artifact”与“正式落库到仓库”区分成两阶段流程

本任务单的目标是：

1. 设计 dense reference 的 GitHub Actions 手动触发工作流；
2. 明确是否以及如何复用现有 `PNJL Phase Diagram` CI；
3. 固化 artifact 命名、保留期、manifest、校验与晋升策略；
4. 为后续更全 `xi` 网格 / 更高成本 reference 生成提供非本机长跑路径。

---

## 2. 现状盘点

### 2.1 现有 workflow 能力

- [x] `.github/workflows/pnjl-phase-diagram.yml` 已具备以下可复用骨架：
  - `workflow_dispatch`
  - Julia 1.12.5 + Python 3.11 环境搭建
  - `julia-actions/cache@v2`
  - artifact 上传
  - step summary
- [x] GitHub 上该 workflow 最近成功运行日期为：
  - `2026-01-04`
  - `2026-01-03`
- [x] 因此它不是废弃 workflow，而是“仍能运行，但目标仍停留在旧相图产物”的基础设施

### 2.2 现有 workflow 与 dense reference 的错位

- [x] 当前 workflow 的核心执行脚本是：
  - `scripts/pnjl/calculate_phase_structure.jl`
- [x] dense reference 需求的核心执行脚本是：
  - `scripts/pnjl/build_dense_phase_reference.jl`
- [x] 当前 workflow 以 `matrix.xi` 为主，按单个 `xi` 拆任务
- [x] dense reference 更自然的单位是“单次配置生成一套聚合 reference + meta + manifest”
- [x] 当前 workflow 会清理：
  - `data/outputs/figures/pnjl/phase_diagrams/*.png`
- [x] 当前 workflow 上传：
  - `data/reference/pnjl/*.csv`
  - `data/outputs/figures/pnjl/phase_diagrams/*.png`
- [x] 这些路径与 dense reference 正式治理并不完全匹配：
  - 缺少 `.meta.json` / manifest 的显式 artifact 契约
  - 上传范围过宽，不利于区分“本次生成了哪些 reference”
  - 没有“生成 artifact”与“晋升到仓库正式文件”的分层

### 2.3 当前脚本能力

- [x] `build_dense_phase_reference.jl` 已支持：
  - `--xi-list`
  - `--xi-min/--xi-max/--xi-step`
  - `--tag`
  - `--output-root`
  - `--reference-root`
  - `--overwrite`
  - `--crossover-only`
  - `--crossover-mu0-only`
  - `--crossover-n-mu`
  - `--crossover-mu-max`
- [x] 脚本已自动写出：
  - `crossover_<tag>.csv`
  - `crossover_<tag>.meta.json`
  - `phase_reference_<tag>_manifest.json`
- [x] 这意味着 dense reference 的 CI 化不需要先发明新输出协议，重点在 workflow 编排与 artifact 治理

---

## 3. 核心判断：能否基于 `PNJL Phase Diagram` CI 实现

### 3.1 结论

- [x] **可以基于它的骨架实现**
- [x] **但不建议直接在 `.github/workflows/pnjl-phase-diagram.yml` 上原地改造为 dense reference 正式流**
- [x] **更推荐新建专用 workflow**，例如：
  - `.github/workflows/pnjl-dense-reference.yml`

### 3.2 为什么不是直接复用原 workflow 文件

- [x] 原 workflow 的用户心智是“相图计算与出图”
- [x] dense reference 的用户心智是“生成 reference artifact，等待人工审阅后再决定是否正式落库”
- [x] 原 workflow 的 artifact 粒度偏粗，且上传路径过宽
- [x] 原 workflow 的输入围绕 `T/rho` 网格，不围绕 `tag / crossover-only / xi dense grid / reference-root`
- [x] 原 workflow 当前还带有 phase diagram figure 生成步骤，这对 dense reference 并非必要
- [x] 若在原文件上继续堆参数，会让 workflow 语义变得混杂，不利于后续维护

### 3.3 哪些部分应该直接借用

- [x] `workflow_dispatch` 触发方式
- [x] Julia / Python runtime setup
- [x] Julia cache 策略
- [x] artifact upload
- [x] summary 输出模式

### 3.4 哪些部分应该重写

- [ ] 输入参数设计
- [ ] 执行脚本
- [ ] artifact 路径白名单
- [ ] summary 字段
- [ ] 结果晋升策略

---

## 4. 建议范围与非范围

### 4.1 本期范围

- [ ] 设计专用 dense reference workflow_dispatch 入口
- [ ] 设计 artifact 命名与 retention 策略
- [ ] 设计基于 `tag` 的输出目录与 manifest 汇总
- [ ] 设计“生成 artifact”与“人工晋升到仓库 reference 文件”的两阶段流程
- [ ] 设计最小校验步骤（文件存在、行数、meta/manifest 一致性）
- [ ] 明确如何复用现有 `pnjl-phase-diagram.yml` 的公共骨架

### 4.2 本期不覆盖

- [ ] 本轮不直接扩 full dense `xi` reference
- [ ] 本轮不要求 CI 自动 commit / 自动开 PR
- [ ] 本轮不强行把 phase diagram workflow 与 dense reference workflow 合并成一个大一统 YAML
- [ ] 本轮不新增新的 PNJL 物理公式或新 reference 字段

---

## 5. 推荐方案

### 5.1 Workflow 形态

建议新增：

- [ ] `.github/workflows/pnjl-dense-reference.yml`

入口目标：

- [ ] 通过 `workflow_dispatch` 手动触发
- [ ] 接收 dense reference 所需最小输入
- [ ] 在 CI 中运行 `scripts/pnjl/build_dense_phase_reference.jl`
- [ ] 上传精确 artifact
- [ ] 在 summary 中写清本次 reference 的配置、路径与 provenance

### 5.2 建议输入参数

第一版建议至少支持：

- [ ] `tag`
  - 例如 `dense`, `dense_xi21`, `dense_mu0_only`
- [ ] `xi_values`
  - 允许显式 csv 列表
- [ ] `xi_min`
- [ ] `xi_max`
- [ ] `xi_step`
- [ ] `T_min`
- [ ] `T_max`
- [ ] `rho_min`
- [ ] `rho_max`
- [ ] `rho_step`
- [ ] `crossover_only`
- [ ] `crossover_mu0_only`
- [ ] `crossover_n_mu`
- [ ] `crossover_mu_max`
- [ ] `overwrite`

当前建议：

- [x] 第一版以 `crossover_only=true` 为主
- [x] 先服务 `crossover_dense.csv` 扩展场景
- [x] 暂不要求第一版覆盖完整 boundary / spinodal / CEP 全链正式产物晋升

### 5.3 建议 artifact 契约

第一版至少上传：

- [ ] `data/reference/pnjl/crossover_<tag>.csv`
- [ ] `data/reference/pnjl/crossover_<tag>.meta.json`
- [ ] `data/reference/pnjl/phase_reference_<tag>_manifest.json`

如开启完整 phase pipeline，再额外上传：

- [ ] `boundary_<tag>.csv`
- [ ] `cep_<tag>.csv`
- [ ] `spinodals_<tag>.csv`

artifact 命名建议：

- [ ] `pnjl-dense-reference-<tag>`

retention 建议：

- [ ] 第一版 `30` 或 `60` 天
- [ ] 不建议默认走过长 retention，避免 artifact 膨胀

### 5.4 建议 summary 内容

至少写清：

- [ ] workflow 输入参数
- [ ] `tag`
- [ ] `xi` 覆盖
- [ ] `mu_q` 覆盖
- [ ] 生成文件列表
- [ ] manifest 路径
- [ ] 对应 git commit SHA
- [ ] 是否为 `crossover_only`

### 5.5 建议晋升流程

推荐明确分成两段：

1. CI 生成 artifact  
   - [ ] workflow 仅负责生成与上传
2. 人工审阅后晋升到仓库正式 reference  
   - [ ] 下载 artifact
   - [ ] 本地核验 `meta.json` / manifest / 行数 / spot-check
   - [ ] 再提交到仓库

当前建议：

- [x] 第一版不要在 CI 内自动改仓库文件
- [x] 正式 reference 仍通过人工审阅后提交，避免 silent promotion

---

## 6. 任务分解

### 6.1 阶段 A：workflow 设计与骨架复用

- [ ] 盘点 `pnjl-phase-diagram.yml` 中可直接复用的 step
- [ ] 提炼 dense reference 专用输入集合
- [ ] 设计新的 artifact 名称、路径白名单与 retention
- [ ] 设计 summary 模板

验收：

- [ ] 给出新 workflow 草案结构
- [ ] 明确“复用哪些 step / 重写哪些 step”

### 6.2 阶段 B：CI 执行链接线

- [ ] 新建 `pnjl-dense-reference.yml`
- [ ] 接线 `scripts/pnjl/build_dense_phase_reference.jl`
- [ ] 把输出落到稳定临时目录后上传 artifact
- [ ] 确保 `.meta.json` 与 manifest 被一并上传

验收：

- [ ] 手动触发一次 smoke 级 run 成功
- [ ] artifact 内文件集与预期一致

### 6.3 阶段 C：最小校验与治理

- [ ] 增加文件存在性与最小行数检查
- [ ] 增加 `meta.json` / csv / manifest 的一致性检查
- [ ] 在 summary 中输出关键 provenance

验收：

- [ ] CI 日志可直接判断产物是否完整
- [ ] artifact 可在不打开仓库工作树的情况下完成初步审阅

### 6.4 阶段 D：晋升到正式 reference 的操作说明

- [ ] 补充 developer-facing README / 任务单说明
- [ ] 明确“何时下载 artifact、何时提交到仓库、何时更新 sidecar/meta”

验收：

- [ ] 仓库内有清晰的“CI 生成 -> 人工审阅 -> 正式落库”说明

---

## 7. 测试与验收标准

### 7.1 最小验证

- [ ] 手动触发 `pnjl-dense-reference.yml`
- [ ] 以小规模参数完成一次 smoke run
- [ ] 成功上传 csv + meta + manifest artifact
- [ ] summary 中能看到配置与产物路径

### 7.2 收尾验证

- [ ] 至少完成一次接近真实 use-case 的 `crossover_only` run
- [ ] 下载 artifact 后可无歧义识别：
  - `git commit`
  - 生成时间
  - `xi` 覆盖
  - `mu_q` 覆盖
  - dense 含义

### 7.3 DoD

- [ ] dense reference 有独立 workflow
- [ ] 与旧 phase diagram workflow 的职责边界清晰
- [ ] artifact 契约稳定
- [ ] 晋升流程明确且不依赖本机长时间运行

---

## 8. 风险与回退

### 8.1 主要风险

- [ ] GitHub-hosted runner 时间较长，可能触发超时或 cache 不稳定
- [ ] 若直接沿用旧 workflow 路径白名单，可能上传过宽文件集，增加审阅噪音
- [ ] 若一开始就加入自动提交/自动 PR，会把 reference 晋升治理复杂化

### 8.2 回退策略

- [ ] 若完整 dense run 成本过高，先只支持 `crossover_only + mu0` smoke
- [ ] 若 GitHub-hosted runner 时长不足，先保留“CI 生成轻量版 + 本机/专机生成完整版”的双轨策略
- [ ] 若新 workflow 一时不能稳定，可先在现有 `pnjl-phase-diagram.yml` 上做临时实验分支，但不作为正式长期方案

---

## 9. 当前建议结论

- [x] 可以借 `PNJL Phase Diagram` CI 的骨架
- [x] 不建议直接把它改造成 dense reference 正式 workflow
- [x] 更推荐新增专用 workflow，并保持：
  - 旧 workflow 继续服务 phase diagram
  - 新 workflow 专门服务 dense reference artifact
- [x] dense reference 的正式落库应保持“artifact 先生成，人工审阅后再提交”的治理口径

