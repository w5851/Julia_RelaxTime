# 介子数密度 regression 兜底与人工趋势对图方案

更新日期：2026-05-06

当前状态：Phase R0/R1 已完成，Phase R2 首个 `plot_review` case 已落盘，Phase R3 已形成“双层治理：core 保留轻量 bundle consistency regression，full 新增 freeze-out `phase_shift_gbu` 全路径 rerun regression”的执行结论；用于承接 `F4.6` 在文献 direct validation 难以稳定落地时的工程兜底路线。

---

## 1. 背景与目标

当前 `介子数密度 / K/pi / BU` 主线已经具备：

1. stable
2. strict BW
3. current BU
4. generalized BU reference

四条 workflow 与对应扫描入口。

但围绕 external literature validation，当前持续存在三类困难：

1. 文献路径定义不完整
2. 当前项目主拉氏量与文献模型不完全同构
3. `BU / phase-shift` 与 `KMT / six-fermion` 两层背书尚未合并成单篇稳定 target

因此需要一条明确兜底路线：

- 当 direct literature validation 迟迟无法形成稳定 `tests/validation/` gate 时，
  将本方向先固化为：
  - 项目自产结果
  - baseline
  - regression
  - 图形人工趋势核对

本方案的目标不是替代文献验证，而是：

1. 先确保项目内 workflow 不漂移；
2. 先保留趋势可视化资产；
3. 避免因为外部 target 不稳而让主线完全失去治理抓手。

---

## 2. 范围与非范围

### 2.1 本方案覆盖

- [x] baseline 资产该放哪里
- [x] regression 测试该锁哪些量
- [x] plot 产物该如何组织
- [x] 用户人工对图时看什么
- [x] 何时允许继续只做 regression、不升格 validation

### 2.2 本方案不覆盖

- [ ] 不直接新增数值实现
- [ ] 不直接决定哪条 workflow 是最终物理真值
- [ ] 不把人工对图结果伪装成自动 validation
- [ ] 不替代后续继续寻找更强文献 target 的工作

---

## 3. 现状盘点

### 3.1 已有资产

- [x] baseline CSV：
  - `tests/baselines/relaxtime/baseline_meson_density_regimes_v1.csv`
- [x] regression 测试：
  - `tests/regression/relaxtime/test_meson_density_regimes_regression.jl`
- [x] stage-level 结果 bundle：
  - `data/outputs/results/relaxtime/meson_density/phase_f_stage_v1/`
- [x] promotion 脚本：
  - `scripts/relaxtime/promote_meson_density_phase_f_bundle.py`
- [x] charged / freeze-out reproduction 输出目录：
  - `data/outputs/results/relaxtime/meson_density/freezeout_validation/`

### 3.2 当前缺口

- [ ] 没有一套“当文献 validation 不稳时就退回 regression”的正式治理文档
- [ ] 没有专门为人工趋势对图准备的 plot 产物目录与命名契约
- [ ] 没有把“人工趋势核对”与“validation gate”做清晰切分
- [ ] 没有明确哪些结果允许进 regression、哪些只保留为分析图

---

## 4. 方案设计

### 4.1 分层原则

将结果分成三层：

1. **Regression layer**
   - 目标：锁定 workflow 行为不漂移
   - 入口：`tests/regression/`
   - 数据源：项目自产固定 baseline

2. **Plot-review layer**
   - 目标：给用户人工看趋势，不做自动 gate
   - 入口：脚本产图 + 结果目录 README
   - 数据源：项目自产扫描结果，必要时可叠文献曲线作视觉辅助

3. **Validation layer**
   - 目标：对外部 reference 做自动断言
   - 入口：`tests/validation/`
   - 前提：路径定义充分、口径清楚、主模型差异可解释

只有第 3 层需要 external truth；第 1、2 层不要求。

### 4.2 regression 兜底适用条件

当满足以下任一条件时，允许先退回 regression 方案：

- [ ] 外部 target 路径无法参数化重建
- [ ] 文献模型缺少当前项目关键项（例如六费米子项）
- [ ] 文献只给结果图，不给上游公式或路径定义
- [ ] direct validation 需要大量人工解释，难以形成稳定自动门禁

### 4.3 regression 层最小锁定对象

当前建议继续锁以下三类量：

- [ ] point-level fixed baselines
  - `n_pi`
  - `n_K`
  - `kpi_ratio`
- [ ] regime separation
  - stable
  - strict BW
  - current BU
  - generalized BU reference
- [ ] path-level smoke bundles
  - `mu_B = 0` temperature scan
  - one crossover-path scan
  - one freeze-out-path reproduction sample

说明：

- `tests/regression/relaxtime/test_meson_density_regimes_regression.jl`
  继续承担 point-level 锁定；
- 若后续扩展 path-level regression，优先新增独立 regression 文件，不污染现有 point-level 文件。

### 4.4 plot-review 层最小产物

每个需要人工趋势核对的方向，都至少产出：

- [ ] 一份原始扫描 CSV
- [ ] 一份 plot-ready 对照 CSV
- [ ] 一张 overlay 图
- [ ] 一张 residual / ratio / relative-difference 辅助图
- [ ] 一页 README，写清“这不是 validation gate”

建议目录：

- `data/outputs/results/relaxtime/meson_density/plot_review/<case_slug>/`

建议 `case_slug`：

- `mu0_regime_compare`
- `freezeout_kminus_piminus_mu_pi_100`
- `crossover_kplus_piplus_mu055`
- `phase_shift_vs_strictbw_window_208_220`

当前已落地的首个正式 case：

- `data/outputs/results/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100/`
  - 来源脚本：`scripts/relaxtime/build_meson_density_plot_review_case.py`
  - 来源样本：`data/outputs/results/relaxtime/meson_density/freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default/`
  - 当前最小产物：
    - `workflow_scan.csv`
    - `comparison_vs_target.csv`
    - `plot_review_comparison.csv`
    - `overlay_kminus_piminus_mu_pi_100.png`
    - `residual_kminus_piminus_mu_pi_100.png`
    - `README.md`

### 4.5 图形人工验收口径

人工看图时，不要求逐点数值一致；只要求以下趋势层判断：

- [ ] 单调性是否保留
- [ ] 各 regime 的相对层级是否保留
- [ ] 关键窗口的拐点/抬升/压低是否保留
- [ ] `stable / BW / current BU / gbu` 的分工是否翻转
- [ ] charged / freeze-out 路径上是否出现明显非物理跳变

不允许把“图形看起来差不多”写成 validation 通过。

### 4.6 validation 禁入条件

以下情形下，结果不得升格进 `tests/validation/`：

- [ ] 只靠人工 digitized 曲线，但路径未定义
- [ ] 文献模型与当前项目主拉氏量关键项不一致，且差异无法单独因子化
- [ ] 结果只能做 qualitative comparison
- [ ] 输出图只适合作为趋势核对

---

## 5. 可执行任务单

### Phase R0：治理文档落盘

- [x] 在活动文档中正式引用本方案
  - 验收：`2026-04-30` 与 `2026-05-06` 任务单都链接到本页

### Phase R1：plot-review 目录契约

- [x] 固定 `plot_review/` 目录结构与命名规则
  - 验收：文档写明目录、case slug、最小文件集合

### Phase R2：首个 plot-review case

- [x] 选择一个已有样本做首个正式 plot-review case
  - 推荐：`freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_*`
  - 验收：生成单独 `plot_review/<case_slug>/` 目录，含 CSV + 图 + README
  - 当前产物：
    - `data/outputs/results/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100/`

### Phase R3：path-level regression 评估

- [x] 判断是否需要新增 path-level regression 测试
  - 验收：给出“新增 / 暂不新增”的明确结论与理由
  - 当前结论：
    - 不把 freeze-out `phase_shift_gbu` 的重计算型 rerun 放入 `core`
    - 但将其正式纳入 `full`
    - 当前双层落地：
      - `core`：
        - canonical `plot_review` bundle consistency regression
        - 基线：`tests/baselines/relaxtime/baseline_meson_density_plot_review_case_v1.csv`
        - 测试：`tests/regression/relaxtime/test_meson_density_plot_review_case_regression.jl`
      - `full`：
        - freeze-out `phase_shift_gbu` 全路径 48 点 rerun regression
        - 基线：`tests/baselines/relaxtime/baseline_meson_density_freezeout_phase_shift_gbu_path_v1.csv`
        - 测试：`tests/regression/relaxtime/test_meson_density_freezeout_phase_shift_gbu_path_regression.jl`
    - `core` 当前锁定：
      - `comparison_vs_target.csv` 与 `plot_review_comparison.csv` 的逐行一致性
      - `README.md` 中的摘要指标与“非 validation gate”标注
      - overlay / residual PNG 产物存在且非空
    - `full` 当前锁定：
      - 同一 freeze-out profile / path profile / meson chemical profile / `phase_shift_gbu` 入口下的 48 点 `workflow_scan` 重算一致性
      - 逐行比较 `T_MeV`、`muB_MeV`、`n_pi`、`n_K`、`kpi_ratio`、`equilibrium_converged`

### Phase R4：validation 升格门槛

- [ ] 为未来 direct validation 升格写清触发条件
  - 验收：至少明确：
    - 路径可参数化
    - 文献模型差异已说明
    - 数值比较口径可自动化

---

## 6. 验收标准

- [x] 已有一页明确的 regression fallback 方案文档
- [x] 已定义 plot-review 目录与最小产物合同
- [x] 已明确 regression / plot-review / validation 三层边界
- [x] 已给出首个推荐 case
- [x] 已明确何时不能把结果纳入 validation

---

## 7. 风险与回退

- [ ] 风险：把 regression 当成 validation 替身
  - 回退：文档和 README 中强制标注“非 external validation gate”
- [ ] 风险：plot 资产过多、无法维护
  - 回退：每次只保留少数 canonical case
- [ ] 风险：baseline 漂移后误以为物理改进
  - 回退：baseline 更新必须附带原因、图、对比说明

---

## 8. DoD

满足以下条件即可视为本方案设计完成：

- [x] 方案文档已落到 `docs/dev/active/`
- [x] 方案已被主任务单引用
- [x] 首个推荐 plot-review case 已明确
- [x] regression 兜底与 validation 边界已写清
