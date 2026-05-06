---
title: charged / freeze-out validation 最小契约设计
archived: true
original: docs/dev/active/2026-05-03_charged_freezeout_validation最小契约设计.md
archived_date: 2026-05-06
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# charged / freeze-out validation 最小契约设计

更新日期：2026-05-03

目的：为介子数密度后续的 charged / freeze-out 子链固定一条不会过早绑定物理结论、但又能逐步落地的验证路径。

## 0. 当前治理定位

截至当前阶段，这条基于 `Blaschke:2019col` Figure 4 right panel 的链路，**不再定性为正式 validation**，而应定性为：

- **原文所用 workflow / path semantics 的最小复现样本**
- **charged / freeze-out 计算链与生产链路的打通样本**

原因很直接：

1. 原文没有给出 right-panel 路径的完整显式参数化；
2. 当前只能构造语义接近的 proxy，而不能做 strict literature reconstruction；
3. 因此当前不应对这条链设置 point-by-point 相对误差门槛，也不应把它纳入 regression gate。

当前这条链的职责应固定为：

1. 证明 workflow 入口是可重复运行的；
2. 证明 charged / freeze-out / flavor / path-profile / meson-density 输出链已打通；
3. 为后续寻找“可严格复现的文献数据”之前，先把生产链路治理稳定。

## 1. 当前判断

### 1.1 Figure 1 是否可以由项目自行重建

可以，但要分成两部分理解：

1. **相变线 / 手征转变层**
   - 这部分原则上可以由项目现有 `T-\mu` 扫描能力自行生成；
   - 当前仓库已有：
     - `Models.run_tmu_scan`
     - phase-aware 扫描与相图相关输出
   - 因此论文里那类 phase-transition / critical-line 背景层，不需要依赖论文图像反向 digitize。

2. **freeze-out line 层**
   - 这部分不能仅靠当前主线结果自动推出；
   - 需要单独固定：
     - freeze-out 参数化来源
     - `sqrt(s_NN) -> (T, \mu_B)` 映射
     - 是否还附带 `\mu_\pi(s)`、`\mu_s(s)` 的现象学规则
   - 因此当前可以说：
     - Figure 1 的相图背景层可自绘；
     - 但真正用于文献对照的 freeze-out 轨迹仍需单独定义契约。

### 1.2 当前外部 target 的角色

当前 `data/outputs/results/relaxtime/literature/meson_density_targets/` 下这批 CSV：

- 可以作为**结果层 validation target**
- 还不能作为**路径层 regression truth**

原因是：

1. 曲线值已经有了；
2. 但曲线背后的 stitched path 还未被参数化复现。

---

## 2. 建议的分层验证顺序

### V1：结果 target 结构验证

当前就应该固定，并且已经适合进入 smoke / core validation。

验证对象：

1. target 文件存在
2. 命名符合约定
3. `x_value` 严格递增
4. `y_value >= 0`
5. `.meta.md` 存在且字段齐全
6. import report 存在

角色：

- 防止人工 digitize 结果被后续误覆盖、误命名、误清洗

### V2：相图背景层自绘验证

这是下一层，但不应与 current `K/\pi` 主线回归混在一起。

验证对象：

1. 项目自绘 phase-transition / critical-line 图层
2. freeze-out 线参数化实现
3. path 输出字段与来源说明

最小产物建议：

1. 一份 `freezeout_path.csv`
2. 一份 `phase_background.csv`
3. 一页说明：
   - 数据来源
   - 参数化公式
   - 适用文献

### V3：charged/freeze-out workflow 结果验证

只有在 workflow 真正支持以下输入后，才进入这一层：

1. charged ratio
2. freeze-out path scanning
3. `\mu_\pi`
4. `\mu_s`

届时再把当前 digitized 文献曲线接入数值比较。

当前比较建议：

1. 先从 `K^-/\pi^-`, `\mu_\pi=100 MeV` 单曲线开始
2. 再加入 `K^-/\pi^-`, `\mu_\pi=134.5 MeV`
3. 最后再处理 `K^+/\pi^+` 的 band-like 结果层

但要补一条治理边界：

- 对 `Blaschke:2019col` right-panel 当前这批样本，**V3 当前只做 workflow reproduction，不做正式数值验收**。

---

## 3. 当前不应做的事

当前不建议：

1. 把 `Blaschke:2019col` Figure 4 曲线直接当成 strict regression 真值
2. 在 freeze-out 路径未定义前，强行做 point-by-point 数值断言
3. 把 `K^+/\pi^+` band 简化成单曲线目标

这些都会让“结果 target”和“路径实现”两层混淆。

---

## 4. 当前最小落地

当前已落地：

1. WPD 导出规范化脚本
   - `scripts/analysis/relaxtime/normalize_wpd_meson_density_targets.py`
2. literature target 正式落盘
3. smoke validation
   - `tests/validation/relaxtime/test_meson_density_literature_targets_smoke.jl`

这意味着 charged / freeze-out validation 已经完成了第一层治理：

- **target 已固定**
- **结构契约已固定**
- **物理 workflow 尚未过早锁死**

---

## 5. 后续推荐顺序

1. 先设计 freeze-out 线的来源与参数化契约
2. 用现有相图/相变数据生成项目自绘 Figure 1 背景层
3. 再把 charged / freeze-out meson-density workflow 挂到统一 workflow 入口
4. 最后接入当前这批 digitized curves 做数值对照

### 5.1 solver / workflow 设计建议

从当前 `Models` 约束模式架构看，后续若加入 freeze-out 路径求解，更合理的做法是：

1. 新增一个**薄模式**或等价 path-spec：
   - 输入不是独立 `(T,\mu_B)` 二元网格，
   - 而是一维路径参数（如 `\sqrt{s_{NN}}` 或等价扫描参数）；
2. 在该模式内部先把路径参数映射到 `(T,\mu_B)`；
3. 底层仍复用 `FixedMu` 主链，而不是复制一套新的 gap solver。

也就是说，freeze-out 模式在架构上应更像：

```text
path parameter -> (T, mu_B) -> FixedMu solve
```

而不是：

```text
path parameter -> new independent solver core
```

这样做的优点是：

1. 和当前 `FixedMu`/`FixedRho`/`FixedSigma` 的 mode-contract 风格一致；
2. 不会复制 solver 内核；
3. 更容易在后续扩展：
   - actual freeze-out baseline
   - pseudo-critical proxy
   - stitched critical + constant-T path

## 6. 当前结论

当前最稳妥的治理方式是：

1. **Figure 1 背景层可以项目自绘**
2. **freeze-out 线已完成 baseline 契约化，并已具备最小 path-scan 入口**
3. **当前文献 CSV 先作为结果 target 管理**
4. **验证先做结构 smoke，再做物理值对照**

截至 2026-05-03，本仓库已新增：

- `config/physics/freezeout/default.toml`
- `config/physics/freezeout/cleymans2006.toml`
- `Models.run_freezeout_fixedmu_scan`
- `Models.FreezeoutScanConfig`
- `config/physics/meson_chemical/*.toml`
- `Models.run_freezeout_meson_density_scan`
- `scripts/relaxtime/run_freezeout_meson_density_scan.jl`

这意味着当前已经具备：

```text
sqrt(s_NN) -> freeze-out profile -> (T, mu_B) -> FixedMu solve
```

以及：

```text
sqrt(s_NN) -> freeze-out profile -> meson chemical profile -> meson-density workflow
```

这条 charged / freeze-out 子链的最小 workflow 闭环。

当前仍未完成的是：

1. more complete `\mu_s` / `\mu_Q` / stitched path 现象学回填
2. stitched / critical-proxy path strategy
3. 对更多 charged targets 的批量数值对照与分层回归准入

---

## 7. 当前新增进展：首条结果层 workflow 对照已跑通

截至 `2026-05-03`，本仓库已进一步新增一条可重复运行的 charged / freeze-out literature target 对照脚本：

- `scripts/relaxtime/run_literature_target_freezeout_validation.jl`

当前默认样本为：

- target：`Blaschke:2019col` Figure 4 right panel
- 曲线：`K^-/\pi^-`
- `\mu_\pi = 100 MeV`
- workflow 口径：
  - `freezeout profile = default`
  - `meson chemical profile = blaschke2019_mu_pi_100`
  - `regime = phase_shift_gbu`

正式输出目录：

- `data/outputs/results/relaxtime/meson_density/freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default/`

其中包含：

1. `workflow_scan.csv`
2. `comparison_vs_target.csv`
3. `README.md`

当前阶段性结论：

1. charged / freeze-out workflow 到 literature target 的结果层比较链已经闭环可运行；
2. 当前 baseline workflow 给出的 `K^-/\pi^-` 量级显著低于该文献 target；
3. 但新增多 regime 对照后，当前更明确的根因判断是：
   - `stable` 与 `strict_bw_stage2` 在同一 freeze-out baseline 下已经能靠近该 target；
   - `phase_shift_current` 与 `phase_shift_gbu` 才是当前显著压低比值的主因；
4. 因此当前 validation 的大相对误差主要是 regime / formula effect，而不是 freeze-out baseline path 自身导致；
5. `mu_s` 与 stitched / critical-proxy path 仍值得继续做，但它们现在更像二阶修正项，而不是当前这条样本的第一根因。

补充说明：

原文中 `Figure 4` 右 panel 对应的扫描路径并不是“全程 actual freeze-out line”。
本地论文文本已明确写出：

1. actual chemical freeze-out line 的结果描述较差；
2. 右 panel 使用的是 `phase transition + straight constant temperature line`；
3. 因此当前项目实现中必须把：
   - freeze-out baseline
   - higher-level path strategy
   显式分层，而不能把 Figure 4 target 简化理解成“整条 freeze-out”。

同时也要明确：

1. 原文没有给出 right-panel 路径的完整显式参数化；
2. 当前可以实现的是 `critical_line_plus_constT_proxy` 这类**语义接近原文**的 path profile；
3. 但它仍不是 strict literature reconstruction，必须把假设写在 profile 名称或 source tag 里。

因此当前最终治理结论应再收紧一句：

- **这条链当前只用于“原文 workflow / 生产链路复现”，不用于正式 validation pass/fail。**

而不是回退当前 workflow 入口或重新分叉脚本拼链。
