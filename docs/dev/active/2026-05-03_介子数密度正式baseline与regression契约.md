# 介子数密度正式 baseline 与 regression 契约

更新日期：2026-05-03

目的：把介子数密度主线从阶段性 `scan/` / `analysis/` 产物，推进到可治理的正式 baseline / regression 契约层。

相关文件：

- 正式结果目录起点：
  - `data/outputs/results/relaxtime/meson_density/phase_f_stage_v1/`
- baseline 文件：
  - `tests/baselines/relaxtime/baseline_meson_density_regimes_v1.csv`
- regression 测试：
  - `tests/regression/relaxtime/test_meson_density_regimes_regression.jl`
- promotion 脚本：
  - `scripts/relaxtime/promote_meson_density_phase_f_bundle.py`

## 1. 当前契约要解决什么

当前仓库里已经有：

1. stable workflow
2. current BU workflow
3. generalized BU reference workflow
4. strict BW Stage2 workflow

但这些结果此前主要停留在：

1. `scan/`
2. `analysis/`
3. 阶段性说明页

它们还没有正式进入：

1. 可固定目录
2. baseline CSV
3. regression 入口

因此当前这一步的目标不是直接固定最终论文真值，而是先建立一层**正式工程契约**。

## 2. 正式结果目录约定

当前先引入：

- `data/outputs/results/relaxtime/meson_density/phase_f_stage_v1/`

它的角色是：

1. 把当前四条阶段性结果从 `scan/` 中提升出来；
2. 形成介子数密度主线自己的正式结果目录起点；
3. 作为后续 baseline / regression / 文档引用的统一上游。

当前仍保留 `scan/` 的作用：

1. 原始扫描输出；
2. 脚本调试与阶段性试算；
3. 中间产物。

但后续正式引用应逐步转向 `meson_density/` 目录。

## 3. 当前 baseline 的最小范围

当前 baseline 先固定为四条链在两个代表温点上的最小固定点：

1. `stable`
2. `current BU`
3. `gbu_reference`
4. `strict BW Stage2`

温点：

1. `T = 208 MeV`
2. `T = 220 MeV`

固定字段：

1. `n_pi`
2. `n_K`
3. `kpi_ratio`

当前这样设计的目的，是先锁：

1. workflow 入口没有漂移；
2. 默认参数没有漂移；
3. 四条链的数量级分工没有漂移。

## 4. 当前 regression 的解释边界

当前 regression 测试的作用是工程治理，不是物理论文定论。

它当前只回答：

1. 这四条 workflow 是否还能在相同参数下复现相同量级与分工；
2. 是否出现了无意的数值漂移或入口契约漂移。

它当前不回答：

1. 哪条链是最终物理真值；
2. 文献口径是否已经完全对齐；
3. `current BU` 与 `gbu_reference` 的最终升格关系。

## 5. 与文献验证目标的关系

后续若要引入本地参考文献库结果作为正式验证目标，建议在这套契约上继续加一层，而不是另起体系。

推荐顺序：

1. 先维持当前“项目内四条链”的 regression 契约；
2. 再新增“文献对照目标”清单；
3. 对每个文献目标明确：
   - 是 `K^+/π^+` 还是总 `K/π`
   - 是 `μ_B = 0` 温扫，还是 freeze-out / critical line
   - 是否包含 `μ_π` / `μ_s` / generalized BU 口径

## 6. 本页结论

当前已经把介子数密度主线推进到：

**有正式结果目录起点、有 baseline CSV、有 regression 测试骨架** 的状态。

下一步可以在这个基础上继续做两件事：

1. 扩充正式 `meson_density/` 目录下的生产结果与说明；
2. 评估本地文献库结果中哪些可以提升为外部验证目标。
