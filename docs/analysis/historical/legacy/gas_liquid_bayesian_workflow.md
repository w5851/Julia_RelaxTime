# Gas-Liquid Bayesian Workflow Legacy Notes

本文档记录旧 `w5851/Julia` 仓库中 `Rotation_PNJL/Gas_Liquid` 的贝叶斯优化 workflow，以便以后判断是否值得在 `Julia_RelaxTime` 主线中重做。本文只做能力沉淀与迁移边界说明，不表示这些旧脚本已经成为当前主项目的有效入口。

## 来源与状态

- 隔离区来源：`D:\Desktop\_cleanup_quarantine\2026-05-27\Julia\Rotation_PNJL`
- 远端来源：`https://github.com/w5851/Julia.git`
- 相关 dirty 目录：
  - `Rotation_PNJL/src/Gas_Liquid/Advanced_Bayesian.py`
  - `Rotation_PNJL/src/Gas_Liquid/Advanced_FindTforDiff.jl`
  - `Rotation_PNJL/examples/Gas_Liquid/`
  - `Rotation_PNJL/scripts/Gas_Liquid/`
  - `Rotation_PNJL/results/output/Gas_Liquid/`
  - `Rotation_PNJL/results/figures/Gas_Liquid/`
- 安全导出：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\Julia`

补充来源：

- 隔离区来源：`D:\Desktop\_cleanup_quarantine\2026-05-27\PNJL_Simulation_new\PNJL_Simulation`
- 远端来源：`https://github.com/w5851/PNJL_Simulation.git`
- 相关本地提交：`origin/main..main` 有 13 个 ahead commit，主题集中在 `Advanced_KappaforS.jl`、`create_kappa_objective_from_s`、Python `scikit-optimize` demo、`s_scan_kappa` 与冻结线扫描输出。
- 审计期间安全导出（已删除）：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\PNJL_Simulation_new`
- 完整工作树归档（已删除）：`D:\Desktop\_cleanup_quarantine\2026-05-27\_safety_exports\PNJL_Simulation_new_working_tree_full.zip`
- 完整归档 SHA-256：`344BFA5BA569811BBD46B4D7516C79EB98CB9ED60F09A585B50D4C71F49AA4F1`

该 dirty 工作区包含一次未完成的脚本整理：部分脚本从 `scripts/Gas_Liquid` 复制到 `examples/Gas_Liquid`，旧位置改成 moved notice；输出路径从 `output/` 迁向 `results/output/`。但整理没有完成，不能直接作为主线迁移素材。

## Workflow 意图

旧 workflow 的目标是用 Python 贝叶斯优化库驱动 Julia 物理计算，寻找一组 RMF/Walecka gas-liquid 参数，使若干实验或目标输入的涨落比值在温度反查后尽量一致。

核心变量：

- 优化参数：`(rho0, B_A, K, m_ratio, E_sym)`
- 固定实验输入：多组 `(kappa3/kappa1, kappa4/kappa2)` 与对应 `mu_B`
- 温度扫描窗口：`T_min` 到 `T_max`
- 目标函数：每组涨落比值分别反查温度后，最小化温度差平方和；无法找到温度时加惩罚项

概念流程：

1. Python 侧使用 `scikit-optimize` 的 `gp_minimize` 产生候选参数。
2. Python 通过 PyJulia 初始化 Julia 环境并加载 `Advanced_FindTforDiff.jl`。
3. Julia 侧把候选参数转换成 gas-liquid/RMF 耦合常数。
4. Julia 侧计算给定 `mu_B` 下的 `kappa3/kappa1` 与 `kappa4/kappa2` 温度反查。
5. Julia 侧返回温度差平方和作为目标函数值。
6. Python 侧保存优化历史、最优参数、收敛图和参数轨迹。

`PNJL_Simulation_new` 中的后续 ahead commit 将该思路扩展为另一条 objective：

1. 从实验表读取 `(kappa3/kappa1, kappa4/kappa2, s)`。
2. 用经验冻结线关系 `muB(s)` 和 `T(muB)` 把碰撞能量点映射到 `(T, muB)`。
3. 对每个候选参数向量 `[rho0, B_A, K, m_ratio, E_sym, a, b, c, d, e]` 同时改变 RMF 参数与冻结线参数。
4. 计算模型涨落比值，返回全部实验点上的偏差平方和。

## 旧代码组件

### Julia 目标函数层

`Advanced_FindTforDiff.jl` 是旧 workflow 的物理目标函数边界，关键函数包括：

- `find_temperature_for_kappa_ratios_with_optimization_params(...)`
- `batch_find_temperatures_with_optimization_params(...)`
- `calculate_temperature_difference_sum_of_squares(...)`
- `calculate_temperature_difference_sum_of_squares_with_weights(...)`
- `create_temperature_difference_objective(...)`
- `create_weighted_temperature_difference_objective(...)`

其中 `create_temperature_difference_objective(...)` 是最值得保留的设计：先把实验输入、`mu_B` 数组、温度窗口、扫描步长和惩罚策略封装成闭包，然后优化器只需要传入候选 `(rho0, B_A, K, m_ratio, E_sym)`。

### Python 优化层

`Advanced_Bayesian.py` 是旧 workflow 的 Python orchestration 层，主要职责包括：

- 通过 PyJulia 启动 Julia 并加载目标函数。
- 用 `scikit-optimize` 管理参数空间、采集函数和迭代历史。
- 支持从已有 CSV 继续优化。
- 输出 CSV 和可视化图。
- 曾经考虑对比 `scikit-optimize`、Optuna、Hyperopt、`bayes_opt` 等库。

旧文档中推荐 `scikit-optimize` 作为第一阶段工具，因为它在 5 参数科学优化上足够稳定，并且和 `numpy/pandas/matplotlib/scipy/sklearn` 生态配合自然。Optuna 被定位为后续并行或大规模优化备选。

`PNJL_Simulation_new/examples/Gas_Liquid/demo_python_kappa_objective.py` 是更具体的 Python demo：

- 从 `data/raw/Gas_Liquid/muB_kappa_pairs.csv` 读取列 `kappa3_1`, `kappa4_2`, `s`。
- 从 `src/Gas_Liquid/Parameter_Gas_Liquid.py` 读取参数边界。
- 通过 PyJulia 加载 `src/Gas_Liquid/Advanced_KappaforS.jl`。
- 调用 Julia 的 `create_kappa_objective_from_s(...)` 创建闭包。
- 用 `gp_minimize` 优化，并把结果写到 `results/output/Gas_Liquid/demo_kappa_bayes_result.json` 与 `demo_kappa_bayes_trace.csv`。

这个 demo 的路径推导比 `Julia/Rotation_PNJL` 中的半迁移脚本更清楚：脚本位于 `examples/Gas_Liquid`，使用 `Path(__file__).resolve().parents[2]` 定位仓库根目录。

### 碰撞能量 / 冻结线 objective

`PNJL_Simulation_new/src/Gas_Liquid/Advanced_KappaforS.jl` 中的关键公式与入口：

- `T_from_μB(μB; a=0.166, b=0.139, c=0.053)`：冻结线温度关系，输入输出按 GeV 解释。
- `μB_from_S(s; d=1.308, e=0.273)`：把 `s` 视为原公式中的 `sqrt(s)`，输出 `μB`。
- `temperature_from_s(s)`：组合得到 `T(s)`。
- `create_kappa_objective_from_s(kappa_pairs, s_values; ...)`：创建面向优化器的闭包。
- `kappa_from_s_and_optimization(s, nodes; ...)`：给定 `s` 和参数后输出 `kappa1..kappa4` 与两个比值。
- `fluctuation_deviation_from_s_and_optimization(...)`：单点偏差平方和 helper。

该分支比“温度反查温度差平方和”更直接：它不再反查温度一致性，而是直接沿冻结线预测涨落比值并拟合实验点。

### 结果与实验文件

dirty 工作区包含：

- `results/output/Gas_Liquid/forwarddiff_optimization_params_scan.csv`
- `results/figures/Gas_Liquid/individual_kappas_temperature_scan.png`
- `results/figures/Gas_Liquid/kappa_ratios_temperature_scan.png`
- `PNJL_Simulation_new/data/raw/Gas_Liquid/muB_kappa_pairs.csv`
- `PNJL_Simulation_new/results/output/Gas_Liquid/demo_kappa_bayes_result.json`
- `PNJL_Simulation_new/results/output/Gas_Liquid/demo_kappa_bayes_trace.csv`
- `PNJL_Simulation_new/results/output/Gas_Liquid/s_scan_kappa.csv`
- `PNJL_Simulation_new/results/output/Gas_Liquid/scan_tmu_for_s_*.csv`

这些结果可作为旧实验参考，但不应直接进入当前主线 baseline。它们缺少当前 `Julia_RelaxTime` 的 provenance、配置指纹和回归契约。

## 当前主项目对应能力

`Julia_RelaxTime` 已经有更规范的主线能力：

- `Models` 统一入口与 gas-liquid 变体：
  - `src/models/variants/gas_liquid/GasLiquidModel.jl`
  - `src/models/variants/gas_liquid/core/EquationSet.jl`
  - `src/models/variants/gas_liquid/core/Thermodynamics.jl`
  - `docs/api/models/variants/gas_liquid/`
- 静态涨落与组合量：
  - `chi_B`, `chi1_B`, `chi2_B`, `chi3_B`, `chi4_B`
  - `baryon_Ssigma`
  - `baryon_kappa_sigma2`

因此，旧 workflow 不应按文件迁入。若未来要恢复这类能力，应围绕当前 `Models` 入口重新设计：Python 或 Julia 优化器只能调用主线稳定 API，不能继续依赖旧 `Rotation_PNJL/src/Gas_Liquid` 的 include 驱动脚本。

## 不建议直接迁移的原因

- 旧 dirty 整理未完成：例如 `scripts/Gas_Liquid/demo_python_bayesian.py` 当前语法不完整，会触发 `IndentationError`。
- `examples/Gas_Liquid` 中部分 Python 脚本仍按旧 `scripts/Gas_Liquid` 层级推导项目根目录，迁移后路径语义可能错误。
- 旧 workflow 的 `kappa` 是静态涨落比值语境，当前主项目还同时存在 transport `kappa_*`，命名必须避免混淆。
- 旧输出缺少主项目要求的 provenance、参数指纹、测试分层和可复现入口。
- 旧参数反演逻辑与当前 `GasLiquidCoreParams` 配置层没有建立清晰映射契约。
- `PNJL_Simulation_new` 的 13 个 ahead commit 未同步到远端；它们已在本轮完整工作树归档中核验，但语义上仍是本地实验历史，不应当作远端归档仓库的稳定状态。
- `PNJL_Simulation_new` 的 dirty 状态还混有 `.vscode`、`__pycache__`、结果 CSV/PNG/JSON 与当前代码改动，需要先区分“实验结果刷新”和“可复用 API 变化”。

## `PNJL_Simulation_new` dirty 项判定

第三个 dirty 项和前两个不同：它不是单纯的旧代码残留，而是一个本地 ahead 13 个提交、且工作区继续有未完成实验改动的嵌套仓库。因此清理优先级应低于 `Transport_Julia_old_version` 这类已被主线明确吸收的旧仓库。

当前可保留的设计点：

- `create_kappa_objective_from_s(...)` 把实验碰撞能量 `s`、冻结线 `muB(s)` / `T(muB)` 和 gas-liquid 参数反演合并成一个 objective，这是后续主线重做时最值得参考的边界。
- Python `demo_python_kappa_objective.py` 展示了“Python 负责贝叶斯优化编排，Julia 负责物理目标函数”的合理分工。
- `s_scan_kappa.jl` 的结果表扩展到 `kappa1..kappa4`，说明以后如果主线化，诊断输出不应只保存两个比值。

当前不宜吸收的 dirty 改动：

- `demo_python_kappa_objective.py` 把 `if __name__ == "__main__"` 保护注释掉，导致 import 时也会直接执行 `main()` 与 `gp_minimize`，不适合作为可复用示例入口。
- `s_scan_kappa.jl` 同样把脚本保护去掉；此外成功路径返回 7 列 named tuple，但失败 fallback 仍只返回 `(s, kappa3_over_kappa1, kappa4_over_kappa2)`，存在结构不一致风险。
- 参数边界与示例参数的修改看起来是一次具体拟合实验，而不是稳定 API 变更。
- 大量 `results/output` 与 `results/figures` 改动属于实验输出刷新，不能作为当前主项目 baseline 或 regression target。

Rotation 相关 dirty 改动应单独立项审计：

- `src/Rotation/Function_Rotation.jl` 修改了 Bessel 权重公式，从旧的 `besselj(n, ptr2)` 形式改为平方 Bessel 项，并移除了 `/ (4π²)` 因子。这是物理公式级差异，不能在清理旧项目时顺手吸收。
- `examples/Rotation/susceptibility_scan.jl` 增加多初值求解，并按最大 pressure 选择分支。这个思路和当前 `Julia_RelaxTime` 中已有的 pressure-max-under-constraints 候选治理概念相近，但必须先确认 rotation 模型的物理约束、归一化和回归目标。
- 相关输出 `susceptibility_scan.csv`、`tmu_rotation.csv`、`trho_rotation.csv` 只能作为旧实验痕迹，不应直接进入主线文档或测试。

## 未来重做建议

如果以后需要把这套 workflow 主线化，建议按以下顺序做一个新任务：

1. 先写 `docs/dev/active/` 任务单，明确目标是“gas-liquid 参数反演与贝叶斯优化 workflow”，不是旧脚本搬迁。
2. 在 Julia 主线先提供稳定目标函数，例如：
   - `gas_liquid_temperature_difference_objective(...)`
   - `gas_liquid_kappa_objective_from_freezeout_s(...)`
   - 输入：实验 `kappa` 比值、`muB_MeV` 数组、温度窗口、参数候选、扫描配置。
   - 输出：标量 objective、逐点诊断表、失败原因。
3. 参数契约必须显式映射：
   - `(rho0, B_A, K, m_ratio, E_sym)` 如何转换为 `GasLiquidCoreParams`
   - 冻结线参数 `(a, b, c, d, e)` 的单位、边界和文献来源
   - 单位、边界、失败惩罚、是否启用 `delta` 通道
4. Python 层只做优化编排和可视化：
   - 推荐先保留 `scikit-optimize`，后续再评估 Optuna。
   - Python 通过固定 CLI 或 JSON/CSV 合约调用 Julia，而不是直接 eval 大段 Julia 字符串。
5. 输出必须有 provenance：
   - git commit
   - model profile
   - 参数边界
   - 实验输入
   - objective 版本
   - 失败点和惩罚统计
6. 最小测试应覆盖：
   - 参数映射函数 unit test
   - 小网格 objective smoke
   - 失败惩罚行为
   - Python/CLI contract mock 或最小 integration

## 处理建议

对已清理隔离区 `Julia`：

- 不迁移 dirty 脚本到主项目。
- 不把旧结果作为主项目 baseline。
- 原隔离区和对应 `_safety_exports` 已按此前审计结论清理；本文件仅保留 workflow 方法层记录。

对当前隔离区 `PNJL_Simulation_new`：

- 不迁移 dirty 工作区、13 个 ahead commit 或实验输出。
- `create_kappa_objective_from_s` 的冻结线 objective 思想已沉淀在本文件；未来若重做，应按当前 `Models` 合同重新设计。
- 完整工作树归档已包含 tracked dirty 内容和 untracked 文件，SHA-256 为上文所列值；审计没有发现需要长期保留的代码、结果或未完成主线依赖。
- 本轮已完成彻底清理：删除 `PNJL_Simulation_new` 原目录、轻量 safety export、完整工作树归档和审计解包目录；主项目只保留本 legacy 审计记录。
