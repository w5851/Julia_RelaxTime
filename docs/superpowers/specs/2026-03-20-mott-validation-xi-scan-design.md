# Mott 验证与 xi 扫描设计说明

## 1. 目标

在现有仓库已具备 Mott 判据与介子质量求解实现的前提下，新增两类能力：

1. 各向同性（`xi=0`, `mu=0`）下的 Mott 温度与介子质量验证测试，基于主工作区 Fortran 参考数据。
2. 各向异性（`xi!=0`, 暂仅 `mu=0`）下的 Mott 温度与介子质量扫描脚本，输出结构化结果用于分析与回归。

## 2. 现状确认结论（已核对）

### 2.1 公式文档已存在

- `docs/reference/formula/relaxtime/propagator/MottTransition.md`
  - 给出 Mott 判据：`M_M = M_q1 + M_q2`。
- `docs/reference/formula/relaxtime/propagator/MesonMass_RPA_Pole.md`
  - 给出 RPA 极点质量方程、混合通道表达、Mott 后宽度非零解释。

### 2.2 src 代码实现已存在

- `src/relaxtime/MottTransition.jl`
  - `mott_threshold_mass` / `mott_gap` / `is_mott_point`
  - `mott_threshold_masses` / `mott_gaps`（混合介子）
- `src/relaxtime/MesonMass.jl`
  - `meson_mass_equation` / `solve_meson_mass`
  - 覆盖 `pi/K/sigma_*` 及 `eta/eta'/sigma/sigma'`。
- `src/models/workflows/MesonMassWorkflow.jl`
  - `solve_gap_and_meson_point`，串联平衡态 -> 介子质量 -> Mott 阈值/gap。

结论：核心物理计算链条已具备，本期是“验证与扫描工程化”任务。

补充（老程序输出盘点）：

- 目前 Fortran 已有落盘以平衡态与输运系数为主，未发现稳定可复用的“介子质量/Mott 温度”直接导出文件。
- 因此“按老程序可导出字段确定验证通道”在当前快照下会退化为：先补导出流程，再做 Mott/介子质量对照。

补充（通道规则对齐）：

- 验证通道改为按“老程序理论可计算能力”确定，不再依赖“当前是否已有导出文件”。
- 按现有传播子/散射实现路径，目标通道集为：
  - `pi`, `K`, `eta`, `eta_prime`, `sigma_pi`, `sigma_K`, `sigma`, `sigma_prime`。

## 3. 范围与非范围

### 3.1 范围

- 构建参考数据映射层（Fortran -> Julia 对照格式）。
- 新增 isotropic 验证测试（Mott 温度 + 介子质量）。
- 新增 anisotropic 扫描脚本（`mu=0`，按 `xi` 网格扫描）。
- 输出 CSV 与最小元数据（参数、单位、状态标志）。

本轮对齐补充：

- 验证通道不预设固定列表，改为“以 Fortran 老程序当前可导出的介子质量字段为准”。
- 各向异性默认网格锁定：`xi ∈ [-0.4, 0.4]`, `step=0.05`。

### 3.2 非范围

- 不扩展 `mu!=0` 的各向异性扫描。
- 不重写核心 Mott/介子质量求解器。
- 不在本期解决全部不收敛点，仅记录并分级。

## 4. 设计方案（推荐）

采用“单一 workflow 入口 + 薄脚本编排 + 数据映射层”的最小侵入设计。

### 4.1 参考数据映射层

在 `scripts/analysis/` 下引入轻量映射模块/函数，职责仅包含：

- 读取 Fortran 参考数据（CSV/TSV 或可解析文本）。
- 统一字段命名（如 `T_MeV`, `m_pi_MeV`, `Tmott_pi_MeV`）。
- 做单位转换和有效性校验（NaN/Inf/缺失列）。
- 输出 Julia 测试可直接消费的数据结构。

在老程序尚未直接导出介子量的阶段，映射层需兼容两类输入：

1. 现有落盘（平衡态/输运）
2. 本期新增导出的介子量文件（mass/threshold/gap/mott_flag）

### 4.2 各向同性验证测试

优先放在 `tests/validation/`，必要时补一个 `tests/regression/` 的固定点基线。

测试逻辑：

1. 读取映射后的参考数据。
2. 调用 `solve_gap_and_meson_point(...; xi=0.0)` 生成 Julia 结果。
3. 用双阈值比较（`atol + rtol`），失败时输出通道/温度点/差值。
4. 对 Mott 温度采用“gap 过零插值”或“离散最小绝对 gap”统一口径。

Mott 温度提取策略建议（结合当前代码现状）：

- 现状：`src/relaxtime/MottTransition.jl` 已提供阈值与 gap 判据，但没有独立的统一 `T_Mott` 求根 API。
- 推荐：采用“分层混合策略”而非二选一：
  1. 主路径：先在离散温度网格上寻找 gap 符号变化区间并做线性插值（稳定、便于批量）。
  2. 精化路径：若存在可靠括号区间且用户要求更高精度，再对该区间做二分迭代（收敛条件 `|gap| <= tol`）。
  3. 退化路径：无过零时取 `|gap|` 最小点并标注 `approx`，避免误报精确 Mott 点。

### 4.3 各向异性扫描脚本

脚本放在 `scripts/analysis/`，参数包括：

- `xi` 网格（默认对称区间）
- `xi` 网格默认值：`[-0.4, 0.4]`，步长 `0.05`
- 温度网格或温度搜索区间
- 介子通道集合
- 输出路径

输出字段建议最小集：

- `xi`
- `meson`
- `T_Mott_MeV`
- `mass_at_refT_MeV`（若存在固定参考温度）
- `status`（success/fallback/fail）
- `residual`

### 4.4 异常与可追溯

扫描不应因单点失败中断：

- 单点失败写 `status=fail` 与错误摘要。
- fallback（如换初值重试）写 `status=fallback`。
- 正常点写 `status=success`。

并在同目录额外输出 run metadata（时间、commit、参数）。

## 5. 风险清单与缓解

1. **判据口径差异风险**：不同实现对 Mott 温度取值方式不同。
   - 缓解：在测试中固定判据函数并写入文档。
2. **单位差异风险**：Fortran 与 Julia 可能混用自然单位/MeV。
   - 缓解：先做单位对照表，再做误差比较。
3. **数值连续性风险**：某些 `xi` 区域不收敛。
   - 缓解：重试与状态标志并存，不静默丢点。

## 6. 验收标准

- 能运行新增验证测试并输出通道级误差诊断。
- 能运行 `mu=0` 的 `xi` 扫描并生成结构化 CSV。
- 开发文档记录清楚：输入、输出、单位、状态含义、已知限制。

## 7. 实施优先级

1. 先完成映射层与 isotropic 验证测试（确保口径可对齐）。
2. 再完成 anisotropic 扫描脚本。
3. 最后补充回归与文档闭环。

## 8. 前置依赖（新增）

在执行映射层与验证测试前，先按导出规范完成 legacy 数据补齐：

- `docs/dev/active/2026-03-20_legacy_meson_export_spec.md`

若该前置未完成，Mott/介子质量跨语言验证只能先跑 Julia 自洽链路，不能做最终 Fortran 对照验收。
