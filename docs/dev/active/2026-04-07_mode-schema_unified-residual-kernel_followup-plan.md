# ModeSchema 与统一残差内核最小兼容方案（后续任务草案）

## 背景

- 在 PR #56 / #57 基线下，solver 主链已经收敛到 `ProblemSpec` + governed root path。
- 当前下一阶段关注点是：不同求解模式（`FixedMu`、`FixedRho`、`FixedEntropy`、`FixedSigma`、`FixedAsymmetricRho`）自变量类型与维度不一致，如何在不复制求解内核的前提下最小兼容。

## 设计结论（本轮讨论共识）

采用“统一向量残差内核 + 模式适配层”策略：

1. 仅保留一个 base 残差内核
   - 统一签名：`residual_vec!(F, x_vec, theta_vec, cfg)`。
   - 内核只关心数值维度与稳定性，不承载 mode 语义。

2. mode 差异通过 Schema 适配表达
   - 每个 mode 仅声明 `x_dim`、`theta_dim`、`pack/unpack`、`constraint_builder`。
   - 命名输入在边界层完成 `named <-> vec` 映射。

3. 模式差异降为“附加约束块”
   - `gap_core` 由统一内核计算。
   - mode 仅追加自身约束分量。
   - 总残差形式：`[gap_core; mode_constraints]`。

4. 统一问题对象，解耦求解器实现
   - `ProblemSpec` / `ImplicitProblem` 共享向量契约。
   - solver 选择（root policy / linearsolver）进入配置对象，不污染问题定义。

5. 强维度契约检查
   - 入口必须校验 `length(x)==x_dim`、`length(theta)==theta_dim`。
   - schema 与 mode 不一致时立即报错，禁止 silent fallback。

## 与现有代码的映射建议

- `src/models/solver/Conditions.jl`
  - 保留核心 gap 残差计算路径，逐步抽出统一 `gap_core_residual!`。
- `src/models/solver/ProblemSpec.jl`
  - 收口 mode 适配信息为 schema 回调，避免在 `forward_solve` 中硬编码分支。
- `src/models/solver/ImplicitAdapters.jl`
  - 继续复用 `solve_gap/gap_residual`，作为 `ImplicitProblem` 的第一层实例化入口。

## 最小实施任务（后续）

- [ ] 定义 `ModeSchema` 最小字段（维度 + pack/unpack + constraints）。
- [ ] 抽出统一 `gap_core` 入口，并保证与现有回归一致。
- [ ] 为 `FixedMu` / `FixedRho` 先实现 schema 化适配样例。
- [ ] 将 `ProblemSpec` 的 mode 分支逐步替换为 schema 驱动。
- [ ] 增加维度契约失败测试（参数长度不匹配应报 `ArgumentError`）。

## 非目标

- 不在本任务中改动物理公式。
- 不一次性重写全 mode 路径。
- 不在本任务中做大规模 API 清理。

## 验收标准（草案）

- `FixedMu` 与 `FixedRho` 在 schema 路径下结果与当前主线等价（允许数值容差内差异）。
- 统一残差内核可复用于 solver 与 implicit 接口。
- 维度错误输入能够稳定、明确地失败。
