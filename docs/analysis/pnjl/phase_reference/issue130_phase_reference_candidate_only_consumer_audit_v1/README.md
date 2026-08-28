# Issue #130 candidate-only consumer audit v1

这是阶段 B 的 solver-free consumer contract 审计；不运行 PNJL/equilibrium/transport，
不写 candidate 或 legacy reference。verdict：`candidate_only_contract_supported`。

- candidate calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- candidate-only rows：`11410`（全部 certified）
- candidate unresolved rows（Stage A）：`4074`
- legacy fallback keys（Stage A）：`382`
- physical deletion：`False`

`consumer_matrix.csv` 区分 candidate-only、显式 legacy rollback 和仍需迁移的 validator。
当前默认 `runtime` fallback 不变；只有后续请求键覆盖和消费者迁移完成后才可评估物理清理。
