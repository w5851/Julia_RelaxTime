# PNJL legacy phase-reference retirement audit

## 结论

本次 verdict 为 `retirement_inconclusive`。当前仍有 `382` 个 legacy key 未被 certified candidate 精确覆盖，
并记录 `24` 个 active consumer/rollback blocker，因此本次不具备 physical deletion 条件。

## 保留与回退

legacy snapshot 继续作为 candidate 未认证/缺键行的逐键 fallback、显式 legacy rollback 和历史复现输入。
这与 RS prod_v1 的物理删除是不同任务；RS allowlist 不得用于 PNJL 文件。

## 证据边界

候选数据通过现有 Python phase-reference adapter 读取，审计过程设置 `solver_called=false`、
`reference_write=false`。本包只保存摘要、key 覆盖和引用矩阵，不复制全量曲线。

## 停止条件

若 hash/bytes 不一致、发现未知 active consumer、candidate certified coverage 下降，或 rollback 不可用，
停止 path retirement 和物理删除，恢复到当前 immutable snapshot。
