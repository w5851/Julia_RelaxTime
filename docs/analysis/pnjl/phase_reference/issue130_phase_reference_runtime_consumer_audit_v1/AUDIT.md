# Issue #130 runtime consumer compatibility audit

## 结论

- verdict: `candidate_isolation_confirmed_requires_explicit_adapter`
- 现有 consumer 仍通过旧 `data/reference/pnjl` 文件名/路径解析；未发现 candidate sibling 的隐式路径。
- strict/derived/render 表不满足旧 boundary/CEP/spinodal/crossover schema，不能直接替换。
- 下一步只能是独立的显式 adapter/consumer contract 设计；它必须保留 strict unresolved、derived `interpolated_noncertified`、CEP bracket 和 first-order/crossover 互斥语义。

## 边界

本审计为静态 source/schema 检查，不执行 Julia include，不调用 equilibrium solver，不运行 transport。因此它确认的是路径隔离与字段契约，不是 transport 数值正确性，也不是 reference promotion。
