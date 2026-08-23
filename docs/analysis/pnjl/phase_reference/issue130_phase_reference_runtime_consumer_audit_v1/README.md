# Issue #130 phase-reference runtime consumer audit v1

这是一个 solver-free 静态兼容性审计。它核验 versioned candidate 的表结构、旧 reference 默认路径和 phase/transport/paper 入口，确认 candidate 不会被隐式消费。

当前 verdict：`candidate_isolation_confirmed_requires_explicit_adapter`。candidate 的 `runtime_consumption=false` 保持不变；candidate 表与旧 consumer schema 不同，任何 runtime 切换都必须另行设计并审核显式 adapter。

本包不修改旧 reference、不运行 PNJL、不运行 transport，也不把 strict unresolved 或 derived non-certified 行静默升级为 production certificate。
