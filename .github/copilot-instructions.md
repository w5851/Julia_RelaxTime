# Julia_RelaxTime 工作区指令

- 与用户交流默认使用中文。
- 主要语言为 Julia；新增核心模块时同时补充对应单元测试。若模块属于公开稳定入口，同时补充 `docs/api/` 文档。
- 物理量命名遵循自然单位制：
  - fm⁻¹ 使用 `_inv_fm` 后缀，例如 `T_inv_fm`、`μ_inv_fm`
  - 允许使用领域内常见 Unicode 符号，如 `μ`、`Φ`
  - 更高阶单位按现有约定命名，如 `σ_fm4`、`coupling_inv_fm4`
- 模型参数配置仅从 `config/models/<model>/` 读取；共享物理常数位于 `config/physics/`。
- 测试分层遵循：`tests/unit/`、`tests/integration/`、`tests/regression/`、`tests/validation/`、`benchmark/`。
- 非测试脚本不要放入 `tests/`，分析脚本放 `scripts/analysis/`，性能探针放 `scripts/perf/`。
- smoke 测试必须保持确定性、快速、无外部依赖；性能比较应放入 `benchmark/` 而不是 smoke。
- 性能关键路径优先关注：避免热路径动态分配、避免不必要重复计算、在优化前先有 profiling 或基准依据。