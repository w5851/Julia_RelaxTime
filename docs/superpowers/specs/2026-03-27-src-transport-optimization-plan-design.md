# src 侧输运链路优化设计（仅设计文档）

## 1. 背景与目标

针对 `src/` 当前输运相关主链路（`TransportCoefficients`、`TransportWorkflow`、`ConfigLoader`）的盘点结论，目标是在**不改变物理结果与公开行为**前提下，优先提升以下工程质量：

- 可维护性：减少重复输入校验与重复默认配置构造。
- 运行稳定性：降低不必要分配与 GC 噪音。
- 演进可控性：为后续 `TransportWorkflow` 拆分与配置缓存优化保留清晰边界。

本设计文档仅输出方案与实施计划，不包含代码改动。

## 2. 约束与非目标

### 2.1 约束

- 保持 include-driven 架构与统一入口（`Models` / `src/models/entrypoints.jl`）不变。
- 不破坏既有调用方式（含 `NamedTuple` 入参、兼容 provider 行为）。
- 任何性能优化均以现有 smoke 测试与关键回归测试作为守门。
- 变更应可拆分为小 PR，便于审阅与回滚。
- 进入实现阶段时必须先创建功能分支，保持 `main` 干净（不直接在主分支实现）。

### 2.2 非目标

- 本轮不直接执行 `TransportWorkflow.jl` 大拆分（711 行拆分属于中期主题）。
- 本轮不引入新的数值模型或物理公式。
- 本轮不推进 `Main.*` 兼容桥移除，只在设计中定义后续治理接口。

## 3. 候选方案（2-3 方案对比）

### 方案 A：A+B 先行（推荐）

先落地：

1. A：`TransportCoefficients` 校验收敛（统一 helper）。
2. B：`TransportWorkflow` 默认 physics 配置常量化（去重复 Dict 构造）。

验证稳定后，再独立推进 C（`ConfigLoader` 缓存）。

优点：

- 风险最低，行为面最可控。
- PR 粒度清晰，回归定位成本低。
- 能快速获得维护性收益。

缺点：

- 配置加载分配优化收益延后。

### 方案 B：A+B+C 一次完成

在同一迭代中同时做校验收敛、常量化与配置缓存。

优点：

- 一次性释放更多性能收益。
- 减少跨 PR 的中间状态。

缺点：

- 改动面变大，定位回归更复杂。
- `ConfigLoader` 缓存需额外考虑 profile 切换与并发一致性，增加验证负担。

### 方案 C：仅结构文档先行

仅输出任务单/DoD，不做代码。

优点：

- 风险最小，便于先统一团队共识。

缺点：

- 收益延后，技术债不减少。

### 推荐结论

推荐采用**方案 A 路线（分两阶段）**：先 A+B，再 C。该路线在当前仓库“数值正确性优先”的治理原则下，投入产出比最高。

## 4. 目标架构与职责边界

### 4.1 `TransportCoefficients` 校验层（A）

新增 2-3 个内部 helper，替代分散重复校验：

- `_validate_tau_namedtuple(tau)`：校验 6 species 的存在性、有限性、非负性。
- `_validate_quark_thermo_inputs(quark_params, thermo_params; provider)`：统一质量/化学势/T/Polyakov 项校验。
- `_validate_bulk_coeffs_isentropic(coeffs)`：集中校验长度、有限性与物理警告规则。

并将 `_validate_transport_inputs(...)` 改为编排器，仅路由到这些 helper，保留现有错误语义。

### 4.2 `TransportWorkflow` 默认配置常量层（B）

在模块顶层定义：

- `const DEFAULT_PHYSICS_FALLBACK_NT = (...)`（不可变 `NamedTuple`）
- `const DEFAULT_A_BUILDER_FALLBACK = (p_nodes=..., p_max=..., cos_nodes=..., use_aniso=...)`

`_default_prefer_energy_aniso_from_toml` 与 `_default_a_builder_config_from_toml` 改为复用常量；仅在调用 `load_config` 前将 `NamedTuple` 转为局部 `Dict` 副本，避免全局可变状态污染。
若 fallback 未来包含嵌套可变值，转换阶段必须使用 `deepcopy`；否则约束 fallback 仅包含标量不可变值。

### 4.3 `ConfigLoader` 缓存层（C，第二阶段）

新增 profile 级缓存（建议模块内私有 cache + lock）：

- key 维度：`(dir, profile, inherited_hash, default_hash, profile_file_hash)`。
  - `inherited_hash`：`inherited_configs` 规范化后序列化并计算内容哈希。
  - `default_hash`：`default_config` 规范化后序列化并计算内容哈希。
  - `profile_file_hash`：`<dir>/<profile>.toml` 文件内容哈希；文件不存在时使用固定哨兵值（如 `"__missing__"`）。
- fingerprint 规范固定为：`SHA-256 + canonical JSON`。
  - canonical JSON 规则：`Dict` 键按字典序排序；序列化保留类型标签（区分 `Int/Float64/Bool/String/Nothing/Array/Dict`）；浮点采用稳定文本表示。
- canonical JSON 落地细则：
  - `dir` 先做 `abspath(normpath(dir))` 再参与 key 构造。
  - 键统一为 UTF-8 字符串并按字典序排序。
  - 浮点格式固定为 `%.17g`；`-0.0` 归一化为 `0.0`。
  - `NaN/Inf` 作为非法输入处理（直接抛错，不进入哈希）。
- 示例（示意）：`key_material = {"dir":"...","profile":"default","default":...}` -> `sha256(canonical_json(key_material))`。
- value：cache 内存储 merge 结果的不可变快照；`load_config` 对外返回防御性深拷贝，阻断调用方写污染。

并提供清理接口（仅测试/调试使用），避免长 session 中 profile 切换污染。

## 5. 数据流与调用流

### 5.1 A/B 后调用流

1. workflow 读取 TOML 默认值时命中常量 fallback。
2. workflow 构建 transport 调用参数进入 `transport_coefficients`。
3. `transport_coefficients` 在 `_validate_transport_inputs` 内分发到统一 helper。
4. 原积分/求解主流程不变。

### 5.2 C 加入后调用流

1. `load_config(...)` 先查缓存 key。
2. 命中则返回缓存副本（防止调用方写污染）。
3. 未命中执行 deep-merge，写入缓存后返回。
4. profile 切换或测试阶段可显式 reset。

## 6. 错误处理策略

- 继续使用 `ValidationUtils` 风格（`ArgumentError` + 明确参数名）。
- helper 内错误信息格式统一，减少“同类错误不同文案”。
- 对物理上可疑但可继续的情况保留 `@warn`（如负质量场景）且不静默吞错。

## 7. 测试与验证计划

### 7.1 阶段一（A+B）

- 单元 smoke：
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- 集成 smoke：
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- 回归 smoke：
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

验收标准：

- 输运相关 smoke 全绿。
- 关键接口行为不变（返回字段、异常语义、默认值语义）。

### 7.2 阶段二（C）

- 补充 `ConfigLoader` 缓存命中/失效测试：
  - 首次加载 miss，二次加载 hit。
  - profile 切换正确区分 key。
  - 配置文件内容变化可触发失效（同 profile 下 hash 变化）。
  - reset 后重新 miss。
- 并发读场景下无崩溃、无脏写。
- 记录真实性能变化（耗时/分配），不设强制门槛；在 PR 描述附命令与结果摘要。
- 固定性能对比协议：`before/after` 同机同 profile；预热 20 次、正式 1000 次；报告 `median time / allocs / bytes`。
- 固定执行入口：实现阶段新增 `scripts/perf/bench_configloader_cache.jl`，通过统一命令运行（`--threads=1`），并在 PR 附原始输出片段。

验收标准：

- 行为一致性测试通过。
- 缓存相关测试覆盖命中、失效、并发路径。
- 优化前后性能记录完整且可复现。

## 8. 风险与缓解

- 风险 1：校验收敛后错误文案变化导致快照测试波动。
  - 缓解：保留核心语义关键词与参数命名；必要时同步测试断言。
- 风险 2：配置缓存引入“被调用方写污染”。
  - 缓解：固定采用“cache 不可变快照 + 对外深拷贝返回”策略，并新增污染防护测试。
- 风险 3：锁粒度不当影响吞吐。
  - 缓解：只在读写 cache map 时加锁，deep-merge 在锁外执行。

## 9. 实施任务拆分（可直接转实现）

### PR-1（A）

- 提取并接入 `TransportCoefficients` 输入校验 helper。
- 对齐错误信息风格。
- 增补/更新对应 unit 测试。
- 在新功能分支提交（示例：`feat/src-transport-validation-unify`）。

### PR-2（B）

- 提取 `TransportWorkflow` 默认 fallback 常量。
- 将默认值常量固定为不可变 `NamedTuple`，并在 `load_config` 前局部转换 `Dict` 副本。
- 跑 smoke，确认行为不变。
- 在独立功能分支提交或基于 PR-1 分支顺序推进。

### PR-3（C，可选第二阶段）

- 为 `ConfigLoader` 增加 profile 合并缓存与 reset 接口。
- 新增缓存命中/失效/并发/污染防护测试。
- 增加性能观测脚本或测试片段，记录真实耗时与分配变化。
- 跑 smoke + 回归。
- 在独立功能分支提交，避免污染主分支。

## 10. 完成定义（DoD）

- 设计内三项优化均有可追踪任务与验证命令。
- 阶段一改动满足“低风险高回报”目标并通过 smoke。
- 阶段二在不改外部 API 的前提下提供可验证缓存收益。
- 所有新增行为均有测试覆盖与文档说明。
- fallback 门禁：若出现嵌套可变值，必须 `deepcopy` 且有引用隔离单测；若仅标量不可变值，必须有 schema 单测防止引入可变容器。
