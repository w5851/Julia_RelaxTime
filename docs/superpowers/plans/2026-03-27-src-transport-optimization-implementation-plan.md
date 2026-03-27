# src Transport Optimization Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不改变现有物理结果与公开接口行为前提下，分阶段落地 src 输运链路优化（A+B，随后 C），并补齐可复现验证。

**Architecture:** 第一阶段聚焦 `TransportCoefficients` 校验收敛与 `TransportWorkflow` 默认 fallback 常量不可变化，确保低风险可维护性收益；第二阶段在 `ConfigLoader` 引入内容哈希缓存，采用“不可变快照 + 返回深拷贝”隔离策略，并通过固定性能观测协议记录真实变化。全程要求先建功能分支，禁止直接在 `main` 上实现。

**Tech Stack:** Julia 1.10+（CI 1.12.5）、include-driven 模块、`tests/unit` + `tests/integration` + `tests/regression`、`scripts/perf` 性能观测脚本。

---

## File Structure (Planned)

- Modify: `src/relaxtime/TransportCoefficients.jl`
- Modify: `src/models/workflows/TransportWorkflow.jl`
- Modify: `src/config/ConfigLoader.jl`
- Create: `scripts/perf/bench_configloader_cache.jl`
- Modify: `tests/unit/relaxtime/test_transport_coefficients.jl`
- Modify: `tests/unit/models/test_transport_workflow.jl`
- Modify: `tests/unit/config/test_config_loader.jl`

## Chunk 1: A+B（低风险高回报）

### Task 1: 创建实现分支并锁定基线

**Files:**
- Test: `tests/unit/runtests.jl`

- [x] **Step 1: 从最新主线创建功能分支**
  - Run: `git checkout -b feat/src-transport-ab-foundation`
  - 预期：当前分支不为 `main`。

- [x] **Step 2: 记录实现前基线状态**
  - Run: `git status`
  - 预期：确认仅存在预期中的工作区变更。

- [x] **Step 3: 运行 unit smoke 建立基线**
  - Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - 预期：PASS（作为后续行为不变对照）。

- [x] **Step 4: 提交“仅分支与基线记录”说明（可选）**
  - 若团队流程要求，可在执行记录文档中记一条 baseline 结论，不做代码提交。

### Task 2: 先写失败测试，覆盖 A 的校验收敛目标

**Files:**
- Modify: `tests/unit/relaxtime/test_transport_coefficients.jl`
- Test: `tests/unit/relaxtime/test_transport_coefficients.jl`

- [x] **Step 1: 新增 tau 失败测试（命名固定）**
  - 添加 `test_validate_tau_namedtuple_missing_species_throws()`。
  - 最小断言：`@test_throws ArgumentError shear_viscosity(...; tau=(u=1.0, d=1.0))`。

- [x] **Step 2: 新增 bulk_coeffs 失败测试（命名固定）**
  - 添加 `test_validate_bulk_coeffs_isentropic_length_and_finite_throws()`。
  - 最小断言：`@test_throws ArgumentError bulk_viscosity_isentropic(...; bulk_coeffs_isentropic=bad_coeffs)`。

- [x] **Step 3: 运行 selector 验证先失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_transport_coefficients.jl"; include("tests/unit/runtests.jl")'`
  - 预期：FAIL 且包含测试名 `test_validate_tau_namedtuple_missing_species_throws`。

### Task 3: 实现 A（TransportCoefficients 校验 helper 收敛）

**Files:**
- Modify: `src/relaxtime/TransportCoefficients.jl`
- Test: `tests/unit/relaxtime/test_transport_coefficients.jl`

- [x] **Step 1: 提取 `_validate_tau_namedtuple(tau)`**
  - 统一 6 species 的存在性、finite、nonnegative 检查。

- [x] **Step 2: 提取 `_validate_quark_thermo_inputs(...)`**
  - 统一 `T/Φ/Φbar` 与 `m/μ` 校验，并保留 provider 覆盖逻辑。

- [x] **Step 3: 提取 `_validate_bulk_coeffs_isentropic(coeffs)`**
  - 收敛长度、finite 校验与负质量 `@warn` 语义。

- [x] **Step 4: 将 `_validate_transport_inputs(...)` 改为编排器**
  - 仅负责分发 helper，避免重复逻辑。

- [x] **Step 5: 运行单文件测试到通过**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_transport_coefficients.jl"; include("tests/unit/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 6: Commit**
  - `git add src/relaxtime/TransportCoefficients.jl tests/unit/relaxtime/test_transport_coefficients.jl`
  - `git commit -m "refactor: unify transport input validation helpers"`

### Task 4: 先写失败测试，覆盖 B 的 fallback 不可变化目标

**Files:**
- Modify: `tests/unit/models/test_transport_workflow.jl`
- Test: `tests/unit/models/test_transport_workflow.jl`

- [x] **Step 1: 新增默认常量不污染测试（命名固定）**
  - 添加 `test_transport_workflow_default_fallback_is_not_mutated()`。
  - 最小断言：重复调用默认读取后，`prefer_energy_aniso` 一致。

- [x] **Step 2: 新增 schema 门禁测试（命名固定）**
  - 添加 `test_transport_workflow_fallback_schema_scalar_immutable_only()`。
  - 最小断言（唯一门禁）：fallback schema 仅允许标量/`NamedTuple` 等不可变类型，不允许 `Dict/Vector` 可变容器。
  - 若后续确需引入可变容器，必须先在 spec 中显式变更并新增 `deepcopy` 引用隔离测试后再放开。

- [x] **Step 3: 新增局部副本隔离测试（命名固定）**
  - 添加 `test_transport_workflow_local_dict_copy_isolation()`。
  - 最小断言：修改返回配置对象后，模块默认常量值不变。

- [x] **Step 4: 运行 selector 确认先失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_transport_workflow.jl"; include("tests/unit/runtests.jl")'`
  - 预期：FAIL 且包含上述测试名之一。

### Task 5: 实现 B（TransportWorkflow 默认常量不可变化）

**Files:**
- Modify: `src/models/workflows/TransportWorkflow.jl`
- Modify: `tests/unit/models/test_transport_workflow.jl`

- [x] **Step 1: 定义不可变 `NamedTuple` fallback 常量**
  - 引入 `DEFAULT_PHYSICS_FALLBACK_NT` 与 `DEFAULT_A_BUILDER_FALLBACK`。

- [x] **Step 2: 在 `load_config` 前执行局部 Dict 转换**
  - 严格局部副本；若出现嵌套可变值则 `deepcopy`。

- [x] **Step 3: 确保 `_default_prefer_energy_aniso_from_toml` 与 `_default_a_builder_config_from_toml` 复用常量**
  - 删除函数内重复临时 Dict 构造。

- [x] **Step 4: 跑 selector 验证通过**
  - Run: `julia --project=. --threads=4 -e 'ENV["UNIT_FILES"]="models/test_transport_workflow.jl"; include("tests/unit/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 5: Commit**
  - `git add src/models/workflows/TransportWorkflow.jl tests/unit/models/test_transport_workflow.jl`
  - `git commit -m "refactor: freeze transport workflow fallback defaults"`

### Task 6: 阶段一回归守门（A+B）

**Files:**
- Test: `tests/unit/runtests.jl`
- Test: `tests/integration/runtests.jl`
- Test: `tests/regression/runtests.jl`

- [x] **Step 1: 运行 unit smoke**
  - Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 2: 运行 integration smoke**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 3: 运行 regression smoke**
  - Run: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 4: 记录 A+B 完成结论**
  - 2026-03-27 已完成并提交：`7edd9c6 feat(relaxtime): unify transport validation and freeze workflow fallbacks`。
  - 在 PR 描述写明：接口行为、异常语义、默认值语义不变。

## Chunk 2: C（ConfigLoader 缓存与可复现性能记录）

### Task 7: 先写失败测试，定义 C 的行为契约

**Files:**
- Modify: `tests/unit/config/test_config_loader.jl`
- Test: `tests/unit/config/test_config_loader.jl`

- [x] **Step 1: 新增缓存命中/失效测试**
  - 添加 `test_config_loader_cache_hit_miss_reset()`，覆盖首次 miss、二次 hit、reset 后 miss。

- [x] **Step 2: 新增内容哈希失效测试**
  - 添加 `test_config_loader_cache_invalidate_on_profile_content_change()`。

- [x] **Step 3: 新增污染防护测试**
  - 添加 `test_config_loader_cache_returns_defensive_copy()`。

- [x] **Step 3.1: 新增 fingerprint 细则失败测试（键排序与路径归一化）**
  - 添加 `test_config_loader_fingerprint_key_sort_and_normalized_dir()`。

- [x] **Step 3.2: 新增 fingerprint 细则失败测试（浮点与特殊值）**
  - 添加 `test_config_loader_fingerprint_float_negzero_and_naninf_policy()`。
  - 最小断言：`-0.0` 归一化等价；`NaN/Inf` 触发异常。

- [x] **Step 4: 新增并发读稳定性测试**
  - 添加 `test_config_loader_cache_concurrent_reads_no_crash_no_dirty_write()`。
  - 最小断言：多线程并发读取不抛异常，且结果一致。

- [x] **Step 5: 运行 selector 确认先失败**
  - 已执行并观察到预期先失败（实现前阶段），后续实现后已转为 PASS。
  - Run: `julia --project=. --threads=4 -e 'ENV["UNIT_FILES"]="config/test_config_loader.jl"; include("tests/unit/runtests.jl")'`
  - 预期：FAIL 且包含上述测试名之一。

### Task 8: 实现 C（ConfigLoader 内容哈希缓存）

**Files:**
- Modify: `src/config/ConfigLoader.jl`
- Modify: `tests/unit/config/test_config_loader.jl`

- [x] **Step 1: 实现缓存存储与锁保护**
  - 模块内私有 cache + lock，仅 map 读写时加锁。

- [x] **Step 2: 固化 fingerprint 规范**
  - `SHA-256 + canonical JSON`；
  - `dir = abspath(normpath(dir))`；
  - 键 UTF-8 字符串按字典序；
  - 浮点 `%.17g`，`-0.0 -> 0.0`；
  - `NaN/Inf` 非法输入直接抛错。

- [x] **Step 3: 固化 value 策略**
  - cache 存不可变快照；`load_config` 返回防御性深拷贝。

- [x] **Step 4: 暴露测试可用 reset 接口**
  - 仅用于测试/调试，避免生产误用。

- [x] **Step 5: 跑单测到通过**
  - Run: `julia --project=. --threads=4 -e 'ENV["UNIT_FILES"]="config/test_config_loader.jl"; include("tests/unit/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 6: Commit**
  - `git add src/config/ConfigLoader.jl tests/unit/config/test_config_loader.jl`
  - `git commit -m "feat: add content-hash cache for config loader"`
  - 2026-03-27 已完成并提交：`558810c feat(config): add content-hash cache with defensive copies`。

### Task 9: 增加固定性能观测脚本（只记录真实变化）

**Files:**
- Create: `scripts/perf/bench_configloader_cache.jl`
- Modify: `tests/unit/config/test_config_loader.jl` (可选：脚本参数校验测试)

- [x] **Step 1: 编写脚本，固定协议参数**
  - 前置检查：若 `scripts/perf/` 不存在则先创建目录。
  - 固定输入样本：指定同一 `config` 路径与 profile（如 `config/physics/default.toml` + `default`）。
  - 固定环境记录：输出 Julia 版本、线程数、profile、样本路径。
  - `before/after` 同机同 profile；预热 20；采样 1000；报告 `median time/allocs/bytes`。
  - 固定输出字段顺序：`scenario, median_time_ns, allocs, bytes`。

- [x] **Step 2: 固定运行命令并验证可执行**
  - 运行结果（原始片段）：`before,520250,1112,56720`；`after,350300,954,44208`。
  - Run: `julia --project=. --threads=1 scripts/perf/bench_configloader_cache.jl`
  - 预期：输出 before/after 三项指标。

- [x] **Step 3: 在 PR 描述附原始输出片段**
  - 不设置硬阈值，只报告真实性能变化。

- [x] **Step 4: Commit**
  - `git add scripts/perf/bench_configloader_cache.jl`
  - `git commit -m "test: add reproducible config loader cache benchmark script"`
  - 2026-03-27 已完成并提交：`3674a37 feat(config): add reproducible config cache benchmark script`。

### Task 10: 阶段二回归守门（C）

**Files:**
- Test: `tests/unit/runtests.jl`
- Test: `tests/integration/runtests.jl`
- Test: `tests/regression/runtests.jl`

- [x] **Step 1: 运行 unit smoke**
  - Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 2: 运行 integration smoke**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 3: 运行 regression smoke**
  - Run: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
  - 预期：PASS。

- [x] **Step 4: 记录 C 完成结论**
  - 在 PR 描述总结缓存行为一致性 + 性能观测结果。
  - 结论：缓存命中/失效、防污染拷贝、并发读稳定性行为与契约一致；性能观测在固定协议下显示 `after` 相比 `before` 中位时间/分配/字节均下降。

## Chunk 3: 合并前收口

### Task 11: 文档与治理一致性收口

**Files:**
- Modify: `docs/superpowers/specs/2026-03-27-src-transport-optimization-plan-design.md`（如需补充实现偏差说明）

- [ ] **Step 1: 对照 spec 检查实现偏差**
  - 若有偏差，补一条“偏差原因+影响+回补计划”。

- [ ] **Step 2: 运行必要治理脚本（如涉及）**
  - 可选 Run: `julia --project=. scripts/dev/check_docs_consistency.jl`

- [ ] **Step 3: 最终提交（仅当有文档改动）**
  - `git add docs/superpowers/specs/2026-03-27-src-transport-optimization-plan-design.md`
  - `git commit -m "docs: align transport optimization spec with implementation details"`

### Task 12: PR 打包与交付

**Files:**
- No repository file changes required

- [x] **Step 1: 准备 PR 摘要**
  - 按 A/B/C 三段说明“为什么改、行为是否变、如何验证”。

- [x] **Step 2: 附验证证据**
  - smoke 命令与结果；缓存行为测试；性能脚本原始输出片段。

- [x] **Step 3: 明确分支策略**
  - 标注实现分支来源与禁止直改 `main` 已遵守。

- [x] **Step 4: 请求代码审阅**
  - 建议至少 1 名熟悉 `TransportWorkflow` 与 1 名熟悉 `ConfigLoader` 的 reviewer。
