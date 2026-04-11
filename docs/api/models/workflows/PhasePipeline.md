# Phase Pipeline（声明式工作流可读视图）

本文档面向“想快速理解当前相图流水线如何声明、如何执行、如何复现”的读者。

## 1. 这是不是 YAML 声明式？

当前实现不是 Snakemake/Nextflow 那种“外部 DSL（YAML/规则文件）驱动”的声明式。

当前实现是**代码内声明式（declarative in code）**：

- 声明“做什么、依赖什么、产出什么”的对象是 `PipelineSpec` 和 `PipelineStage`；
- 真正执行顺序、依赖校验、失败语义由 `PipelineRunner` 统一处理；
- CLI 只负责把参数变成 `PipelineSpec` 所需输入并触发执行。

因此：

- 不需要 YAML 才能做到声明式；
- 但可额外提供 YAML 风格“只读投影”来提升人类可读性（见第 5 节）。

## 2. 组件关系（谁负责什么）

### 2.1 声明层

- `PipelineSpec`：声明 pipeline 名称、版本、模型类型、stage 列表、I/O 契约。
- `PipelineStage`：声明单个 stage 的 `id`、`requires`、`provides`、`run!`。
- `PipelineIOContract`：声明 required inputs/outputs 与 schema 版本。

对应源码：

- `src/models/workflow/PipelineTypes.jl`

### 2.2 编排层

- `run_pipeline(...)` 在执行前做规则校验：
  - stage id 唯一；
  - `provides` 不重复（首版无 merge policy）；
  - 依赖可满足；
  - 无依赖环；
  - `spec.stages` 与 catalog 对齐。
- 运行时策略：fail-fast。
- 统一写 manifest：成功/失败都写，包含 `stage_records`。

对应源码：

- `src/models/workflow/PipelineRunner.jl`

### 2.3 目录层（phase 的标准 stage 集）

- `StageCatalog` 提供 phase pipeline 的 7 个标准 stage。
- `run_phase_pipeline_via_runner` 负责拼装 `PipelineContext`、`PipelineSpec`、`stages` 并调用 `run_pipeline`。

对应源码：

- `src/models/workflow/StageCatalog.jl`

### 2.4 入口层

- `Models.run_phase_pipeline(...)` 是薄封装，直接走 runner 路径。
- CLI 脚本将命令行参数映射到 pipeline 输入。

对应源码：

- `src/models/entrypoints.jl`
- `scripts/pnjl/calculate_phase_structure.jl`

## 3. 7-stage 执行链（当前标准）

执行顺序固定为：

1. `build_model`
2. `prepare_grid`
3. `solve_points`
4. `collect_diagnostics`
5. `analyze_phase`
6. `export_artifacts`
7. `emit_repro_manifest`

`requires/provides` 摘要（按当前实现）：

- `build_model`: requires `[:model_kind]` -> provides `[:phase_model_kind]`
- `prepare_grid`: requires `[:phase_kwargs]` -> provides `[:phase_kwargs_prepared]`
- `solve_points`: requires `[:phase_model_kind, :phase_kwargs_prepared]` -> provides `[:phase_result]`
- `collect_diagnostics`: requires `[:phase_result]` -> provides `[:phase_diagnostics]`
- `analyze_phase`: requires `[:phase_result]` -> provides `[:phase_analysis]`
- `export_artifacts`: requires `[:phase_result]` -> provides `[:phase_artifact_paths]`
- `emit_repro_manifest`: requires `[:phase_result]` -> provides `[:phase_pipeline_output]`

说明：`required_outputs` 当前约束的是 `:phase_pipeline_output`，即只有最终产出存在才算成功路径完整。

## 4. CLI -> Spec -> Runner 的映射

映射链路：

1. CLI 解析命令行到 `PhaseCliConfig`
2. 生成 `T_grid` / `rho_grid` 与 phase kwargs
3. 调用 `Models.run_phase_pipeline(model_kind; kwargs...)`
4. `run_phase_pipeline_via_runner` 构建：
   - `PipelineContext(state=Dict(:model_kind, :phase_kwargs), provenance=...)`
   - `PipelineSpec(name="phase_pipeline_runner", version="v1", stages=[...], io_contract=...)`
   - `PipelineStage[]`（来自 StageCatalog）
5. `run_pipeline` 统一执行并写 `run_manifest.json`

对使用者的含义：

- CLI 不直接编排每个子步骤；
- CLI 只声明参数意图，Runner 负责执行治理。

## 5. YAML 风格只读映射（帮助理解，不参与执行）

下面是当前实现的等价可读投影（示意）：

```yaml
pipeline:
  name: phase_pipeline_runner
  version: v1
  model_kind: PNJL
  io_contract:
    required_inputs: [model_kind, phase_kwargs]
    required_outputs: [phase_pipeline_output]

stages:
  - id: build_model
    requires: [model_kind]
    provides: [phase_model_kind]
  - id: prepare_grid
    requires: [phase_kwargs]
    provides: [phase_kwargs_prepared]
  - id: solve_points
    requires: [phase_model_kind, phase_kwargs_prepared]
    provides: [phase_result]
  - id: collect_diagnostics
    requires: [phase_result]
    provides: [phase_diagnostics]
  - id: analyze_phase
    requires: [phase_result]
    provides: [phase_analysis]
  - id: export_artifacts
    requires: [phase_result]
    provides: [phase_artifact_paths]
  - id: emit_repro_manifest
    requires: [phase_result]
    provides: [phase_pipeline_output]
```

注意：以上 YAML 是**解释层**，不是执行输入文件。

## 6. 失败语义与可复现保证

### 6.1 失败语义

- Runner 采用 fail-fast：某 stage 失败后停止后续执行。
- 失败 stage 记为 `failed`；后续未执行 stage 记为 `skipped`。
- 返回结构含 `failed_stage`、`error_kind`、`error_msg`、`completed_stages`。

### 6.2 Manifest 与可复现

`run_manifest.json` 由 runner/CLI 路径统一维护关键信息：

- pipeline 元信息（name/version/model_kind/run_id/schema）；
- `stage_records`（状态、时间、错误摘要）；
- `config_hash`、`artifact_hash`；
- 产物路径与有效配置快照。

这让“同配置、同阶段定义”具备可追溯比较基础。

## 7. 人类阅读顺序（推荐）

想快速理解“谁依赖谁、谁产出谁”，建议按这个顺序：

1. 本文（概览）
2. `src/models/workflow/StageCatalog.jl`（看 7-stage 与 requires/provides）
3. `src/models/workflow/PipelineRunner.jl`（看校验与失败语义）
4. `src/models/workflow/PipelineTypes.jl`（看契约类型）
5. `scripts/pnjl/calculate_phase_structure.jl`（看 CLI 参数如何投影到 pipeline）

这样读可以把“业务计算逻辑”和“工作流编排逻辑”分开理解。
