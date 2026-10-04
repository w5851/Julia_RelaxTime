## L1 高层架构图（手动）

```mermaid
flowchart LR
  subgraph Data[数据与结果]
    data_raw[data/raw]
    data_processed[data/processed]
    data_outputs[data/outputs]
  end

  subgraph Docs[文档]
    docs_api[docs/api]
    docs_guides[docs/guides]
    docs_arch[docs/architecture]
  end

  subgraph Source[核心源码]
    src_models[src/models]
    src_utils[src/utils]
    src_integration[src/integration]
    src_simulation[src/simulation]
    src_relaxtime[src/relaxtime]
  end

  subgraph Scripts[脚本与服务]
    scripts_server[scripts/server]
    scripts_dev[scripts/dev]
    scripts_relaxtime[scripts/relaxtime]
  end

  subgraph Tests[测试套件]
    tests_unit[tests/unit]
    tests_integration[tests/integration]
    tests_regression[tests/regression]
    tests_validation[tests/validation]
    tests_baselines[tests/baselines]
  end

  subgraph Web[前端]
    web_static[web/*]
  end

  config[config/*]

  src_integration --> src_utils
  src_models --> src_integration
  src_relaxtime --> src_integration
  src_simulation --> src_utils

  scripts_server --> src_simulation
  scripts_server --> src_models
  scripts_server --> web_static
  scripts_relaxtime --> src_models
  scripts_relaxtime --> src_relaxtime
  scripts_dev --> src_models

  docs_api --> src_models
  docs_api --> src_relaxtime
  docs_guides --> scripts_server

  src_models --> data_outputs
  scripts_relaxtime --> data_outputs

  tests_unit --> src_models
  tests_unit --> src_relaxtime
  tests_integration --> src_models
  tests_regression --> src_models
  tests_regression --> tests_baselines
  tests_validation --> src_models

  config --> scripts_relaxtime
```

## L2 Models 职责与调用

本图按当前职责组织；每条实际文件加载边见自动依赖图。模型的 API/capabilities/adapter 锚点继续保留，目录完整性不再依靠 noop 文件。

```mermaid
flowchart TB
  Models["Models.jl / entrypoints.jl"] --> Factory["factory.jl"]
  Factory --> NJL["njl/：NJL 与 NJL2 物理实现"]
  Factory --> PNJL["pnjl_physics/：PNJL 与磁场变体"]
  Factory --> RPNJL["rpnjl/：RPNJL 适配"]
  Factory --> Variants["variants/：rotation 与 gas_liquid"]
  Models --> API["solver/api/SolverAPI.jl"]
  API --> Runtime["solver/runtime/：约束求解"]
  Runtime --> Spec["solver/spec/：约束与 residual"]
  Runtime --> Seeds["solver/orchestrator/：种子策略"]
  Models --> Workflows["workflow_apps/：介子与输运工作流"]
  Models --> Pipeline["workflow_engine/：流程编排"]
  Workflows --> API
  Workflows --> Transport["relaxtime/：传播子、散射与输运"]
```

`solver/compat/` 和 `solver/diff/` 保存适配与导数诊断接口；热力学导数的当前实现与适用后端见 [derivatives API](../api/models/derived/derivatives/README.md)。

## L3 关键计算步骤

下列箭头表示计算结果的流向，不是文件 include 方向。

```mermaid
flowchart LR
  Amplitude["ScatteringAmplitude"] --> Differential["DifferentialCrossSection"]
  Differential --> Total["TotalCrossSection"]
  Total --> Average["AverageScatteringRate"]
  Average --> Relaxation["RelaxationTime"]
```

统一求解由 `Models` 创建或接收模型，通过 `solver/api/` 进入对应约束求解器；工作流再消费平衡态和热力学量。非 FixedMu 联合求解与 mixed-meson 约定由各自 solver/workflow 合同维护，目录清理不改变这些语义。

回归测试在 `tests/regression/` 将当前计算与 `tests/baselines/` 中的内部基线比较，按各测试既定的 `rtol/atol` 判断；外部文献/实现对照属于 `tests/validation/`。
