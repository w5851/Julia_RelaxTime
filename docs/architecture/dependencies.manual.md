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
    src_root[src/ (root)]
    src_utils[src/utils]
    src_integration[src/integration]
    src_simulation[src/simulation]
    src_pnjl[src/pnjl]
    src_relaxtime[src/relaxtime]
  end

  subgraph Scripts[脚本与服务]
    scripts_server[scripts/server]
    scripts_dev[scripts/dev]
    scripts_other[scripts/...]
  end

  subgraph Web[前端]
    web_static[web/*]
  end

  config[config/*]
  tests[tests/*]

  src_utils --> src_integration
  src_integration --> src_relaxtime
  src_relaxtime --> src_pnjl
  src_simulation --> src_utils

  scripts_server --> src_simulation
  scripts_server --> src_pnjl
  scripts_server --> web_static

  scripts_other --> src_pnjl
  scripts_other --> src_relaxtime

  docs_api --> src_pnjl
  docs_api --> src_relaxtime
  docs_guides --> scripts_server

  src_pnjl --> data_outputs
  scripts_other --> data_outputs

  config --> scripts_other
  tests --> src_pnjl
  tests --> src_relaxtime
```

## L3 关键链路补充（手动）

**弛豫时间链路（RTA）**

```mermaid
flowchart LR
  ScatteringAmplitude[ScatteringAmplitude]
    --> DifferentialCrossSection[DifferentialCrossSection]
    --> TotalCrossSection[TotalCrossSection]
    --> AverageScatteringRate[AverageScatteringRate]
    --> RelaxationTime[RelaxationTime]
```

**PNJL 求解链路**

```mermaid
flowchart LR
  SeedStrategies[SeedStrategies]
    --> Solver[Solver]
    --> Scan[Scan (Tmu/Trho/DualBranch)]
```
