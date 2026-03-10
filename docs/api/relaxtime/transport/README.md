# relaxtime Transport 主题入口

本主题承接 `relaxtime` 链路中的 transport 相关 API，重点回答三类问题：

1. 如何直接计算 `eta`、`zeta`、`sigma`。
2. transport workflow 如何与 `Models` 平衡态入口衔接。
3. provider 契约如何把具体模型物理量注入 transport 积分器。

## 主题边界

本主题负责：

- `TransportCoefficients` 的直接调用方式与积分配置
- `RelaxationTime` 与 `AverageScatteringRate` 的 transport 支撑链路
- provider 契约，包括 `default_transport_provider`、`Models.transport_provider`、`Models.TransportProvider`、`Models.prepare_transport_provider`
- 与 `Models` workflow 页之间的桥接关系

本主题不负责：

- `Models.solve_gap_and_transport` 这类统一业务入口的主说明
- `bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p` 的公式与结果结构

后两类内容分别位于：

- `docs/api/models/workflows/TransportWorkflow.md`
- `docs/api/models/derived/derivatives/BulkViscosityDerivatives.md`

## 推荐阅读顺序

### 1. 面向用户入口

- `Overview.md`

适合想快速判断“我应该直接调用什么”的读者。

### 2. 职责核心

- `CoreConcepts.md`

适合需要注入自定义 provider、理解 `prefer_energy_aniso` 路径、或从 `Models` 平衡态结果进入 transport 积分层的读者。

### 3. 细节页

- `TransportCoefficients.md`
- `RelaxationTime.md`
- `AverageScatteringRate.md`
- `generated/Exports.md`

### 4. workflow 细节

- `../workflow/TransportWorkflow.md`

## 最短导航

- 只想算输运系数：看 `Overview.md` 与 `TransportCoefficients.md`
- 想从 `Models` 一键入口进入：看 `docs/api/models/workflows/TransportWorkflow.md`
- 想理解 provider 与缓存：看 `CoreConcepts.md`
- 想理解 τ 来自哪里：看 `RelaxationTime.md` 与 `AverageScatteringRate.md`
- 想看 transport 主题的完整公开接口清单：看 `generated/Exports.md`