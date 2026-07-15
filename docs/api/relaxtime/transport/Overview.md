# Transport 概览

本文档是 `relaxtime` transport 主题的面向用户入口，优先回答“大多数调用方该从哪里开始”。

## 何时从这里开始

当你符合以下任一场景时，应优先阅读本页：

- 你已经有 `quark_params`、`thermo_params` 与 `tau`，想直接计算 `eta`、`zeta`、`sigma`
- 你想比较直接调用 `TransportCoefficients` 与走 `Models` workflow 的差异
- 你想知道 provider 在 transport 链路中的位置，但暂时不需要看全部契约细节

## 首选入口

### 1. 直接 transport 计算

最直接的 transport 入口是：

- `transport_coefficients`
- `shear_viscosity`
- `electric_conductivity`
- `bulk_viscosity`

这些入口定义在 `src/relaxtime/TransportCoefficients.jl`，适合你已经准备好了：

- `quark_params`
- `thermo_params`
- `tau`
- 可选的 `bulk_coeffs`

最短示例：

```julia
using .TransportCoefficients

coeffs = transport_coefficients(
    quark_params,
    thermo_params;
    tau=tau,
    bulk_coeffs=bulk_coeffs,
)
```

### 2. 统一工作流入口

如果你还没有平衡态结果，希望先求平衡态再算 transport，优先使用：

- `Models.solve_gap_and_transport`
- `Models.solve_transport_from_equilibrium`

这两个入口的业务定位见：

- `../../models/workflows/TransportWorkflow.md`
- `../workflow/TransportWorkflow.md`

### 3. provider 入口

如果你要控制分布函数、色散关系或 species 级别的质量/化学势解析，入口是：

- `default_transport_provider()`
- `Models.transport_provider(model)`
- `Models.prepare_transport_provider(provider, equilibrium; ...)`

其详细契约见 `CoreConcepts.md`。

## 三层调用关系

transport 主题建议按三层理解：

### A. 直接积分层

由 `TransportCoefficients` 提供，负责：

- η、ζ、σ 的积分
- `TransportIntegrationConfig` 的积分配置
- provider 的消费

### B. workflow 协调层

由 `workflow/TransportWorkflow.md` 对应的实现承接，负责：

- 从平衡态构造 `quark_params` 与 `thermo_params`
- 按需计算 `tau`
- 按需计算 `bulk_coeffs`
- 将 provider 预备后传给 transport 积分层

### C. 上游导数层

由 `docs/api/models/derived/derivatives/` 承接，负责：

- `bulk_viscosity_coefficients`
- `compute_B_bracket`
- `legacy_transport_c_p`

本主题只消费这些结果，不重新承担它们的主说明。

## 何时用哪条路径

### 只想算 η / σ / ζ

直接看 `TransportCoefficients.md`，通常不需要先读 workflow 页。

### 想在同一平衡态上复用多套 transport 配置

优先走 `Models.solve_transport_from_equilibrium`，然后再看 `CoreConcepts.md` 了解 provider 如何被缓存与覆写。

### 想替换默认 PNJL 分布或质量来源

先看 `CoreConcepts.md`，再回到 `TransportCoefficients.md` 对照 provider 字段被消费的方式。

### 想理解 τ 的数值来源

继续看：

- `RelaxationTime.md`
- `AverageScatteringRate.md`

## 关键输入口径

- 温度、化学势、质量、动量统一使用 fm⁻¹
- `tau` 使用 fm
- `TransportIntegrationConfig` 只负责积分网格与权重，不负责平衡态求解
- `prefer_energy_aniso` 决定 `xi != 0` 时 RS 分布更偏向“$E_{\mathrm{dist}}$ 能量直通”还是 provider 自带的各向异性分布接口；它不改变使用 $E_{\mathrm{kin}}$ 的输运核

## 下一步

- 需要公式和参数细节：看 `TransportCoefficients.md`
- 需要 provider 字段与回退规则：看 `CoreConcepts.md`
- 需要完整工作流输入分层：看 `../workflow/TransportWorkflow.md`
