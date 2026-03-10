# Transport 核心边界

## Provider Contract

本文档解释 transport provider 的职责边界，以及 `relaxtime` 积分层如何消费它。

## 为什么 provider 不放在 Models 主题正文里

`transport_provider`、`TransportProvider`、`prepare_transport_provider` 由 `Models` 导出，但它们解决的问题不是“如何从 `Models` 启动业务流程”，而是：

- transport 积分层需要哪些最小物理函数
- workflow 如何把平衡态结果注入这些函数
- `xi != 0` 时 provider 路径如何选择

因此，它们的语义主归属是 transport 领域主题，而 `Models` workflow 页只负责桥接和跳转。

## 两类 provider 来源

### 1. `default_transport_provider()`

定义在 `src/relaxtime/TransportCoefficients.jl`，返回一个 `NamedTuple` provider。

它适合：

- 直接使用默认 PNJL 分布函数
- 不需要携带额外上下文对象
- 只想通过关键字覆写 `prefer_energy_aniso` 或少量行为

默认字段包括：

- `energy_from_p`
- `energy_from_p_aniso`
- `quark_distribution`
- `antiquark_distribution`
- `quark_distribution_aniso`
- `antiquark_distribution_aniso`
- `prefer_energy_aniso`

### 2. `Models.transport_provider(model)`

定义在 `src/models/transport_provider.jl`，返回 `Models.TransportProvider`。

它适合：

- provider 需要显式携带 models-side 上下文
- 后续还要通过 `prepare_transport_provider` 注入平衡态质量与化学势
- 想把 transport 积分层与具体模型实现解耦

当前已实现的模型入口是 `PNJLModel`。

## `Models.TransportProvider` 的字段职责

`Models.TransportProvider` 是一个对象型 provider，当前字段包括：

- `energy_from_p`
- `energy_from_p_aniso`
- `quark_distribution`
- `antiquark_distribution`
- `quark_distribution_aniso`
- `antiquark_distribution_aniso`
- `prefer_energy_aniso`
- `mass_for_species`
- `mu_for_species`
- `ctx`

其中应重点把握两条规则：

1. transport 层只依赖这些字段函数，不直接依赖 `ctx` 本身。
2. `ctx` 的意义是给 `mass_for_species` / `mu_for_species` 这类函数提供缓存和上下文，而不是让 transport 积分器知道 `Models` 的内部结构。

## 最小契约与可选扩展

### 最小必备能力

对于大多数 transport 计算，provider 至少应能提供：

- `energy_from_p(p, m)`
- `quark_distribution(...)`
- `antiquark_distribution(...)`

### 各向异性相关扩展

当 `xi != 0` 时，推荐提供：

- `energy_from_p_aniso(p, m, xi, costheta)`
- `quark_distribution_aniso(...)`
- `antiquark_distribution_aniso(...)`
- `prefer_energy_aniso`

### species 级别扩展

如果质量或化学势不应简单直接从 `quark_params` 读取，还可以提供：

- `mass_for_species(species, quark_params, thermo_params)`
- `mu_for_species(species, quark_params, thermo_params)`

这正是 `prepare_transport_provider` 会补齐或覆写的部分。

## `prepare_transport_provider` 在做什么

`prepare_transport_provider(provider, equilibrium; quark_params, thermo_params, masses=nothing)` 的职责是：

- 从 `equilibrium.masses` 或显式 `masses` 中提取 `(u, d, s)` 有效质量
- 从 `quark_params.μ` 中提取 `(u, d, s)` 化学势
- 把这些值写入 provider 的 `ctx`
- 重新绑定 `mass_for_species` 与 `mu_for_species`

当前最重要的语义不是“造一个全新算法对象”，而是：

- 把 workflow 手里的平衡态结果，整理成 transport 层可直接消费的 provider 契约

### 对 `TransportProvider` 的行为

如果传入的是 `Models.TransportProvider`：

- 返回同类型的新 provider
- 保留原有分布函数和能量函数
- 更新 `ctx`、`mass_for_species`、`mu_for_species`

### 对 `NamedTuple` provider 的行为

如果传入的是 `NamedTuple`：

- 自动补 `ctx`
- 若缺失则补 `energy_from_p_aniso`
- 若缺失则补 `prefer_energy_aniso=true`
- 若缺失则补 `mass_for_species` 与 `mu_for_species`

这意味着即使你不使用 `Models.TransportProvider` 类型，也可以沿着同一契约进入 transport 层。

## `prefer_energy_aniso` 的实际含义

当 `xi != 0` 时，transport 计算有两条常见路径：

1. 先计算各向异性能量，再把它直接传给分布函数。
2. 直接调用 provider 提供的 `*_distribution_aniso`。

`prefer_energy_aniso` 用来决定更偏向哪条路径：

- `true`：优先复用已经计算出的各向异性能量，减少重复 `sqrt`
- `false`：优先调用 provider 的 `*_distribution_aniso`，适合你自己实现了非平凡的各向异性分布接口

workflow、`transport_kwargs` 与 provider 默认值之间的优先级细节，见 `../workflow/TransportWorkflow.md`。

## 它在 workflow 中如何被消费

在 `src/models/workflows/TransportWorkflow.jl` 中，provider 的典型生命周期是：

1. 若调用方未显式传 provider，workflow 优先尝试 `Models.transport_provider(model)`。
2. 根据 keyword、`transport_kwargs` 和 toml 默认值，确定有效的 `prefer_energy_aniso`。
3. 若 provider 是 `Models.TransportProvider`，则调用 `Models.prepare_transport_provider(...)` 注入平衡态质量与化学势。
4. 将整理后的 provider 传给 `TransportCoefficients.transport_coefficients(...)`。

这说明 provider 主题的重点不是“谁导出了它”，而是“它在 transport 链路里承担什么契约”。

## 何时需要自定义 provider

以下场景通常值得自定义 provider：

- 你想替换默认 PNJL 分布函数实现
- 你需要不同的色散关系
- 你希望 species 质量或化学势来自额外缓存，而不是直接来自 `quark_params`
- 你想显式比较 `prefer_energy_aniso=true/false` 两条路径的数值与性能影响

如果你只是想算一组标准 transport 系数，通常不必先自定义 provider。

## 与其他页面的关系

- 用户入口：看 `Overview.md`
- 公式与积分配置：看 `TransportCoefficients.md`
- τ 与平均散射率：看 `RelaxationTime.md`、`AverageScatteringRate.md`
- `Models` 统一业务入口：看 `../../models/workflows/TransportWorkflow.md`
- 体粘滞上游导数：看 `../../models/derived/derivatives/BulkViscosityDerivatives.md`