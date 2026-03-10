# 状态合同与输入规范

本页集中说明 solver 主题最底层但最稳定的合同：`x_state`、`MeanFieldState` 与 `mu_vec`。

## `MeanFieldState`

`MeanFieldState` 是 `Models` 统一的平均场状态对象，字段为：

- `phi = (φ_u, φ_d, φ_s)`
- `Phi`
- `PhiBar`

设计目的：

- 明确每个分量的物理含义
- 避免直接操作裸向量时的索引错误
- 让模型层、求解层、workflow 层共享同一状态表示

## 可接受输入形式

当前统一规范支持：

- `MeanFieldState` 本身
- 长度为 3 的向量：只含 `φ_u, φ_d, φ_s`
- 长度至少为 5 的向量：`(φ_u, φ_d, φ_s, Φ, Φbar)`
- 含 `:φ` 或 `:phi` 的 `NamedTuple`

这意味着调用方不必一开始就手写 `MeanFieldState`，但主题文档默认以它作为规范表示。

## `meanfield_state`

`meanfield_state(x_state)` 的职责是把各种兼容输入规范化为 `MeanFieldState`。

适用场景：

- 你接收到了 legacy 风格裸向量
- 你要进入统一的 `Models` 状态合同链路

## `state_vector`

`state_vector` 负责把 `MeanFieldState` 统一转换为 5 维向量：

- `(φ_u, φ_d, φ_s, Φ, Φbar)`

它的主要作用不是给终端用户看，而是给求解器、导数链路和与旧向量接口兼容的内部流程使用。

## `normalize_mu_vec`

`mu_vec` 的内部标准形式是三味化学势向量。为此，`normalize_mu_vec` 提供统一规则：

- `Real`：扩展为 `(μ, μ, μ)`
- 长度为 3 的向量：逐味保留

这让 `solve_gap`、`omega`、导数链路与 workflow 可以共享相同的输入口径。

## 为什么这些合同属于主题主文档

如果没有统一的状态合同，其它入口都会各自重新定义：

- 哪些维度表示 `φ`
- `Φ` / `Φbar` 是否存在
- `mu_vec` 是单值还是三味向量

因此，`MeanFieldState`、`meanfield_state`、`state_vector`、`normalize_mu_vec` 不是边角料，而是整个 `Models` 求解主题的地基。