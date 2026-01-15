# 约束模式模块 `PNJL.ConstraintModes`

代码位置：`src/pnjl/solver/ConstraintModes.jl`

本模块定义 PNJL 平衡求解的约束模式（方程组的“参数 θ”和“未知量 x”的维度与含义）。

## 模式

- `FixedMu()`：固定化学势，参数 `θ=[T, μ]`，未知量 5 维。
- `FixedRho(ρ/ρ0)`：固定归一化重子数密度，参数 `θ=[T]`，未知量 8 维。
- `FixedEntropy(s)`：固定熵密度，参数 `θ=[T]`，未知量 8 维。
- `FixedSigma(σ)`：固定比熵，参数 `θ=[T]`，未知量 8 维。

## 实用函数

- `state_dim(mode)`：返回未知量维度
- `param_dim(mode)`：返回参数维度
- `constraint_description(mode)`：返回人类可读描述
