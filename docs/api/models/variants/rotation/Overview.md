# Rotation 主题总览

本页优先回答“如何在 `Models` 入口下使用 rotation 变体”，而不是先展开旋转热力学公式细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 构造旋转变体模型对象
- 在给定 `(T, mu, omega)` 下完成单点求解并获取热力学量
- 在统一 `Models` 工作流下接入旋转分支

## 首选公开入口

优先关注：

- `RotationModel`
- `solve_rotation_point`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

model = Models.create_model(:Rotation)
result = Models.solve_rotation_point(0.80, 0.20; omega=0.05)
```

## 输入与单位口径

- `T_fm`、`mu_fm` 使用 `fm^-1`
- `omega` 使用与模型实现一致的自然单位约定
- 单点 workflow 输出中 `pressure`、`rho`、`entropy`、`energy` 与 `omega_potential` 均为数值量

## 非首页首选入口

`RotationThermo` 中的底层函数（例如残差构造与导数函数）属于职责核心层，不建议作为新调用方首选入口。
