# Derivatives 主题总览

本页优先回答“如何在 `Models` 入口下使用热力学导数、质量导数与 bulk viscosity 组合导数”，而不是先展开底层实现细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 在给定 `(T_fm, mu_fm)` 下获取压力、能量、熵、数密度及其导数
- 计算夸克有效质量对 `T` / `mu` 的一阶或二阶导数
- 读取 `(∂P/∂ε)_n` 与 `(∂P/∂n)_ε` 这类 bulk viscosity 组合导数
- 计算 `bulk_viscosity_coefficients` 或体粘滞公式中的 `B` 项

## 首选公开入口

主题主入口来自 `src/models/Models.jl` 对 `src/models/derivatives/ThermoDerivatives.jl` 的聚合导出。

优先关注：

- `thermo_derivatives`
- `mass_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `bulk_viscosity_coefficients`
- `compute_B_bracket`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

T_fm = 150.0 / Models.Constants_PNJL.ħc_MeV_fm
mu_fm = 320.0 / Models.Constants_PNJL.ħc_MeV_fm

md = Models.mass_derivatives(T_fm, mu_fm; order=1)
td = Models.thermo_derivatives(T_fm, mu_fm)
bulk = Models.bulk_derivative_coeffs(T_fm, mu_fm)
cp = Models.bulk_viscosity_coefficients(T_fm, mu_fm)

d1T = Models.dP_dT(T_fm, mu_fm)
d1mu = Models.dP_dmu(T_fm, mu_fm)
```

## 输入与单位口径

- `T_fm`、`mu_fm` 使用 `fm^-1`
- `order` 当前主要用于 `mass_derivatives` 与 `dP_dT` / `dP_dmu` 的阶数控制
- PNJL derivatives 默认 `derivative_backend=:auto` 会解析为 `:taylordiff`；ThermoDerivatives 主题已下线 `derivative_backend=:forwarddiff`
- 本主题默认采用 `Models` 聚合导出接口，而不是旧 `PNJL.*` 兼容路径

## 低温积分与同支响应

`mass_derivatives`、`thermo_derivatives`、`dP_dT` 和 `dP_dmu` 接受
`thermo_quadrature_kwargs::NamedTuple=(;)` 与 `base_context=nothing`。
空 NamedTuple 使用模型的默认张量 Gauss 路线；需要解析低温费米面时，可显式选用
一维自适应积分。下例复用一个已锁定的平均场和 Jacobian：

```julia
using StaticArrays

model = Models.create_model(:PNJL)
T_fm = 60.0 / 197.327
mu_fm = 2.2
quadrature = (;
    thermo_quadrature_policy=:rs_reduced_adaptive,
    thermo_quadrature_rtol=1e-8,
    thermo_quadrature_atol=1e-11,
    thermo_quadrature_maxevals=200000,
)
context = Models.PNJLChiBTaylorDiff.build_taylor_base_context(
    T_fm, SVector(mu_fm, mu_fm, mu_fm);
    model, p_num=12, t_num=4, thermo_quadrature_kwargs=quadrature,
)
common = (; model, p_num=12, t_num=4, base_context=context)
md = Models.mass_derivatives(T_fm, mu_fm; common...)
td = Models.thermo_derivatives(T_fm, mu_fm; common...)
d4mu = Models.dP_dmu(T_fm, mu_fm; order=4, common...)
```

已有平衡态时，可以在构造 context 时传入 `base_state`；局部 polish 和既有分支检查
继续生效。context 绑定模型、点位、节点配置和积分 NamedTuple，消费者省略积分选项时
继承它的设置，显式传入不同的非空设置会抛出 `ArgumentError`。
残差、Jacobian、Taylor 递推及最终压强均使用该策略。

`dP_dmu(...; order=n)` 沿三味共同的夸克化学势求导，结果单位为 `fm^(n-4)`；
重子方向导数另有 `3^-n` 因子。质量及其导数使用传入模型的质量公式和参数。
自适应控制参数及 `xi > -1` 的适用域见 [phase 积分策略](../../phase/Overview.md#phase-热积分策略)。

上述自洽 Taylor 入口要求 `T_fm > 0`。精确零温可使用
`Models.model_pressure` 与 `Models.number_densities` 的固定态自适应路线，
在远离 `|mu|=M` 的点对解析费米海求化学势导数；五场 PNJL 的 Polyakov 方向在
零温退化。近零温宜先检查有量纲响应，再解释归一化 susceptibility 或密度趋零时的比值。

## 推荐使用顺序

多数业务场景建议按以下顺序进入：

1. 先用 `thermo_derivatives` 获取基础热力学量与一阶导数
2. 需要质量信息时再用 `mass_derivatives`
3. 只关心 bulk viscosity 常用组合量时，用 `bulk_derivative_coeffs`
4. 只有在完整体粘滞公式装配时，再进入 `bulk_viscosity_coefficients` 与 `compute_B_bracket`

## 非首页首选入口

以下接口属于下游公式接口，不应作为 derivatives 首页首选入口：

- `compute_B_bracket`
