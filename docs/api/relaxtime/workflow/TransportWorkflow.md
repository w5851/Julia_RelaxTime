# TransportWorkflow

将各向异性 PNJL 平衡求解（能隙方程 + Polyakov 环）与 RTA 输运系数计算串联。

本文档是 transport workflow 的领域细节页，重点解释内部流程、输入合同与参数优先级。

如果你只是想判断“应该从哪个 `Models` 入口开始”，优先阅读 `../../models/workflows/TransportWorkflow.md`；本页不重复承担入口选择与使用场景导航。

## 入口

### `solve_gap_and_transport`

```julia
solve_gap_and_transport(T_fm, mu_fm; xi=0.0, compute_tau=false, K_coeffs=nothing, tau=nothing, compute_bulk=true, p_num=..., t_num=..., seed_state=..., solver_kwargs=(;), tau_kwargs=(;), transport_kwargs=(;), prefer_energy_aniso=nothing)
```

- **输入单位**
  - `T_fm`, `mu_fm`：fm⁻¹
  - `xi`：无量纲
  - `tau`：fm（且按项目约定与动量 p 无关）
  - `K_coeffs`：用于内部计算 τ 的有效耦合系数
- **流程**
  1. 调用 `ThermoDerivatives.solve_equilibrium_mu(T_fm, mu_fm; xi=xi, ...)` 得到平衡解 `x_state=(φ_u,φ_d,φ_s,Φ,Φbar)` 与热力学量。
  2. 从 `x_state` 提取 `Φ, Φbar`，并计算三味有效质量 `masses`。
  3. 若 `compute_tau=true`，用 `RelaxationTime.relaxation_times` 计算平均散射率与 τ，并返回 `tau/tau_inv/rates`（需要 `K_coeffs`）。
  4. 若 `compute_bulk=true`，用 `ThermoDerivatives.bulk_derivative_coeffs` 生成体粘滞 ζ 所需的导数组合。
  5. 调用 `TransportCoefficients.transport_coefficients(quark_params, thermo_params; tau=..., bulk_coeffs=...)` 返回 `(eta, zeta, sigma)`。

- **关于 `prefer_energy_aniso`（ξ≠0 时的分布计算路径选择）**
  - 背景：当 `xi != 0` 时，输运 integrand 需要各向异性色散能量 $E_\xi(p,\cos\theta)$。
  - 默认行为：如果 provider 支持 `energy_from_p_aniso(p,m,xi,cosθ)`，并且 `prefer_energy_aniso=true`，则在 ξ≠0 时会复用已计算的 $E_\xi$，直接调用 `quark_distribution(E,...) / antiquark_distribution(E,...)`（避免 `*_distribution_aniso` 内部重复 `sqrt`）。
  - 覆写方式（二选一）：
    - 直接传 keyword：`prefer_energy_aniso=false`
    - 或者放在 `transport_kwargs`：`transport_kwargs=(; prefer_energy_aniso=false, ...)`
  - toml 默认值（方案 A）：当两种覆写都未提供时，workflow 会从 `PHYSICS_PARAM_PROFILE` 对应的 `config/physics/<profile>.toml` 的 `[transport_workflow]` 段读取 `prefer_energy_aniso` 作为默认值。
  - 典型用法：当你注入的 provider 实现了“非平凡”的 `*_distribution_aniso(p,m,...)`（不只是 RS 形式）时，可显式设置 `prefer_energy_aniso=false` 以确保走该实现。

- **返回**（NamedTuple）
  - `equilibrium`：平衡求解输出（pressure/energy/rho/.../x_state 等）
  - `quark_params`：`(m=(u,d,s), μ=(u,d,s))`
  - `thermo_params`：`(T, Φ, Φbar, ξ)`
  - `masses`：三味有效质量向量
  - `densities`：六种粒子/反粒子数密度（用于 τ 计算）
  - `tau`, `tau_inv`：弛豫时间及其倒数
  - `rates`：平均散射率（若内部计算 τ 则给出，便于复用/诊断）
  - `bulk_coeffs`：`compute_bulk=true` 时给出
  - `transport`：`(eta, zeta, sigma)`

本页把 `solve_gap_and_transport` 视为 workflow 细节入口；“它适合哪些用户场景”与“何时改用 `solve_transport_from_equilibrium`”的高层说明留给 `Models` 侧入口页。

## Day 1 输入契约冻结（v2026-02-12）

以下约定作为阶段 4 解耦期间的稳定输入合同：

### 输入分层

| 层级 | 关键输入 | 说明 |
|---|---|---|
| 平衡层 | `T_fm`, `mu_fm`, `xi`, `seed_state`, `solver_kwargs` | 只负责求平衡态与热力学状态 |
| τ 计算层（可选） | `compute_tau`, `K_coeffs`, `tau_kwargs` | `compute_tau=true` 时由 workflow 内部计算 `tau/tau_inv/rates` |
| 输运层 | `tau`（显式或内部计算）, `transport_kwargs`, `prefer_energy_aniso` | 委托给 `TransportCoefficients.transport_coefficients` |

### 参数优先级

1. 显式 keyword `prefer_energy_aniso=...`
2. `transport_kwargs.prefer_energy_aniso`
3. `config/physics/<PHYSICS_PARAM_PROFILE>.toml` 中 `[transport_workflow].prefer_energy_aniso`
4. provider 默认值（通常为 `true`）

### 输入契约细则

- `model/x_state(or equilibrium)`：当前入口对外暴露为 `solve_gap_and_transport(T_fm, mu_fm; ...)`，内部产出 `equilibrium.x_state` 后再进入输运层。
- `quark_params` 来源：由 `equilibrium` 与质量计算组装，统一传给 `TransportCoefficients`，避免 workflow 内部重复散落的 species 分支。
- `thermo_params` 来源：统一为 `(T, Φ, Φbar, ξ)`，作为输运层唯一热态输入。
- `transport_kwargs`：只承载输运积分/provider 行为相关键；不应承载平衡求解器参数。

## 示例

```julia
include("src/pnjl/workflows/TransportWorkflow.jl")
using .TransportWorkflow

T = 0.15
mu = 0.0
xi = 0.2

# 例：外部提供常数 τ（fm）
tau0 = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

res = solve_gap_and_transport(
    T,
    mu;
    xi=xi,
    tau=tau0,
    compute_tau=false,
    compute_bulk=false,
    p_num=12,
    t_num=6,
    solver_kwargs=(iterations=40,),
    transport_kwargs=(p_nodes=24, p_max=8.0,),
)

# 例：显式关闭能量直通（强制优先用 provider 的 *_distribution_aniso）
res2 = solve_gap_and_transport(
  T,
  mu;
  xi=xi,
  tau=tau0,
  compute_tau=false,
  compute_bulk=false,
  p_num=12,
  t_num=6,
  solver_kwargs=(iterations=40,),
  transport_kwargs=(p_nodes=24, p_max=8.0, prefer_energy_aniso=false),
)

@show res.thermo_params
@show res.masses
@show res.tau
@show res.transport
```

## 性能提示

- `compute_bulk=true` 会触发多次自动微分与求解（用于导数），通常明显慢于只算 η/σ。
- 扫描任务建议外层脚本自行管理 seed（用上一次点的 `equilibrium.x_state` 作为 `seed_state`）以提高收敛与速度。

## 与 Models 入口页的分工

- `docs/api/models/workflows/TransportWorkflow.md`：负责“选入口”和“理解业务闭环”
- 本页：负责“看内部流程、输入优先级与 workflow 细节”
- `../transport/CoreConcepts.md`：负责 provider 契约与 bridge 语义
