# TransportWorkflow

将各向异性 PNJL 平衡求解（能隙方程 + Polyakov 环）与 RTA 输运系数计算串联。

## 入口

### `solve_gap_and_transport`

```julia
solve_gap_and_transport(T_fm, mu_fm; xi=0.0, compute_tau=false, K_coeffs=nothing, tau=nothing, compute_bulk=true, p_num=..., t_num=..., seed_state=..., solver_kwargs=(;), tau_kwargs=(;), transport_kwargs=(;))
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

@show res.thermo_params
@show res.masses
@show res.tau
@show res.transport
```

## 性能提示

- `compute_bulk=true` 会触发多次自动微分与求解（用于导数），通常明显慢于只算 η/σ。
- 扫描任务建议外层脚本自行管理 seed（用上一次点的 `equilibrium.x_state` 作为 `seed_state`）以提高收敛与速度。
