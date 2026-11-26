# Landau 能级与积分接口

本页说明 magnetic 主题中与 Landau 能级离散化直接相关的低层公开接口。实现位于 [src/models/pnjl_physics/core/MagneticIntegrals.jl](../../../../../src/models/pnjl_physics/core/MagneticIntegrals.jl#L1)。

## 主要导出

- `QUARK_CHARGE_ABS`
- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `resolve_magnetic_nmax`
- `omega0_flavor_landau`
- `omegat_flavor_landau`
- `density_flavor_landau`

## `QUARK_CHARGE_ABS`

三味电荷绝对值常量，按 flavor 顺序组织。它属于公开常量，但更偏实现与推导辅助，而不是首页首选入口。

## `alpha_n`

Landau 简并度：

- `n = 0 -> 1`
- `n > 0 -> 2`

## `energy_landau`

定义单个 flavor、单个 Landau 层上的能谱：

$$E_{f,n} = \sqrt{2 n |q_f| eB + p_z^2 + M_f^2}$$

## `smooth_cutoff`

平滑截断因子，用于 Landau 真空项积分中的 UV 抑制。若只想做常规业务调用，不建议直接把它当作入口。

## `resolve_magnetic_nmax` 与 `resolve_nmax_from_cutoff`

`resolve_magnetic_nmax(T_fm, mu_vec, magnetic; masses=nothing)` 是生产 profile 的
统一入口。显式 `MagneticConfig.n_max` 优先；默认 `:thermal_tail` 使用
`E_tail=max(abs.(mu_vec)) + thermal_tail_factor*T_fm` 估计热尾，并在
`n_max_floor`/`n_max_cap` 预算内返回一次共享层数。`solve_magnetic_gap` 对同一个
`(T,mu,eB)` 点的全部 seed、primary 和 fallback attempt 复用该结果。

`:landau_legacy` 必须显式指定已检查收敛的 `n_max`；缺省时抛出
`ArgumentError`。热尾估计不能控制该路线的 smooth 真空求和尾部。

`resolve_nmax_from_cutoff` 仍根据质量、化学势、电荷与磁场强度给出 legacy
`:vacuum_cutoff` 策略的起点；它不是默认 MFIR 真空项的求和上限，也不替代收敛验证。
真正的主题级收敛检查仍应通过 `magnetic_nmax_convergence_report` 完成。

## `omega0_flavor_landau` / `omegat_flavor_landau`

分别给出单味 Landau 真空项与热项。它们是 thermodynamics 组装层的低级输入，不建议多数用户直接从这里开始。

`omega0_flavor_landau(...; Λ, N_c, cutoff_N, n_max, p_num, pz_max)`
采用 `fΛ²=1/(1+(p/Λ)^(2N))` 的 smooth 真空截断。模型组装层传入本模型的
`Λ_inv_fm` 和 `N_color`。此函数的 `p_num` 是整个 `[0,pz_max]` 区间的
Gauss 节点数；弱磁场下还应分别检查 Landau 层数和纵向求积阶数。

热项使用每层的 `t=asinh(pz/sqrt(M²+2n|q|eB))` 变量，在
`E=|mu|+(-12,-6,0,6,12,40)T` 分段，保留落在积分区间内的边界。
**热项的 `p_num` 是每段节点数**，与旧版整区间节点数的计算量和精度不同。
每层积分终点同时受 `pz_max` 和 `E=|mu|+40T` 限制；自动 Landau 层数
仍由配置的 `thermal_tail_factor` 控制，默认值为 30。移动节点保留在
ForwardDiff 图中，求解及求导时固定离散的 `n_max`。

质量、化学势、温度和动量单位为 `fm^-1`，`eB` 为 `fm^-2`，势密度为
`fm^-4`。热项要求有限 `T>0`、有限质量和化学势、正电荷绝对值及正
`pz_max`、非负整数 `n_max`；不满足时抛出 `ArgumentError`。精确 `T=0`
的填充海应使用独立零温公式，不能把热项置零代替。

MFIR 的 `zeta_prime_minus_one` 在 `0<x<1` 使用精确递推
`ζ′(-1,x)=ζ′(-1,x+1)-x log x`，再对右侧较大参数求积，解析保留小质量对数项。

## `density_flavor_landau`

直接积分 Polyakov 占据数之差，返回单味净密度，单位为 `fm^-3`，等于固定
均场和固定磁场下的 `-∂ΩT/∂mu`。它与热势使用相同分段和节点；颜色迹已在
Polyakov 多项式中。该接口通常与 `calculate_magnetic_number_densities` 配合使用。
