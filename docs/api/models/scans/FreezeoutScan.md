# `Models.run_freezeout_fixedmu_scan`

本页从 `Models` 统一入口说明 freeze-out 路径扫描，而不是从底层模块直接调用。

## 定位

`Models.run_freezeout_fixedmu_scan` 的职责是：

- 接收一维路径参数 `\sqrt{s_{NN}}`
- 先按 freeze-out baseline 参数化映射到 `(T,\mu_B)`
- 再由 path profile 决定是否在 baseline 之上施加更高层路径策略
- 再复用现有 `FixedMu` 主链完成逐点求解

它当前是一个 **path/workflow 层薄封装**，不是新的独立 solver core。

实现转发位于 `src/models/entrypoints.jl`，底层执行位于 `src/models/scans/FreezeoutPathScan.jl`。

## 入口形态

核心关键字参数包括：

- `sqrt_s_NN_values`
- `xi_values`
- `profile_name`
- `path_profile_name`
- `output_path`
- `overwrite`
- `resume`
- `bootstrap_multiseed`
- `solver_backend`
- `auto_pnjl_backend`
- `semantic_mode`
- `selector`
- `model_kind`
- `traversal`
- `p_num`, `t_num`
- `progress_cb`

也支持结构化配置入口：

- `Models.FreezeoutScanConfig`

返回值为命名元组，包含：

- `total`
- `success`
- `failure`
- `skipped`
- `output`
- `freezeout_profile`
- `path_profile`
- `traversal`

## 当前默认物理口径

当前默认 `profile_name="default"`，对应公式页中固定的 Cleymans 型 baseline：

- `\mu_B(\sqrt{s_{NN}}) = d / (1 + e \sqrt{s_{NN}})`
- `T(\mu_B) = a - b \mu_B^2 - c \mu_B^4`

`Models.load_freezeout_profile(profile="andronic2018")` 现在返回
`AndronicFreezeoutProfile`，使用文献默认的 logistic 温度与有理式化学势。
`default` 和 `cleymans2006` 仍返回原来的 `FreezeoutParameterProfile`，原字段与
构造器保留。两者均属于 `AbstractFreezeoutParameterProfile`；扫描和 baseline
path 层按这一公共类型接收 profile，温度公式按具体类型分派。

设 `z = sqrt_s_NN_GeV`，Andronic 2018 的默认值为：

```math
T_f = \frac{0.1584}{1+\exp(2.60-\ln z/0.45)}\;\mathrm{GeV},\qquad
\mu_B = \frac{1.3075}{1+0.288z}\;\mathrm{GeV}.
```

来源是 Nature **561**, 321 (2018), Fig. 3 图注（[arXiv:1710.09425v2](https://arxiv.org/abs/1710.09425v2)，PDF p.6）。
该参数化的高能极限是 158.4 MeV，不能替换为同文 ALICE 单点拟合的 156.5 MeV。
两种配置均为化学冻结线，不包含 crossover 线或 critical-line 代理路径。

```julia
profile = Models.load_freezeout_profile(profile="andronic2018")
pt = Models.freezeout_point_MeV(profile, 7.7) # T ~= 138.428 MeV, muB ~= 406.359 MeV
points = Models.build_freezeout_scan_points([7.7, 19.6, 200.]; profile=profile)
```

映射入口接受有限、非负的碰撞能量；温度函数的化学势参数单位是 GeV。
Andronic 的 `freezeout_temperature_GeV(profile, muB_GeV)` 用同一有理式的逆映射，
定义域是 `0 <= muB_GeV <= d_GeV`，零化学势返回 `T_lim_GeV`。
零能量仅保留解析极限；实际计算的路径必须具有正能量和正温度。
`freezeout_point_MeV` 保留原有 GeV、MeV、`T_fm`、`muB_fm` 字段，后两者单位为 fm^-1。

未知命名 profile、未知 family、损坏的 TOML 或缺失的 Andronic 系数会抛出
`ArgumentError`，不会静默改用 Cleymans。旧 Cleymans 的部分系数覆盖仍可继承
默认值；新 family 必须给出完整系数，以避免混合两篇文献的映射参数。
不同路径的 q=0 外推比较入口和输出见
[化学冻结线比较脚本](../../../guides/scripts/meson_freezeout_comparison.md)。

补充的文献配置 `ihrg2019`、`frg2024_mudep`、`frg2024_muindep`、
`functional_qcd2026_mean` 均可用同一加载接口，返回 `FreezeoutParameterProfile`。
这里的多项式 family 只描述代数形式，不表示来源仍是 Cleymans。归一化曲率式
`T=T0*(1-k2*(muB/T0)^2-k4*(muB/T0)^4)` 以 GeV 为单位时对应
`a=T0,b=k2/T0,c=k4/T0^3`。API 的数学输入检查不替代文献适用能区的审查；
FRG 的 7.7 GeV 外推及守恒荷差异见
[来源与范围](../../../analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/literature_candidates_20261017.md)。

公式与参数来源见：

- [../../../reference/formula/relaxtime/meson_density/MesonDensity_Freezeout路径参数化.md](../../../reference/formula/relaxtime/meson_density/MesonDensity_Freezeout%E8%B7%AF%E5%BE%84%E5%8F%82%E6%95%B0%E5%8C%96.md)

## 与 `run_tmu_scan` 的关系

这条入口不是另一套 T-μ 求解器，而是：

```text
sqrt(s_NN) -> freeze-out profile -> (T, mu_B) -> FixedMu solve
```

因此它与 `Models.run_tmu_scan` 的关系是：

- `run_tmu_scan`：直接在显式 `(T, μ, ξ)` 网格上扫
- `run_freezeout_fixedmu_scan`：先把路径参数映射成 `(T,\mu_B)`，再沿有序路径逐点复用 `FixedMu`

## continuity / resume 契约

当前实现已复用主扫描链上的两类治理：

- continuity seed 复用
- CSV 续扫

其中输出 CSV 前三列固定为：

```text
sqrt_s_NN_GeV,muB_MeV,xi
```

这是为了与 `ScanCommon.load_completed_keys3` 的续扫键兼容。

## traversal 的作用

`traversal` 决定路径点的遍历顺序，因此也会影响 continuation 的传递方向。

当前支持：

- `:as_given`
- `:sqrts_ascending`
- `:sqrts_descending`
- `:muB_descending`
- `:muB_ascending`

默认值为 `:sqrts_ascending`。

如果目标是稳定地沿冻结线做连续性求解，推荐显式固定 traversal，而不要把输入顺序留给调用方偶然决定。

## 当前边界

当前入口现在显式拆成两层：

1. freeze-out baseline profile
2. path profile

默认 `path_profile_name="baseline_freezeout"` 时，它仍然退化为原来的 raw freeze-out baseline。

它还没有把以下更上层路径策略揉进同一入口：

- pseudo-critical proxy
- fixed `\mu/T` mapping
- stitched critical + constant-`T` path
- charged `\mu_\pi` / `\mu_s` 现象学

但当前已经预留了 path-profile 这一层，便于后续把这些策略作为显式对象继续下沉，而不是继续混写在 freeze-out baseline 公式里。

## 相关主题

- [Overview.md](Overview.md)
- [TmuScan.md](TmuScan.md)
- [TrhoScan.md](TrhoScan.md)
