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
