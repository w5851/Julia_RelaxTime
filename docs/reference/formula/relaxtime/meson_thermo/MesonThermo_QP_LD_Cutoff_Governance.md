# 介子热力学：QP / LD 分区公式与 cutoff 治理口径

本文档是 `MesonThermo_BU_EOS_OffShell_LD.md` 的配套第二页，专门固定：

1. quasipole (`QP`) 与 Landau damping (`LD`) 的分区公式；
2. `LD cutoff / threshold` 在当前项目中的治理语义；
3. 未来正式合同与图资产至少应输出哪些字段。

---

## 1. 分区目的

在 `Maslov & Blaschke 2023` 的主线里，`LD` 不是一个可忽略的小修正，而是会显著影响：

- `P_meson(T)`
- `P_total(T)`
- `trace anomaly`
- cutoff sensitivity

因此对当前项目，`QP / LD` 必须作为**一等公民输出**，而不是仅靠总积分结果间接包含。

---

## 2. 基础记号

统一记号：

```math
s = \omega^2 - q^2.
```

其中：

- `s > 0`：timelike 区
- `s < 0`：spacelike 区

当前项目中的介子压强主线仍建立在：

```math
\delta_M(\omega,q;T)
=
\arg \mathcal D_M(\omega+i\eta,q;T)
```

以及：

```math
P_M
=
d_M \int \frac{d^3q}{(2\pi)^3}
\int \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

---

## 3. QP / LD 的正式分区

### 3.1 Landau damping 区

文献强调的 LD 区是：

```math
s < 0
\quad \Longleftrightarrow \quad
\omega^2 - q^2 < 0.
```

在当前 `\omega > 0` 的热力学积分口径下，可以写成：

```math
0 < \omega < q.
```

因此：

```math
P_M^{\mathrm{LD}}
=
d_M \int \frac{d^3q}{(2\pi)^3}
\int_0^{q} \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

### 3.2 Quasipole / timelike 区

与之互补的 `QP` 区可写成：

```math
\omega > q.
```

因此：

```math
P_M^{\mathrm{QP}}
=
d_M \int \frac{d^3q}{(2\pi)^3}
\int_q^{\omega_{\max}} \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

如果实现层对 `\omega` 做有限截断，那么 `QP` 上限是当前数值积分使用的 `\omega_{\max}`。

### 3.3 总介子压强

正式分区后总压强写成：

```math
P_M = P_M^{\mathrm{QP}} + P_M^{\mathrm{LD}}.
```

对多通道总和：

```math
P_{\mathrm{meson}}
=
\sum_M P_M^{\mathrm{QP}}
+ \sum_M P_M^{\mathrm{LD}}.
```

---

## 4. 为什么 LD 在热力学上会被放大

文献的核心物理点之一是：尽管 `LD` 区中的相移或谱函数本身未必很大，但热权重

```math
\frac{1}{e^{\omega/T}-1}
```

在 `\omega \to 0^+` 时增强低频贡献，因此 `LD` 对 pressure 与 trace anomaly 的影响会被显著放大。

对当前项目，这意味着：

1. 不能只看 `\delta_M` 的几何大小判断 `LD` 是否重要；
2. 必须在热权重下单独评估 `LD` 的积分贡献；
3. 若 `LD` 被截断，截断治理必须显式记录。

---

## 5. cutoff / threshold 的正式治理参数

文献讨论了 `LD` 贡献对动量阈值的敏感性。对当前项目，建议把这类参数明确沉淀为治理字段，而不是散落在脚本参数里。

### 5.1 建议参数

- `qmax`
  - 总动量积分上界
- `omega_max`
  - 总频率积分上界
- `ld_cutoff`
  - `LD` 动量截断
- `ld_cutoff_mode`
  - 例如 `:match_model_lambda`、`:match_qmax`、`:fixed_multiple_of_Lambda`、`:explicit`
- `ld_threshold_mode`
  - 例如 `:spacelike_strict`、`:omega_lt_q`
- `scheme`
  - `:current` / `:gbu_reference`

### 5.2 公式层语义

若对 `LD` 只在部分动量区保留，则：

```math
P_M^{\mathrm{LD}}(\Lambda_{\mathrm{LD}})
=
d_M \int_0^{\Lambda_{\mathrm{LD}}}
\frac{d^3q}{(2\pi)^3}
\int_0^q \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

若 `QP` 部分不对动量额外截断，则：

```math
P_M^{\mathrm{QP}}(\Lambda_{\mathrm{QP}})
=
d_M \int_0^{\Lambda_{\mathrm{QP}}}
\frac{d^3q}{(2\pi)^3}
\int_q^{\omega_{\max}} \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1},
```

其中通常可取 `\Lambda_{\mathrm{QP}} = q_{\max}`。

对当前项目，若目标是优先贴近 `Maslov & Blaschke 2023` 的默认口径，则更合适的默认治理是：

```math
\Lambda_{\mathrm{LD}} = \Lambda,
```

其中 `\Lambda` 是 PNJL 模型自身使用的三动量 cutoff，而不是外层数值扫描给出的 `q_{\max}`。

因此当前项目建议：

- `ld_cutoff_mode = :match_model_lambda`
  - 语义：默认把 `\Lambda_{\mathrm{LD}}` 取为 PNJL 模型 cutoff `\Lambda`
- `ld_cutoff_mode = :match_qmax`
  - 语义：legacy 数值口径，令 `\Lambda_{\mathrm{LD}} = q_{\max}`
- `ld_cutoff_mode = :explicit`
  - 语义：由调用方显式给出 `ld_cutoff`

---

## 6. 当前项目建议的最小合同字段

如果后续把 `QP / LD` 提升到正式 CSV / regression / plot-review 资产，建议至少输出：

### 6.1 通道级

- `P_pi_qp`
- `P_pi_ld`
- `P_sigma_pi_qp`
- `P_sigma_pi_ld`

若后续扩到 strange 通道，再对称增加：

- `P_K_qp`, `P_K_ld`
- `P_sigma_K_qp`, `P_sigma_K_ld`

### 6.2 聚合级

- `P_meson_qp`
- `P_meson_ld`
- `P_meson`
- `P_total`
- `trace_anomaly`

### 6.3 治理级

- `phase_shift_variant`
- `qmax`
- `omega_max`
- `ld_cutoff`
- `ld_cutoff_mode`
- `ld_threshold_mode`

---

## 7. 当前文献对齐下的通道优先级

若目标是贴近 `Maslov & Blaschke 2023` 的物理问题，而不是只做一个任意 meson EOS 原型，那么当前项目的优先级应是：

1. `pi`（轻味赝标量）
2. `sigma_pi`（轻味标量伙伴）
3. 再考虑 `K / sigma_K`

原因是该文献的核心不是“任意两种介子相加”，而是：

- 轻味 `pion` 通道的 off-shell / LD 贡献
- 以及与之同一模型中的 `sigma` 相关通道

按当前项目命名，这里的 `sigma` 更接近：

- `sigma_pi`

而不是当前 meson thermo 第一版里已经打通的：

- `K`

因此：

- `pi + K` 不能视为与该文献“同构的正式口径”
- 贴近该文献的第一优先补项应是：
  - `sigma_pi` 的相移相关 pressure / density / workflow 支持

---

## 8. plot-review 与 regression 的最低要求

当 `QP / LD` 正式进入实现后，建议最少做以下结果资产：

1. `P_meson(T)`
2. `P_meson_qp(T)`
3. `P_meson_ld(T)`
4. `trace_anomaly(T)`
5. `ld_cutoff sensitivity`

其中 `ld_cutoff sensitivity` 至少应对两组不同 `\Lambda_{\mathrm{LD}}` 给出可对比结果。

---

## 9. 本页结论

本页给出两个必须遵守的收口结论：

1. `LD` 必须被显式建模为
   ```math
   0 < \omega < q
   ```
   的独立 pressure 贡献，而不是隐式混在总积分里；
2. 若以 `Maslov & Blaschke 2023` 为文献方向，当前项目的近同构优先通道应是
   - `pi`
   - `sigma_pi`
   
而不是把 `pi / K` 直接当作文献对应实现。

关于 `P_M` 并入系统总热力学势后，如何回到项目统一自动微分热力学主线，见：

- [../../models/shared/OmegaTotal_并入介子压强后的统一AD热力学流程.md](../../models/shared/OmegaTotal_并入介子压强后的统一AD热力学流程.md)
