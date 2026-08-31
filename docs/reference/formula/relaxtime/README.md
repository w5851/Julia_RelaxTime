# relaxtime 公式索引

本目录沉淀 `src/relaxtime/` 相关的公式文档，目标是把“传播子 / 极化函数 / 相移 / 散射 / 输运 / 介子热力学”这些实现链路整理成可追溯的公式层说明。

## 当前主题

- `integrals/`：一圈积分与各向异性修正
- `polarization/`：极化函数与缓存
- `propagator/`：RPA 传播子、介子极点、Mott 阈值
- `meson_density/`：稳定粒子、BW、BU 相移数密度主线
- `meson_thermo/`：BU / off-shell / LD 介子压强与 EOS 主线
- `ChargedRPA_BU_ProductionRoute.md`：固定 BQS quark-only 背景上的 charged-RPA/BU candidate 公式闭合包
- `scattering/`：振幅、截面、平均散射率
- `transport/`：弛豫时间与输运系数

## 建议阅读顺序

如果目标是理解介子 EOS 主线，建议按以下顺序阅读：

1. `propagator/Propagator_传播子byPolarization.md`
2. `meson_density/MesonDensity_BU相移公式.md`
3. `meson_thermo/MesonThermo_BU_EOS_OffShell_LD.md`
4. `meson_thermo/MesonThermo_QP_LD_Cutoff_Governance.md`
5. `../models/shared/OmegaTotal_并入介子压强后的统一AD热力学流程.md`

若要审阅带电 `pi^\pm/K^\pm` 路线，先读
[ChargedRPA_BU_ProductionRoute.md](ChargedRPA_BU_ProductionRoute.md)，再回到
`couplings/`、`polarization/`、`propagator/` 和 `meson_density/` 的分层公式。

## 当前重点提醒

- `meson_density` 与 `meson_thermo` 共享同一 `phase shift` 对象，但权重不同。
- `meson_thermo` 的正式派生量口径不在 `relaxtime` 层手工定义，而要并入 `\Omega_total` 后回到 `Models` 的统一 AD 热力学主线。
