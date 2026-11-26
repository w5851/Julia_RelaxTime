# κ/λ 固定点回归与旧目标偏差

2026-11-05，nightly run [41996062882](https://github.com/w5851/Julia_RelaxTime/actions/runs/41996062882)
在 T=210 MeV、μB=0 的 κBB/κQQ 超过旧 5% 容差。旧目标的 `source`
为 `julia_transport_kappa_lambda_fixedpoints_v1`，是内部 Julia 导出，不是独立外部参考。
旧 CSV 原字节保留在 `tests/validation/data/targets/relaxtime/legacy/transport/`。
新消费者位于 `tests/regression/relaxtime/test_kappa_lambda_thermal_weighted_regression.jl`。

## 同点归因

在 quote-remote、Julia 1.12.5 上复用同一平衡态、数密度、A、有效耦合和末端积分，
只改变热平均配置。旧两个采样模块从 `dffe4070eb5c0788e9c8a1e227cf40978b51ce3e`
原样载入独立命名空间，散射核和其余源码保持当前版本。所有配置均使用
Λ=3.0522942112582996 fm⁻¹ 的散射支持 cutoff 和 `cm` 相对速度。
旧 `RelaxationTime` 也显式传入该 cutoff；不是截断范围或速度定义变化造成的偏差。

| T=210 MeV、μB=0 的配置 | κBB | κQQ |
| --- | ---: | ---: |
| 旧目标 | 0.05481230045 | 0.11609896464 |
| 旧 p20/c4/φ8/t6、60 点 PCHIP | 0.05481229863 | 0.11609896010 |
| 仅把旧表增至240点 | 0.05125703638 | 0.10787935686 |
| 仅把旧表增至960点 | 0.05131154221 | 0.10797634982 |
| thermal-weighted、沿用旧外层与t6 | 0.05100161832 | 0.10731512740 |
| thermal-weighted、旧外层与t24 | 0.05099495109 | 0.10730349261 |
| 当前 p48/c36/φ36/t24 | 0.05057040384 | 0.10627614855 |
| 加密 p72/c40/φ48/t48，插值 rtol=2.5e-4 | 0.05057183709 | 0.10628159735 |

旧默认复现原目标的相对差不超过4e-8；仅加密散射表就消除了大部分漂移。
旧960点表与新算法沿用旧网格仍有差异，故不把“960点”称为收敛参考。
新默认相对旧目标的变化为 −7.73895% / −8.46073%；同点加密后的变化
仅为 −0.002834% / −0.005127%。这支持旧散射表未充分收敛是主要原因，
外层网格变化也有贡献；t6→t24 在此点的影响较小。

## v2 的数值契约

[CSV](baseline_kappa_lambda_thermal_weighted_v2.csv)覆盖原四点的全部14个κ/λ值，
保存当前值与加密值。[元数据](baseline_kappa_lambda_thermal_weighted_v2.json)
记录每项旧/新差、源码、配置、哈希及32 KB诊断归档的恢复方式。未复制完整源码或结果包。

- 点位为 T=210/250 MeV × μB=0/100 MeV，ξ=0；λ只在μB=100取值，μB=0不适用。
- 继续使用既有 validation helper 的低阶背景 p12/c6、残差上限1e-4、
  传播子A的p16/c4与[0,20]区间、有效耦合A的p64与[0,10]区间，
  以及末端p64/c8与[0,10]区间。helper显式冻结热平均设置，避免生产默认变化静默改写本契约。
- 当前与加密配置共8次四点计算；另有5次T210/μB0受控比较。全部14项的最大
  加密相对差为0.019539%，对应T250/μB0的κSS。
- 回归当前值使用 `rtol=1e-6, atol=1e-8`；同时要求与记录的加密值相差不超过0.1%。
  前者保护代码结果稳定性，后者给出高于实测0.019539%的精度检查余量，仍比旧5%严格。
- 加密只核对热平均外层、截面t积分及插值。没有独立改变背景、A或末端积分，
  因此本表不是完整物理误差认证，也不替代v15生产配置及1867点的既有资格。

运行：

```powershell
$env:REGRESSION_FILES='relaxtime/test_kappa_lambda_thermal_weighted_regression.jl'
julia --project=. tests/regression/runtests.jl
```

本次物理源码、模型参数、混合介子治理、non-fixedmu联合求解和正式产物保持原值；
变化是测试配置与内部基线消费位置。外部文献及Fortran验证继续沿用各自已有目标。

## v3：`#365` 物理修正后的数值承诺变更（2026-11-26）

`#365`（`54902ee0`）修正了 `src/relaxtime/RelaxationTime.jl` 中
`_species_omega_s` 的轻夸克伙伴通道：`2.0*n.ubar*w.usbar` 改为
`2.0*n.ubar*w.subar`。s 夸克的碰撞伙伴应为 ū，原写法取的是 u–s̄ 通道。
该修正没有同步更新本基线，nightly 自 2026-11-24 起在
`full (regression)` 报失配。

### 判断依据

- **配对正确**：2+1 味同位旋对称下 s 夸克的轻夸克伙伴为 u（`w.us`）与 ū
  （`w.subar`），2.0 因子按简并对加；`_species_omega_sbar` 的
  `2.0*n.u*w.usbar` 同时核对无误。
- **μB=0 零漂移**：精确电荷共轭对称下两通道速率相等，实测 6 条记录漂移
  `< 1e-15`。
- **μB≠0 全体负向漂移**：C 对称破缺叠加质量差使 ŵ(sū) > ŵ(us̄)，
  `omega_s` 增大 → `tau_s` 减小 → 输运系数负向移动；对 s 味最敏感的
  `kappa_ss` 跌幅最大。
- **旁证**：`ward_moller_sphere_20261121_v1/rate_assembly_drift.json`
  显示 μB=0 漂移 1e-10、μB=450 的 `tau_s` 漂移 9.68%，签名一致。

### v2 → v3 逐点差异（`expected` 列，`refined` 列同量级）

| record | v2 | v3 | 相对变化 |
| --- | --- | --- | --- |
| kappa_mu0_t210_bb | 0.05057040 | 0.05057040 | 0.0000% |
| kappa_mu0_t210_qq | 0.10627615 | 0.10627615 | 0.0000% |
| kappa_mu0_t210_ss | 0.12089917 | 0.12089917 | 0.0000% |
| kappa_mu0_t250_bb | 0.06545485 | 0.06545485 | 0.0000% |
| kappa_mu0_t250_qq | 0.13673962 | 0.13673962 | 0.0000% |
| kappa_mu0_t250_ss | 0.16138505 | 0.16138505 | 0.0000% |
| kappa_mu100_t210_bb | 0.05048711 | 0.05025033 | −0.4690% |
| kappa_mu100_t210_qq | 0.10975413 | 0.10944793 | −0.2790% |
| kappa_mu100_t210_ss | 0.12458691 | 0.12204940 | −2.0367% |
| kappa_mu100_t250_bb | 0.06587922 | 0.06572164 | −0.2392% |
| kappa_mu100_t250_qq | 0.14080159 | 0.14061536 | −0.1323% |
| kappa_mu100_t250_ss | 0.16559190 | 0.16400456 | −0.9586% |
| lambda_mu100_t210 | 775.90228 | 772.26336 | −0.4690% |
| lambda_mu100_t250 | 1459.25072 | 1455.76021 | −0.2392% |

### 生成与复现

候选由 `scripts/dev/export_kappa_lambda_thermal_weighted_baseline.jl` 生成，
目标文件必须不存在；晋升由
`scripts/dev/promote_kappa_lambda_thermal_weighted_baseline.py` 完成，写入
`baseline_sha256` 与源文件指纹。逐点比对使用
`scripts/dev/compare_kappa_lambda_baseline_candidate.py`。

- 主机：quote-remote（`W`），Julia 1.12.5，8 线程；
- 复核：候选的 `expected` 列与 CI 报错值逐位一致（`0.050250331752419587`
  等四例全等），本机 Julia 1.12.5 复算结果相同；
- 回归：`test_kappa_lambda_thermal_weighted_regression.jl` **74/74 通过**（8m36s）；
- v2 保留为历史，测试消费者已迁到 v3。

