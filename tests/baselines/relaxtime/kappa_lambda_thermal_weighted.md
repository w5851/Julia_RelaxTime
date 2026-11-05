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
